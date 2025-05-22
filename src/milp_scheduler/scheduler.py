BIG_M = 100000  # Big-M constant for conditional constraints
"""
MILP-based job scheduler for Snakemake.

This module provides an enhanced job scheduler that uses mixed-integer linear programming
to optimize job execution across heterogeneous computing resources.
"""

import os
import sys
import json
import time
import math
import logging
import hashlib
from typing import Dict, List, Set, Tuple, Optional, Any, Union

# Optional imports with graceful fallbacks
try:
    import pulp
    PULP_AVAILABLE = True
except ImportError:
    PULP_AVAILABLE = False

try:
    import networkx as nx
    NETWORKX_AVAILABLE = True
except ImportError:
    NETWORKX_AVAILABLE = False

from snakemake.jobs import Job
from snakemake.scheduler import JobScheduler

logger = logging.getLogger("snakemake.scheduler")

class MILPJobScheduler(JobScheduler):
    """Enhanced job scheduler using MILP optimization for heterogeneous resources."""

    def __init__(self, workflow, executor_plugin):
        """Initialize the MILP job scheduler.
        
        Args:
            workflow: The Snakemake workflow
            executor_plugin: The executor plugin to use
        """
        # Initialize the base scheduler first
        super().__init__(workflow, executor_plugin)
        
        # Initialize MILP-specific attributes
        self.node_assignments = {}  # Stores job-to-node assignments
        self._job_start_times = {}  # For tracking job execution times
        self._finished_jobs_history = set()  # Track recorded job history
        
        # Critical path tracking (for extended MILP)
        self.dag_graph = None
        self.critical_path = None
        self.critical_path_length = 0
        
        # Override job selector based on scheduler type
        if self.workflow.scheduling_settings.scheduler == "milp":
            if not PULP_AVAILABLE:
                logger.warning(
                    "Falling back to greedy scheduler because PuLP package is not available. "
                    "Install with 'pip install pulp'."
                )
            else:
                self._job_selector = self.job_selector_milp
                
        elif self.workflow.scheduling_settings.scheduler == "milp-ext":
            if not PULP_AVAILABLE:
                logger.warning(
                    "Falling back to greedy scheduler because PuLP package is not available. "
                    "Install with 'pip install pulp'."
                )
            elif not NETWORKX_AVAILABLE:
                logger.warning(
                    "Falling back to regular MILP scheduler because NetworkX is not available. "
                    "Install with 'pip install networkx'."
                )
                self._job_selector = self.job_selector_milp
            else:
                self._job_selector = self.job_selector_milp_ext
    
    def run(self, jobs, executor=None):
        """Run jobs with tracking of execution times.
        
        Args:
            jobs: List of jobs to run
            executor: Executor to use (optional)
        """
        # Record start times for all jobs
        current_time = time.time()
        for job in jobs:
            self._job_start_times[job.jobid] = current_time

        # Call the original run method
        super().run(jobs, executor)
    
    def _proceed(self, job):
        """Handle job completion with history recording.
        
        Args:
            job: The completed job
        """
        with self._lock:
            self._tofinish.append(job)

            # Record job history if job was previously started and not yet recorded
            if job.jobid in self._job_start_times and job.jobid not in self._finished_jobs_history:
                try:
                    # Calculate runtime
                    start_time = self._job_start_times[job.jobid]
                    end_time = time.time()
                    runtime = end_time - start_time

                    # Get node assignment if available
                    node_assignment = self.get_node_assignment(job) or {
                        "cluster": "unknown",
                        "node": "unknown",
                        "start_time": start_time,
                        "end_time": end_time
                    }

                    # Add actual runtime
                    node_assignment["actual_runtime"] = runtime

                    # Record in job history
                    self._record_job_history(job, node_assignment)

                    # Mark as recorded
                    self._finished_jobs_history.add(job.jobid)

                    # Clean up
                    del self._job_start_times[job.jobid]
                except Exception as e:
                    logger.debug(f"Failed to record job history: {e}")

            # Original _proceed code from parent class
            if self.dryrun:
                if len(self.running) - len(self._tofinish) - len(self._toerror) <= 0:
                    self._open_jobs.release()
            else:
                self._open_jobs.release()
    
    def job_selector_milp(self, jobs):
        """MILP-based job selector with resource optimization.
        
        Args:
            jobs: Set of jobs to select from
            
        Returns:
            Set of selected jobs
        """
        if len(jobs) == 1:
            # For single job, just return it directly
            return set(jobs)

        with self._lock:
            if not self.resources["_cores"]:
                if self.resources.get('_cores', 0) <= 0:
                    raise ValueError('Total cores must be > 0 for MILP scheduler')
                return set()

            try:
                # Convert jobs set to list for indexing
                jobs_list = list(jobs)

                # 1. Load scheduler configuration
                config = self._load_scheduler_config()
        logger.debug(f'Loaded scheduler config: {config}')

                # 2. Load system resource profile
                system_profile = self._load_system_profile(config)
        logger.debug(f'Loaded system profile: {list(system_profile.keys())}')

                # 3. Load historical job performance data
                historical_data = self._load_historical_data(jobs_list, config) if \
                    config["scheduler"]["estimation"]["history"]["enabled"] else {}

                # 4. Extract available nodes from system profile
                available_nodes = self._extract_available_nodes(system_profile)

                # 5. Check for sufficient nodes
                if not available_nodes:
                    logger.warning("No nodes available in system profile. Falling back to greedy scheduler.")
                    return self.job_selector_greedy(jobs)

                # 6. Log available resources
                self._log_available_resources(available_nodes)

                # 7. Process job requirements and estimate runtimes
                job_requirements = {}
                for job in jobs_list:
                    job_requirements[job.jobid] = self._extract_job_requirements(job, config)
                    self._estimate_job_io_sizes(job_requirements[job.jobid], job, config)
                    self._estimate_job_runtime(job_requirements[job.jobid], historical_data, config)

                # 8. Create MILP problem
                prob = self._create_milp_problem(jobs_list, job_requirements, available_nodes, config)

                # 9. Solve the MILP problem
                solution = self._solve_milp_problem(prob, config)

                # 10. Check solution status
                if solution["status"] != "Optimal":
                    logger.warning(
                        f"Failed to find optimal solution with MILP scheduler ({solution['status']}). "
                        f"Falling back to {config['scheduler']['optimization']['fallback']} scheduler."
                    )
                    if config["scheduler"]["optimization"]["fallback"] == "greedy":
                        return self.job_selector_greedy(jobs)
                    else:
                        return self.job_selector_ilp(jobs)

                # 11. Process solution and extract job assignments
                selected_jobs, node_assignments = self._process_milp_solution(
                    solution, jobs_list, job_requirements, available_nodes
                )

                # 12. Store node assignments for executor
                self.node_assignments = node_assignments

                # 13. Update resources
                self.update_available_resources(selected_jobs)

                # 14. Return selected jobs
                return selected_jobs

            except Exception as e:
                logger.warning(f"Error in MILP scheduler: {str(e)}")
                logger.warning("Falling back to greedy scheduler.")
                return self.job_selector_greedy(jobs)
    
    def job_selector_milp_ext(self, jobs):
        """Enhanced MILP scheduler with critical path optimization.
        
        Args:
            jobs: Set of jobs to select from
            
        Returns:
            Set of selected jobs
        """
        logger.info("Using enhanced MILP scheduler with critical path analysis")

        # Check if we need to use fallback from previous attempt
        if hasattr(self, '_use_fallback') and self._use_fallback and hasattr(self, '_fallback_jobs'):
            fallback_jobs = self._fallback_jobs
            fallback_type = getattr(self, '_fallback_type', 'greedy')

            # Reset fallback flags
            self._use_fallback = False
            self._fallback_jobs = None
            self._fallback_type = None

            # Use appropriate fallback scheduler
            if fallback_type == 'ilp':
                logger.info("Using ILP fallback scheduler after previous MILP failure")
                return self.job_selector_ilp(fallback_jobs)
            else:
                logger.info("Using greedy fallback scheduler after previous MILP failure")
                return self.job_selector_greedy(fallback_jobs)

        if len(jobs) == 1:
            # Calculate critical path even for single job
            all_jobs = list(self.workflow.dag.needrun_jobs())
            logger.info(f"Calculating critical path for all {len(all_jobs)} jobs")

            dag_graph = self._build_dag_graph_ext(all_jobs)
            if dag_graph:
                # Calculate and log critical path
                critical_path, critical_path_length = self._calculate_critical_path_ext(dag_graph)

                # Single job is already implicitly optimal, so report its runtime as makespan
                job = list(jobs)[0]
                logger.info(f"Selected single job with makespan: {job.resources.get('runtime', 1)}")

            return set(jobs)

        with self._lock:
            if not self.resources["_cores"]:
                if self.resources.get('_cores', 0) <= 0:
                    raise ValueError('Total cores must be > 0 for MILP scheduler')
                return set()

            try:
                # Calculate critical path with ALL jobs, not just ready jobs
                all_jobs = list(self.workflow.dag.needrun_jobs())
                logger.info(f"Calculating critical path for all {len(all_jobs)} jobs")

                # Build DAG graph with optimized runtimes
                dag_graph = self._build_dag_graph_ext(all_jobs)

                # Calculate critical path
                critical_path = None
                critical_path_length = 0

                if dag_graph:
                    # Calculate the critical path
                    critical_path, critical_path_length = self._calculate_critical_path_ext(dag_graph)

                if critical_path_length <= 0:
                    # If we couldn't calculate the critical path, use a reasonable default
                    logger.warning("Could not calculate critical path length, using default")
                    critical_path_length = max([job.resources.get("runtime", 1) for job in jobs], default=1)

                # Convert jobs set to list for indexing
                jobs_list = list(jobs)

                # 1. Load scheduler configuration
                config = self._load_scheduler_config()
        logger.debug(f'Loaded scheduler config: {config}')

                # 2. Load system resource profile
                system_profile = self._load_system_profile(config)
        logger.debug(f'Loaded system profile: {list(system_profile.keys())}')

                # 3. Load historical job performance data
                historical_data = self._load_historical_data(jobs_list, config) if \
                    config["scheduler"]["estimation"]["history"]["enabled"] else {}

                # 4. Extract available nodes from system profile
                available_nodes = self._extract_available_nodes(system_profile)

                # 5. Check for sufficient nodes
                if not available_nodes:
                    logger.warning("No nodes available in system profile. Falling back to greedy scheduler.")
                    return self.job_selector_greedy(jobs)

                # 6. Process job requirements and estimate runtimes
                job_requirements = {}
                for job in jobs_list:
                    job_requirements[job.jobid] = self._extract_job_requirements(job, config)
                    self._estimate_job_io_sizes(job_requirements[job.jobid], job, config)
                    self._estimate_job_runtime(job_requirements[job.jobid], historical_data, config)

                # 7. Create the MILP problem
                import pulp
                prob = pulp.LpProblem("EnhancedScheduler", pulp.LpMinimize)

                # 8. Create variables
                # Job assignment variables
                x = {}
                for job in jobs_list:
                    job_id = job.jobid
                    x[job_id] = {}
                    for node in available_nodes:
                        x[job_id][node["name"]] = pulp.LpVariable(f"job_{job_id}_node_{node['name']}", cat="Binary")

                # Time variables
                start_time = {job.jobid: pulp.LpVariable(f"start_{job.jobid}", lowBound=0, cat="Continuous")
                              for job in jobs_list}

                end_time = {job.jobid: pulp.LpVariable(f"end_{job.jobid}", lowBound=0, cat="Continuous")
                            for job in jobs_list}

                # Makespan variable
                makespan = pulp.LpVariable("makespan", lowBound=0, cat="Continuous")

                # 9. Define objective: minimize makespan
                prob += makespan

                # 10. Constraints

                # Each job must be assigned to exactly one node
                for job in jobs_list:
                    prob += pulp.lpSum(
                        [x[job.jobid][node["name"]] for node in available_nodes]) == 1, f"Job_{job.jobid}_assignment"

                # Node capacity constraints
                for node in available_nodes:
                    node_data = node["data"]

                    # CPU cores constraint
                    prob += pulp.lpSum([
                        x[job.jobid][node["name"]] * job_requirements[job.jobid]["cores"]
                        for job in jobs_list
                    ]) <= node_data["resources"]["cores"], f"Node_{node['name']}_cores"

                    # Memory constraint
                    prob += pulp.lpSum([
                        x[job.jobid][node["name"]] * job_requirements[job.jobid]["memory_mb"]
                        for job in jobs_list
                    ]) <= node_data["resources"]["memory_mb"], f"Node_{node['name']}_memory"

                # Feature compatibility constraints
                for job in jobs_list:
                    job_id = job.jobid
                    job_features = job_requirements[job_id]["features"]

                    for node in available_nodes:
                        node_features = node["data"].get("features", [])

                        # If job requires features not present on node, prevent assignment
                        if job_features and not all(feature in node_features for feature in job_features):
                            prob += x[job_id][node["name"]] == 0, f"Job_{job_id}_incompatible_with_{node['name']}"

                # Runtime calculation and constraints
                for job in jobs_list:
                    job_id = job.jobid

                    for node in available_nodes:
                        node_name = node["name"]
                        node_data = node["data"]

                        # Calculate runtime on this node
                        runtime = self._calculate_runtime_on_node(job_requirements[job_id], node_data)

                        # Using big-M method for conditional constraints
                        M = BIG_M  # Big number

                        # If job assigned to node, end_time >= start_time + runtime
                        prob += end_time[job_id] >= start_time[job_id] + runtime * x[job_id][node_name] - M * (
                                1 - x[job_id][node_name]), f"Runtime_{job_id}_on_{node_name}"

                # Dependency constraints
                for job in jobs_list:
                    # Get job dependencies
                    dependencies = [dep.jobid for dep in self.workflow.dag.dependencies[job]]

                    # For each dependency, ensure this job starts after dependency ends
                    for dep_id in dependencies:
                        if dep_id in end_time:  # Only if the dependency is in our job set
                            prob += start_time[job.jobid] >= end_time[dep_id], f"Dependency_{dep_id}_before_{job.jobid}"

                # Makespan constraint (makespan >= end time of all jobs)
                for job in jobs_list:
                    prob += makespan >= end_time[job.jobid], f"Makespan_for_{job.jobid}"

                # Critical path constraint with a small relaxation factor to avoid infeasibility
                time_limit = config["scheduler"]["optimization"]["time_limit_seconds"]
                base_relaxation = 0.95
                relaxation_factor = base_relaxation

                while True:
                    # Create a fresh problem with relaxed constraint
                    min_makespan = critical_path_length * relaxation_factor
                    prob += makespan >= min_makespan, "Critical_path_constraint"
                    logger.info(
                        f"Trying critical path constraint: makespan >= {min_makespan} (relaxation: {relaxation_factor})")

                    # Try to solve
                    solver = pulp.PULP_CBC_CMD(msg=False, timeLimit=time_limit)
                    prob.solve(solver)

                    # If solution found or no more relaxation possible, break
                    if pulp.LpStatus[prob.status] == "Optimal" or relaxation_factor < 0.5:
                        break

                    # Relax further and retry
                    relaxation_factor -= 0.1
                    prob.constraints.pop("Critical_path_constraint")
                
                logger.info(
                    f"Added critical path constraint: makespan >= {min_makespan} (relaxed from {critical_path_length})")

                # 11. Solve the problem
                solver = pulp.PULP_CBC_CMD(msg=False, timeLimit=time_limit)
                prob.solve(solver)

                # 12. Check solution status
                if pulp.LpStatus[prob.status] != "Optimal":
                    logger.warning(f"Failed to find optimal solution. Status: {pulp.LpStatus[prob.status]}")

                    # Set a flag and save the jobs for fallback processing
                    self._use_fallback = True
                    self._fallback_jobs = jobs
                    self._fallback_type = config["scheduler"]["optimization"]["fallback"]

                    # Return empty set from MILP scheduler
                    # The next scheduler call will use the fallback method
                    return set()

                # 13. Process solution
                selected_jobs = set()
                node_assignments = {}

                for job in jobs_list:
                    for node in available_nodes:
                        node_name = node["name"]
                        var_name = f"job_{job.jobid}_node_{node_name}"

                        # Find the variable
                        var = next((v for v in prob.variables() if v.name == var_name), None)

                        if var and pulp.value(var) > 0.5:
                            selected_jobs.add(job)

                            # Find start and end times
                            start_var = next((v for v in prob.variables() if v.name == f"start_{job.jobid}"), None)
                            end_var = next((v for v in prob.variables() if v.name == f"end_{job.jobid}"), None)

                            start_time_val = pulp.value(start_var) if start_var else 0
                            end_time_val = pulp.value(end_var) if end_var else 0

                            # Store assignment
                            node_assignments[job.jobid] = {
                                "cluster": node["cluster"],
                                "node": node_name,
                                "start_time": start_time_val,
                                "end_time": end_time_val
                            }

                # 14. Store node assignments for executor
                self.node_assignments = node_assignments

                # Log the makespan
                makespan_value = pulp.value(makespan)
                logger.info(f"MILP solution found with makespan: {makespan_value}")
                logger.info(f"Selected {len(selected_jobs)} jobs")

                # 15. Return selected jobs
                if not selected_jobs:
                    logger.warning("No jobs selected. Falling back to greedy scheduler.")
                    return self.job_selector_greedy(jobs)

                # Update available resources
                self.update_available_resources(selected_jobs)
                
                return selected_jobs

            except Exception as e:
                logger.warning(f"Error in enhanced MILP scheduler: {str(e)}")
                logger.warning("Falling back to greedy scheduler.")
                return self.job_selector_greedy(jobs)

    def _load_scheduler_config(self):
        """Load scheduler configuration settings with improved path resolution."""
        # Try multiple locations in order
        possible_paths = [
            # 1. User-specified path (if any)
            getattr(self.workflow, "scheduler_config_path", None),

            # 2. Current working directory
            os.path.join(os.getcwd(), "scheduler_config.yaml"),

            # 3. Same directory as Snakefile
            os.path.join(os.path.dirname(getattr(self.workflow, "snakefile", ".")), "scheduler_config.yaml"),

            # 4. Config directory in Snakemake package
            os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                         "config", "scheduler_config.yaml"),

            # 5. User's Snakemake config directory
            os.path.expanduser("~/.snakemake/scheduler_config.yaml")
        ]

        import yaml
        # Try each path
        for path in possible_paths:
            if path and os.path.exists(path):
                try:
                    with open(path, 'r') as f:
                        config = yaml.safe_load(f)
                    logger.debug(f"Loaded scheduler config from {path}")
                    return config
                except Exception as e:
                    logger.warning(f"Failed to load scheduler config from {path}: {e}")

        # If no file found, use default config
        logger.warning("Using default configuration (no config file found)")
        return self._get_default_config()

    def _load_system_profile(self, config):
        """Load system resource profile with improved path resolution."""
        profile_path = config["scheduler"]["paths"]["system_profile"]

        # Try multiple locations in order
        possible_paths = [
            # 1. Absolute path
            profile_path if os.path.isabs(profile_path) else None,

            # 2. Relative to working directory
            os.path.join(os.getcwd(), profile_path),

            # 3. Relative to Snakefile directory
            os.path.join(os.path.dirname(getattr(self.workflow, "snakefile", ".")), profile_path),

            # 4. User's Snakemake config directory
            os.path.expanduser(f"~/.snakemake/{profile_path}")
        ]

        # Try each path
        for path in [p for p in possible_paths if p]:
            if os.path.exists(path):
                try:
                    with open(path, 'r') as f:
                        profile = json.load(f)
                    logger.debug(f"Loaded system profile from {path}")
                    return profile
                except Exception as e:
                    logger.warning(f"Failed to load system profile from {path}: {e}")

        # Default profile with local resources if file doesn't exist
        logger.warning(f"System profile not found. Using default profile.")
        return self._get_default_system_profile()

    def _get_default_system_profile(self):
        """Create a default system profile based on local resources."""
        return {
            "clusters": {
                "local": {
                    "nodes": {
                        "default": {
                            "resources": {
                                "cores": self.resources["_cores"],
                                "memory_mb": self.resources.get("mem_mb", 1000),
                                "gpu_count": 0,
                                "local_storage_mb": 10000
                            },
                            "features": ["cpu", "x86_64"],
                            "properties": {
                                "cpu_flops": 100000000000,
                                "memory_bandwidth_mbps": 10000,
                                "read_mbps": 1000,
                                "write_mbps": 800
                            }
                        }
                    }
                }
            }
        }

    def _extract_available_nodes(self, system_profile):
        """Extract available nodes from system profile."""
        available_nodes = []

        for cluster_name, cluster in system_profile.get("clusters", {}).items():
            for node_name, node_data in cluster.get("nodes", {}).items():
                available_nodes.append({
                    "cluster": cluster_name,
                    "name": node_name,
                    "data": node_data
                })

        return available_nodes

    def _load_historical_data(self, jobs, config):
        """Load historical performance data for jobs."""
        historical_data = {}
        history_dir = config["scheduler"]["paths"]["job_history"]

        if not os.path.exists(history_dir):
            os.makedirs(history_dir, exist_ok=True)
            return historical_data

        # Try to load history for each job
        for job in jobs:
            job_hash = self._get_job_hash(job)
            history_file = os.path.join(history_dir, f"{job_hash}.json")

            if os.path.exists(history_file):
                try:
                    with open(history_file, 'r') as f:
                        job_history = json.load(f)
                    historical_data[job.jobid] = job_history
                except Exception as e:
                    logger.debug(f"Failed to load history for job {job.jobid}: {e}")

        return historical_data

    def _get_job_hash(self, job):
        """Generate a unique hash for a job based on rule and params."""
        import hashlib

        # Combine rule name and input/output patterns
        rule_name = job.rule.name
        input_patterns = str(sorted([str(i) for i in job.input]))
        output_patterns = str(sorted([str(o) for o in job.output]))

        # Include params hash if available
        params_hash = ""
        if hasattr(job.rule, "params"):
            params_str = str(job.rule.params)
            params_hash = hashlib.md5(params_str.encode()).hexdigest()

        # Create combined hash
        combined = f"{rule_name}_{input_patterns}_{output_patterns}_{params_hash}"
        return hashlib.md5(combined.encode()).hexdigest()

    def _extract_job_requirements(self, job, config):
        """Extract job requirements from Snakefile job specification."""
        # Start with standard resources
        requirements = {
            "jobid": job.jobid,
            "rule_name": job.rule.name,
            "cores": job.resources.get("_cores", 1),
            "memory_mb": job.resources.get("mem_mb", 1000),  # 1 GiB default
            "runtime_minutes": job.resources.get("runtime", 30),  # 30 minutes default
            "disk_mb": job.resources.get("disk_mb", 0),

            # Initialize extended specification fields
            "features": [],
            "resources": {},
            "properties": {}
        }

        # Extract extended specification if available
        if hasattr(job.rule, "params") and hasattr(job.rule.params, "job_specification"):
            job_spec = job.rule.params.job_specification

            # Extract features
            if "features" in job_spec:
                requirements["features"] = job_spec["features"]

            # Extract additional resources
            if "resources" in job_spec:
                requirements["resources"].update(job_spec["resources"])

            # Extract properties
            if "properties" in job_spec:
                requirements["properties"].update(job_spec["properties"])

            # Validate job specification against nomenclature
            valid, warnings = self.validate_job_specification(job_spec, config)
            if warnings:
                logger.debug(f"Job {job.jobid} ({job.rule.name}) has {len(warnings)} nomenclature warnings")

        return requirements

    def _estimate_job_io_sizes(self, job_requirements, job, config):
        """Estimate I/O sizes for job if not explicitly specified."""
        # Check if auto-estimation is enabled
        if not config["scheduler"]["estimation"]["auto_estimate_file_sizes"]:
            return

        # If sizes are already specified, don't override
        if "input_size_mb" in job_requirements["resources"] and "output_size_mb" in job_requirements["resources"]:
            return

        # Estimate input size from actual files
        if "input_size_mb" not in job_requirements["resources"]:
            input_size_bytes = 0
            for input_file in job.input:
                try:
                    input_path = str(input_file)
                    if os.path.exists(input_path):
                        input_size_bytes += os.path.getsize(input_path)
                except Exception as e:
                    logger.debug(f"Error estimating size of {input_file}: {e}")

            # Convert to MB
            job_requirements["resources"]["input_size_mb"] = max(1, input_size_bytes / (1024 * 1024))

        # For output, use historical data or simple heuristic if available
        if "output_size_mb" not in job_requirements["resources"]:
            # Default to same as input or at least 1MB
            job_requirements["resources"]["output_size_mb"] = max(1,
                                                                 job_requirements["resources"].get("input_size_mb", 1))

    def _estimate_job_runtime(self, job_requirements, historical_data, config):
        """Estimate job runtime based on requirements, system properties, and history."""
        # Start with baseline runtime from resources
        baseline_runtime = job_requirements["runtime_minutes"]

        # Apply historical correction if available
        job_id = job_requirements["jobid"]
        if job_id in historical_data and "statistics" in historical_data[job_id]:
            hist_runtime = historical_data[job_id]["statistics"].get("avg_runtime_minutes", 0)
            if hist_runtime > 0:
                # Blend historical with baseline
                weight = config["scheduler"]["estimation"]["history"]["adaptation_weight"]
                job_requirements["estimated_runtime"] = weight * hist_runtime + (1 - weight) * baseline_runtime
        else:
            # Apply uncertainty factor
            uncertainty = config["scheduler"]["runtime"]["uncertainty_factor"]
            job_requirements["estimated_runtime"] = baseline_runtime * (1 + uncertainty)

    def _create_milp_problem(self, jobs, job_requirements, available_nodes, config):
        """Create MILP problem for job scheduling."""
        import pulp

        # Create problem
        prob = pulp.LpProblem("HeterogeneousScheduler", pulp.LpMinimize)

        # Create variables
        # Job assignment variables (1 if job j assigned to node n, 0 otherwise)
        x = {}
        for j in range(len(jobs)):
            job = jobs[j]
            job_id = job.jobid
            x[job_id] = {}
            for n in range(len(available_nodes)):
                node = available_nodes[n]
                x[job_id][node["name"]] = pulp.LpVariable(f"job_{job_id}_node_{node['name']}", cat="Binary")

        # Time variables
        start_time = {job.jobid: pulp.LpVariable(f"start_{job.jobid}", lowBound=0, cat="Continuous")
                     for job in jobs}

        end_time = {job.jobid: pulp.LpVariable(f"end_{job.jobid}", lowBound=0, cat="Continuous")
                   for job in jobs}

        # Makespan variable
        makespan = pulp.LpVariable("makespan", lowBound=0, cat="Continuous")

        # Energy variable (if energy optimization enabled)
        if config["scheduler"]["optimization"]["objective_weights"]["energy"] > 0:
            total_energy = pulp.LpVariable("total_energy", lowBound=0, cat="Continuous")

            # Calculate maximum possible energy for normalization
            max_possible_energy = sum(
                job_requirements[job.jobid].get("resources", {}).get("max_watts", 100) *
                job_requirements[job.jobid]["estimated_runtime"]
                for job in jobs
            )
        else:
            total_energy = None
            max_possible_energy = 1

        # Objective function
        makespan_weight = config["scheduler"]["optimization"]["objective_weights"]["makespan"]
        energy_weight = config["scheduler"]["optimization"]["objective_weights"]["energy"]

        # Maximum possible makespan for normalization
        max_possible_makespan = sum(job_requirements[job.jobid]["estimated_runtime"] for job in jobs)

        # Define objective
        if total_energy:
            prob += (makespan_weight * (makespan / max_possible_makespan) +
                    energy_weight * (total_energy / max_possible_energy))
        else:
            prob += makespan

        # Constraints

        # 1. Each job must be assigned to exactly one node
        for job in jobs:
            prob += pulp.lpSum(
                [x[job.jobid][node["name"]] for node in available_nodes]) == 1, f"Job_{job.jobid}_assignment"

        # 2. Node capacity constraints
        for n, node in enumerate(available_nodes):
            node_data = node["data"]

            # CPU cores constraint
            prob += pulp.lpSum([
                x[job.jobid][node["name"]] * job_requirements[job.jobid]["cores"]
                for job in jobs
            ]) <= node_data["resources"]["cores"], f"Node_{node['name']}_cores"

            # Memory constraint
            prob += pulp.lpSum([
                x[job.jobid][node["name"]] * job_requirements[job.jobid]["memory_mb"]
                for job in jobs
            ]) <= node_data["resources"]["memory_mb"], f"Node_{node['name']}_memory"

            # GPU memory constraint (if used by any job)
            if any("gpu_memory_mb" in job_requirements[job.jobid].get("resources", {}) for job in jobs):
                gpu_memory = node_data["resources"].get("gpu_memory_mb", 0)
                if gpu_memory > 0:
                    prob += pulp.lpSum([
                        x[job.jobid][node["name"]] * job_requirements[job.jobid].get("resources", {}).get(
                            "gpu_memory_mb", 0)
                        for job in jobs
                    ]) <= gpu_memory, f"Node_{node['name']}_gpu_memory"

        # 3. Feature compatibility constraints
        for j, job in enumerate(jobs):
            job_id = job.jobid
            job_features = job_requirements[job_id]["features"]

            for n, node in enumerate(available_nodes):
                node_features = node["data"].get("features", [])

                # If job requires features not present on node, prevent assignment
                if not all(feature in node_features for feature in job_features):
                    prob += x[job_id][node["name"]] == 0, f"Job_{job_id}_incompatible_with_{node['name']}"

        # 4. Runtime calculation and constraints
        for j, job in enumerate(jobs):
            job_id = job.jobid

            for n, node in enumerate(available_nodes):
                node_name = node["name"]
                node_data = node["data"]

                # Calculate runtime on this node
                runtime = self._calculate_runtime_on_node(job_requirements[job_id], node_data)

                # Link start and end times with assignment
                # Using big-M method: If job j is assigned to node n, enforce runtime constraint
                M = BIG_M  # Big number

                # If job assigned to node, end_time >= start_time + runtime
                prob += end_time[job_id] >= start_time[job_id] + runtime * x[job_id][node_name] - M * (
                        1 - x[job_id][node_name])

                # If job not assigned to node, constraint is relaxed
                prob += end_time[job_id] >= start_time[job_id]

        # 5. Dependency constraints
        for j, job in enumerate(jobs):
            # Get job dependencies
            dependencies = [dep.jobid for dep in self.workflow.dag.dependencies[job]]

            # For each dependency, ensure this job starts after dependency ends
            for dep_id in dependencies:
                if dep_id in end_time:  # Only if the dependency is in our job set
                    prob += start_time[job.jobid] >= end_time[dep_id]

        # 6. Makespan constraint (makespan >= end time of all jobs)
        for job in jobs:
            prob += makespan >= end_time[job.jobid]

        # 7. Energy constraints (if enabled)
        if total_energy:
            energy_expressions = []
            for j, job in enumerate(jobs):
                job_id = job.jobid
                job_req = job_requirements[job_id]

                for n, node in enumerate(available_nodes):
                    node_name = node["name"]
                    node_data = node["data"]

                    # Calculate power for this job on this node
                    power = job_req.get("resources", {}).get("max_watts", 100)

                    # Calculate runtime for this job on this node
                    runtime = self._calculate_runtime_on_node(job_req, node_data)

                    # Energy = Power * Time * Assignment
                    energy_expressions.append(power * runtime * x[job_id][node_name])

            # Total energy is sum of all job energies
            prob += total_energy == pulp.lpSum(energy_expressions)

        return prob

    def _calculate_runtime_on_node(self, job_req, node_data):
        """Calculate the runtime of a job on a specific node."""
        # Start with estimated runtime
        runtime = job_req.get("estimated_runtime", job_req["runtime_minutes"])

        # If job specifies computational requirements and node has matching properties
        if "properties" in job_req and "properties" in node_data:
            job_props = job_req["properties"]
            node_props = node_data["properties"]

            # Skip FLOPS-based calculation and use the specified runtime
            # This ensures we use the actual runtime from the Snakefile

            # Only calculate I/O time to add to the base runtime
            io_time = 0
            if "input_size_mb" in job_req.get("resources", {}):
                input_size = job_req["resources"]["input_size_mb"]

                # Find appropriate bandwidth
                read_bandwidth = node_props.get("read_mbps", 0)
                if read_bandwidth > 0:
                    # Calculate I/O time but scale it down to be a small fraction of runtime
                    io_time += min(1.0, input_size / read_bandwidth * 0.1)  # Scale factor 0.1

            if "output_size_mb" in job_req.get("resources", {}):
                output_size = job_req["resources"]["output_size_mb"]

                # Find appropriate bandwidth
                write_bandwidth = node_props.get("write_mbps", 0)
                if write_bandwidth > 0:
                    # Calculate I/O time but scale it down to be a small fraction of runtime
                    io_time += min(1.0, output_size / write_bandwidth * 0.1)  # Scale factor 0.1

            # Total time = base runtime + small I/O adjustment
            calculated_runtime = runtime + io_time

            # Use calculated runtime
            return calculated_runtime

        # Otherwise use the estimated runtime
        return runtime

    def _solve_milp_problem(self, prob, config):
        """Solve the MILP problem with time limit."""
        import pulp

        # Get time limit from config
        time_limit = config["scheduler"]["optimization"]["time_limit_seconds"]

        # Setup solver with time limit
        solver = pulp.PULP_CBC_CMD(msg=False, timeLimit=time_limit)

        # Solve the problem
        prob.solve(solver)

        # Return solution info
        return {
            "status": pulp.LpStatus[prob.status],
            "prob": prob,
            "makespan": pulp.value(prob.objective)
        }

    def _process_milp_solution(self, solution, jobs, job_requirements, available_nodes):
        """Process MILP solution and extract job assignments."""
        import pulp

        prob = solution["prob"]
        selected_jobs = set()
        node_assignments = {}

        # Extract variables from problem
        variables = prob.variables()

        # Get job assignment variables
        job_vars = {v.name: v for v in variables if v.name.startswith("job_")}

        # Extract assignments
        for j, job in enumerate(jobs):
            job_id = job.jobid

            for node in available_nodes:
                node_name = node["name"]
                var_name = f"job_{job_id}_node_{node_name}"

                if var_name in job_vars and pulp.value(job_vars[var_name]) > 0.5:
                    selected_jobs.add(job)

                    # Get start and end times
                    start_var = next((v for v in variables if v.name == f"start_{job_id}"), None)
                    end_var = next((v for v in variables if v.name == f"end_{job_id}"), None)

                    start_time = pulp.value(start_var) if start_var else 0
                    end_time = pulp.value(end_var) if end_var else 0

                    # Store assignment
                    node_assignments[job_id] = {
                        "cluster": node["cluster"],
                        "node": node_name,
                        "start_time": start_time,
                        "end_time": end_time
                    }

                    break  # Job can only be assigned to one node

        # Log assignments
        logger.debug(f"MILP solution found with makespan: {solution['makespan']}")
        logger.debug(f"Selected {len(selected_jobs)} jobs")

        return selected_jobs, node_assignments

    def _log_available_resources(self, available_nodes):
        """Log available resources for debugging."""
        logger.debug(f"Available resources: {len(available_nodes)} nodes")

        for node in available_nodes:
            node_data = node["data"]
            cores = node_data["resources"].get("cores", 0)
            memory = node_data["resources"].get("memory_mb", 0)
            features = node_data.get("features", [])
            logger.debug(
                f"Node {node['cluster']}/{node['name']}: {cores} cores, {memory / 1024:.1f}GB memory, features: {features}")

    def _get_default_config(self):
        """Get default scheduler configuration."""
        return {
            "scheduler": {
                "type": "milp",
                "paths": {
                    "system_profile": "system_profile.json",
                    "job_history": "~/.snakemake/job_history",
                    "gnnrl_models": "~/.snakemake/gnnrl_models"
                },
                "nomenclature": {
                    "features": {
                        "computing": [
                            "cpu", "x86_64", "arm64", "gpu", "cuda",
                            "opencl", "fpga", "avx512", "avx2"
                        ],
                        "gpu_versions": [
                            "cuda_10", "cuda_11", "cuda_12", "a100",
                            "v100", "rtx_3090", "rtx_4090"
                        ],
                        "storage": [
                            "ssd", "nvme", "hdd", "fast_storage", 
                            "scratch", "shared_fs", "local_fs"
                        ],
                        "network": [
                            "infiniband", "ethernet_10g", "ethernet_100g",
                            "omni_path", "high_bandwidth"
                        ],
                        "memory": ["high_memory", "numa"],
                        "partitions": [
                            "partition_standard", "partition_gpu",
                            "partition_bigmem", "partition_io"
                        ]
                    },
                    "resources": {
                        "compute": ["cores", "gpu_count", "gpu_memory_mb"],
                        "storage": ["memory_mb", "local_storage_mb", "shared_storage_mb"],
                        "data": ["input_size_mb", "output_size_mb", "temp_size_mb"],
                        "power": ["max_watts"]
                    },
                    "properties": {
                        "computational": ["cpu_flops", "gpu_flops", "integer_ops", "vector_ops"],
                        "memory": ["memory_bandwidth_mbps"],
                        "storage": ["read_mbps", "write_mbps", "iops", "io_pattern"],
                        "network": ["network_bandwidth_mbps", "network_latency_ms"],
                        "power": ["idle_watts", "power_efficiency"]
                    }
                },
                "estimation": {
                    "auto_estimate_file_sizes": True,
                    "history": {
                        "enabled": True,
                        "adaptation_weight": 0.7
                    }
                },
                "runtime": {
                    "default_minutes": 30,
                    "uncertainty_factor": 0.2,
                    "max_timeout_factor": 1.5,
                    "min_runtime_minutes": 1
                },
                "optimization": {
                    "objective_weights": {
                        "makespan": 0.8,
                        "energy": 0.1,
                        "utilization": 0.1
                    },
                    "constraints": {
                        "penalty_scale": 50.0
                    },
                    "time_limit_seconds": 10,
                    "fallback": "greedy"
                },
                "milp": {
                    "solver": "CBC",
                    "threads": 4,
                    "gap_tolerance": 0.05,
                    "formulation": {
                        "use_big_m": True,
                        "big_m_value": 100000
                    }
                },
                "updates": {
                    "check_interval_seconds": 300,
                    "dynamic_resources": False
                },
                "validation": {
                    "enforce_nomenclature": False,
                    "warn_on_unknown_terms": True
                }
            }
        }

    def validate_system_profile(self, profile, config):
        """
        Validate system profile against standardized nomenclature.

        Args:
            profile: System profile dictionary
            config: Scheduler configuration with nomenclature definitions

        Returns:
            Tuple of (valid, warnings)
        """
        valid = True
        warnings = []

        # Skip validation if not enabled
        if not config.get("scheduler", {}).get("validation", {}).get("enforce_nomenclature", False):
            return True, []

        # Get nomenclature definitions
        nomenclature = config.get("scheduler", {}).get("nomenclature", {})
        warn_on_unknown = config.get("scheduler", {}).get("validation", {}).get("warn_on_unknown_terms", True)

        # Flatten allowed features
        allowed_features = []
        for category, features in nomenclature.get("features", {}).items():
            allowed_features.extend(features)

        # Flatten allowed resources
        allowed_resources = []
        for category, resources in nomenclature.get("resources", {}).items():
            allowed_resources.extend(resources)

        # Flatten allowed properties
        allowed_properties = []
        for category, properties in nomenclature.get("properties", {}).items():
            allowed_properties.extend(properties)

        # Check each cluster and node
        for cluster_name, cluster in profile.get("clusters", {}).items():
            for node_name, node in cluster.get("nodes", {}).items():
                # Validate features
                for feature in node.get("features", []):
                    if feature not in allowed_features:
                        msg = f"Unknown feature '{feature}' in node {cluster_name}/{node_name}"
                        if warn_on_unknown:
                            warnings.append(msg)
                            logger.warning(msg)
                        else:
                            valid = False
                            logger.error(msg)

                # Validate resources
                for resource in node.get("resources", {}):
                    if resource not in allowed_resources:
                        msg = f"Unknown resource '{resource}' in node {cluster_name}/{node_name}"
                        if warn_on_unknown:
                            warnings.append(msg)
                            logger.warning(msg)
                        else:
                            valid = False
                            logger.error(msg)

                # Validate properties
                for prop in node.get("properties", {}):
                    if prop not in allowed_properties:
                        msg = f"Unknown property '{prop}' in node {cluster_name}/{node_name}"
                        if warn_on_unknown:
                            warnings.append(msg)
                            logger.warning(msg)
                        else:
                            valid = False
                            logger.error(msg)

        return valid, warnings

    def validate_job_specification(self, job_spec, config):
        """
        Validate job specification against standardized nomenclature.

        Args:
            job_spec: Job specification dictionary from Snakefile
            config: Scheduler configuration with nomenclature definitions

        Returns:
            Tuple of (valid, warnings)
        """
        valid = True
        warnings = []

        # Skip validation if not enabled
        if not config.get("scheduler", {}).get("validation", {}).get("enforce_nomenclature", False):
            return True, []

        # Get nomenclature definitions
        nomenclature = config.get("scheduler", {}).get("nomenclature", {})
        warn_on_unknown = config.get("scheduler", {}).get("validation", {}).get("warn_on_unknown_terms", True)

        # Flatten allowed features
        allowed_features = []
        for category, features in nomenclature.get("features", {}).items():
            allowed_features.extend(features)

        # Flatten allowed resources
        allowed_resources = []
        for category, resources in nomenclature.get("resources", {}).items():
            allowed_resources.extend(resources)

        # Flatten allowed properties
        allowed_properties = []
        for category, properties in nomenclature.get("properties", {}).items():
            allowed_properties.extend(properties)

        # Validate features
        for feature in job_spec.get("features", []):
            if feature not in allowed_features:
                msg = f"Unknown feature '{feature}' in job specification"
                if warn_on_unknown:
                    warnings.append(msg)
                    logger.warning(msg)
                else:
                    valid = False
                    logger.error(msg)

        # Validate resources
        for resource in job_spec.get("resources", {}):
            if resource not in allowed_resources:
                msg = f"Unknown resource '{resource}' in job specification"
                if warn_on_unknown:
                    warnings.append(msg)
                    logger.warning(msg)
                else:
                    valid = False
                    logger.error(msg)

        # Validate properties
        for prop in job_spec.get("properties", {}):
            if prop not in allowed_properties:
                msg = f"Unknown property '{prop}' in job specification"
                if warn_on_unknown:
                    warnings.append(msg)
                    logger.warning(msg)
                else:
                    valid = False
                    logger.error(msg)

        return valid, warnings

    def _record_job_history(self, job, node_assignment):
        """Record job execution history to a file."""
        # Create a history directory if it doesn't exist
        history_dir = os.path.expanduser("~/.snakemake/job_history")
        os.makedirs(history_dir, exist_ok=True)

        # Generate a job identifier
        job_hash = self._get_job_hash(job)
        history_file = os.path.join(history_dir, f"{job_hash}.json")

        # Create the history record
        timestamp = time.time()
        history_entry = {
            "job_id": job.jobid,
            "rule_name": job.rule.name,
            "timestamp": timestamp,
            "date": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(timestamp)),
            "resources": {
                "cores": job.resources.get("_cores", 1),
                "mem_mb": job.resources.get("mem_mb", 1000),
                "runtime": job.resources.get("runtime", 0)
            },
            "assignment": node_assignment,
            "actual_runtime": node_assignment.get("actual_runtime", 0),
            "features": [],
            "properties": {}
        }

        # Add job specification if available
        if hasattr(job.rule, "params") and hasattr(job.rule.params, "job_specification"):
            job_spec = job.rule.params.job_specification
            if "features" in job_spec:
                history_entry["features"] = job_spec["features"]
            if "properties" in job_spec:
                history_entry["properties"] = job_spec["properties"]

        # Load existing history if available
        history = []
        if os.path.exists(history_file):
            try:
                with open(history_file, 'r') as f:
                    history = json.load(f)
            except Exception as e:
                logger.debug(f"Error loading job history: {e}")

        # Add new entry
        if not isinstance(history, list):
            history = []
        history.append(history_entry)

        # Calculate statistics
        if len(history) > 0:
            runtimes = [entry.get("actual_runtime", 0) for entry in history if "actual_runtime" in entry]
            if runtimes:
                avg_runtime = sum(runtimes) / len(runtimes)
                history_stats = {
                    "statistics": {
                        "avg_runtime_minutes": avg_runtime / 60,
                        "count": len(history),
                        "last_updated": timestamp
                    }
                }
                # Add statistics to a separate file for quick lookup
                stats_file = os.path.join(history_dir, f"{job_hash}.stats.json")
                try:
                    with open(stats_file, 'w') as f:
                        json.dump(history_stats, f, indent=2)
                except Exception as e:
                    logger.debug(f"Error saving job history stats: {e}")

        # Save updated history
        try:
            with open(history_file, 'w') as f:
                json.dump(history, f, indent=2)
            logger.debug(f"Job history saved to {history_file}")
        except Exception as e:
            logger.debug(f"Error saving job history: {e}")

    def get_node_assignment(self, job):
        """Get node assignment for a job if available."""
        if hasattr(self, "node_assignments") and job.jobid in self.node_assignments:
            return self.node_assignments[job.jobid]
        return None
        
    def _build_dag_graph_ext(self, all_jobs=None):
        """Build a directed graph representing the DAG of all jobs with optimized runtimes."""
        if not NETWORKX_AVAILABLE:
            return None

        if all_jobs is None:
            all_jobs = list(self.workflow.dag.needrun_jobs())

        # Create directed graph
        G = nx.DiGraph()

        # Add all jobs as nodes with optimized runtimes
        for job in all_jobs:
            # For "all" rule or collector rules, use minimal runtime
            if job.rule.name == "all" or (len(job.input) > 0 and len(job.output) == 0):
                runtime = 0.1  # Small but non-zero
            else:
                # Calculate optimal runtime based on node assignments
                runtime = self._get_optimal_runtime(job)

            # Add job to graph
            G.add_node(job.jobid, runtime=runtime, job=job, rule_name=job.rule.name)

        # Add dependencies as edges
        for job in all_jobs:
            for dep_job in self.workflow.dag.dependencies[job]:
                if dep_job.jobid in G.nodes and job.jobid in G.nodes:
                    G.add_edge(dep_job.jobid, job.jobid)

        return G

    def _get_optimal_runtime(self, job):
        """
        Calculate the optimal (minimum) runtime for a job across all possible nodes.
        Considers both specified runtimes and computational characteristics.
        """
        # Get the base runtime from resources
        base_runtime = job.resources.get("runtime", None)

        # If no job specification, return the base runtime or a default
        if (not hasattr(job.rule, "params") or
                not hasattr(job.rule.params, "job_specification")):
            return base_runtime if base_runtime is not None else 1

        # Get the job specification
        job_spec = job.rule.params.job_specification

        # If no special properties or features, return the base runtime
        if (not job_spec or
                ("properties" not in job_spec and "features" not in job_spec)):
            return base_runtime if base_runtime is not None else 1

        # Load system profile to check available nodes
        try:
            config = self._load_scheduler_config()
        logger.debug(f'Loaded scheduler config: {config}')
            system_profile = self._load_system_profile(config)
        logger.debug(f'Loaded system profile: {list(system_profile.keys())}')
            available_nodes = self._extract_available_nodes(system_profile)
        except Exception as e:
            logger.debug(f"Error loading system profile: {e}")
            return base_runtime if base_runtime is not None else 1

        # If no nodes available, return the base runtime
        if not available_nodes:
            return base_runtime if base_runtime is not None else 1

        # Get job requirements for runtime calculation
        job_features = job_spec.get("features", [])

        # Create a job requirements structure for runtime calculation
        job_req = {
            "jobid": job.jobid,
            "rule_name": job.rule.name,
            "runtime_minutes": base_runtime if base_runtime is not None else 1,
            "cores": job.resources.get("_cores", 1),
            "memory_mb": job.resources.get("mem_mb", 1000),
            "features": job_features,
            "properties": job_spec.get("properties", {}),
            "resources": job_spec.get("resources", {})
        }

        # Track the minimum runtime across all compatible nodes
        min_runtime = float('inf')

        # Check each available node for compatibility and calculate runtime
        for node in available_nodes:
            node_data = node["data"]
            node_features = node_data.get("features", [])

            # Skip if node doesn't have required features
            if not all(feature in node_features for feature in job_features):
                continue

            # Calculate runtime on this node using the full calculation logic
            runtime = self._calculate_runtime_on_node(job_req, node_data)

            # Update minimum runtime
            if runtime < min_runtime:
                min_runtime = runtime

        # If we couldn't find a compatible node, use the base runtime
        if min_runtime == float('inf'):
            return base_runtime if base_runtime is not None else 1

        # Ensure we're using a reasonable value
        return max(0.1, min_runtime)
        
    def _calculate_critical_path_ext(self, dag_graph=None):
        """Calculate critical path using dynamic programming approach."""
        if not NETWORKX_AVAILABLE:
            return None, 0

        if dag_graph is None:
            dag_graph = self._build_dag_graph_ext()

        if dag_graph is None or len(dag_graph.nodes()) == 0:
            return None, 0

        # Store the graph
        self.dag_graph = dag_graph

        try:
            # Get topological sort
            topo_sort = list(nx.topological_sort(dag_graph))

            # Initialize longest path lengths
            longest_path = {node: 0 for node in dag_graph.nodes()}

            # Process nodes in topological order
            for node in topo_sort:
                # Process each successor
                for succ in dag_graph.successors(node):
                    node_runtime = dag_graph.nodes[node]['runtime']
                    if longest_path[succ] < longest_path[node] + node_runtime:
                        longest_path[succ] = longest_path[node] + node_runtime

            # Find sink nodes (no outgoing edges)
            sink_nodes = [n for n in dag_graph.nodes() if dag_graph.out_degree(n) == 0]
            if not sink_nodes:
                sink_nodes = list(dag_graph.nodes())

            # Find maximum path length to any sink node
            max_path_length = 0
            critical_end_node = None

            for node in sink_nodes:
                node_runtime = dag_graph.nodes[node]['runtime']
                path_length = longest_path[node] + node_runtime
                if path_length > max_path_length:
                    max_path_length = path_length
                    critical_end_node = node

            # Reconstruct the critical path
            critical_path = []
            if critical_end_node is not None:
                # Start with the end node
                critical_path = [critical_end_node]

                # Trace backward to find the path
                current = critical_end_node
                while True:
                    pred_candidates = []
                    for pred in dag_graph.predecessors(current):
                        pred_runtime = dag_graph.nodes[pred]['runtime']
                        if abs(longest_path[current] - longest_path[pred] - pred_runtime) < 1e-6:
                            pred_candidates.append(pred)

                    if not pred_candidates:
                        break

                    # Choose predecessor with highest runtime as tiebreaker
                    pred = max(pred_candidates, key=lambda x: dag_graph.nodes[x]['runtime'])
                    critical_path.insert(0, pred)
                    current = pred

            # Record results
            self.critical_path = critical_path
            self.critical_path_length = max_path_length

            # Log critical path info
            if critical_path:
                path_info = []
                for node_id in critical_path:
                    node_data = dag_graph.nodes[node_id]
                    rule_name = node_data.get('rule_name', 'unknown')
                    runtime = node_data['runtime']
                    path_info.append(f"{node_id}({rule_name}, {runtime})")

                # Create path_str by joining path_info elements
                path_str = " → ".join(path_info)
                logger.info(f"Critical path: {path_str}")
                logger.info(f"Critical path length: {max_path_length}")

            return critical_path, max_path_length

        except Exception as e:
            logger.warning(f"Error calculating critical path: {str(e)}")
            # Use maximum runtime as fallback
            try:
                max_runtime = max([dag_graph.nodes[n]['runtime'] for n in dag_graph.nodes()], default=0)
                return None, max_runtime
            except:
                return None, 0