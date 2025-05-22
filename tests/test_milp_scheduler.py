"""Unit tests for the MILP scheduler."""

import os
import sys
import json
import pytest
import tempfile
from unittest.mock import MagicMock, patch

# Check if PuLP is installed 
try:
    import pulp
    PULP_AVAILABLE = True
except ImportError:
    PULP_AVAILABLE = False

# Skip tests if PuLP is not available
pytestmark = pytest.mark.skipif(not PULP_AVAILABLE, reason="PuLP library is required for these tests")

# Import the scheduler
from snakemake.scheduler_plugins.milp_scheduler.scheduler import MILPJobScheduler

# Test fixture for a minimal system profile
@pytest.fixture
def system_profile_file():
    profile = {
        "clusters": {
            "local": {
                "nodes": {
                    "default": {
                        "resources": {
                            "cores": 4,
                            "memory_mb": 8192,
                            "gpu_count": 0
                        },
                        "features": ["cpu", "x86_64"],
                        "properties": {
                            "cpu_flops": 100000000000,
                            "read_mbps": 500,
                            "write_mbps": 400
                        }
                    }
                }
            }
        }
    }
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
        json.dump(profile, f)
        profile_path = f.name
    
    yield profile_path
    
    # Cleanup temp file
    os.unlink(profile_path)

# Test fixture for a mock workflow
@pytest.fixture
def mock_workflow():
    workflow = MagicMock()
    workflow.scheduling_settings.scheduler = "milp"
    workflow.global_resources = {"_cores": 4, "mem_mb": 8192}
    workflow.dryrun = False
    workflow.touch = False
    workflow.output_settings.quiet = True
    workflow.execution_settings.keep_going = False
    return workflow

# Test fixture for a mock executor plugin
@pytest.fixture
def mock_executor_plugin():
    executor = MagicMock()
    return executor

# Test initialization of MILP scheduler
def test_milp_scheduler_init(mock_workflow, mock_executor_plugin):
    scheduler = MILPJobScheduler(mock_workflow, mock_executor_plugin)
    assert hasattr(scheduler, "node_assignments")
    assert hasattr(scheduler, "_job_start_times")
    assert hasattr(scheduler, "_finished_jobs_history")

# Test default system profile creation
def test_default_system_profile(mock_workflow, mock_executor_plugin):
    scheduler = MILPJobScheduler(mock_workflow, mock_executor_plugin)
    profile = scheduler._get_default_system_profile()
    
    # Verify structure
    assert "clusters" in profile
    assert "local" in profile["clusters"]
    assert "nodes" in profile["clusters"]["local"]
    assert "default" in profile["clusters"]["local"]["nodes"]
    
    # Verify resources
    node = profile["clusters"]["local"]["nodes"]["default"]
    assert "resources" in node
    assert "cores" in node["resources"]
    assert "memory_mb" in node["resources"]
    assert "features" in node
    assert "properties" in node

# Test job hash generation
def test_job_hash_generation(mock_workflow, mock_executor_plugin):
    scheduler = MILPJobScheduler(mock_workflow, mock_executor_plugin)
    
    # Create a mock job
    job = MagicMock()
    job.rule.name = "test_rule"
    job.input = ["input1.txt", "input2.txt"]
    job.output = ["output.txt"]
    job.rule.params = MagicMock(spec_set=True)
    
    # Generate hash
    job_hash = scheduler._get_job_hash(job)
    
    # Verify it's a valid MD5 hash
    assert len(job_hash) == 32
    assert all(c in "0123456789abcdef" for c in job_hash)
    
    # Verify same job gets same hash
    job2 = MagicMock()
    job2.rule.name = "test_rule"
    job2.input = ["input1.txt", "input2.txt"]
    job2.output = ["output.txt"]
    job2.rule.params = MagicMock(spec_set=True)
    
    job_hash2 = scheduler._get_job_hash(job2)
    assert job_hash == job_hash2

# Test job requirements extraction
def test_extract_job_requirements(mock_workflow, mock_executor_plugin):
    scheduler = MILPJobScheduler(mock_workflow, mock_executor_plugin)
    
    # Create config
    config = scheduler._get_default_config()
    
    # Create a mock job with extended specifications
    job = MagicMock()
    job.jobid = "test_job_1"
    job.rule.name = "test_rule"
    job.resources = {"_cores": 2, "mem_mb": 4096, "runtime": 60, "disk_mb": 1000}
    
    # Define job specification
    job_spec = {
        "features": ["cpu", "avx2"],
        "resources": {"input_size_mb": 500, "output_size_mb": 100},
        "properties": {"cpu_flops": 5000000000}
    }
    job.rule.params.job_specification = job_spec
    
    # Extract requirements
    requirements = scheduler._extract_job_requirements(job, config)
    
    # Verify basic fields
    assert requirements["jobid"] == "test_job_1"
    assert requirements["rule_name"] == "test_rule"
    assert requirements["cores"] == 2
    assert requirements["memory_mb"] == 4096
    assert requirements["runtime_minutes"] == 60
    assert requirements["disk_mb"] == 1000
    
    # Verify extended fields
    assert requirements["features"] == ["cpu", "avx2"]
    assert requirements["resources"]["input_size_mb"] == 500
    assert requirements["resources"]["output_size_mb"] == 100
    assert requirements["properties"]["cpu_flops"] == 5000000000

# Test runtime calculation on node
def test_calculate_runtime_on_node(mock_workflow, mock_executor_plugin):
    scheduler = MILPJobScheduler(mock_workflow, mock_executor_plugin)
    
    # Create job requirements
    job_req = {
        "jobid": "test_job_1",
        "rule_name": "test_rule",
        "runtime_minutes": 30,
        "resources": {
            "input_size_mb": 1000,
            "output_size_mb": 500
        }
    }
    
    # Create node data
    node_data = {
        "resources": {
            "cores": 4,
            "memory_mb": 8192
        },
        "properties": {
            "read_mbps": 100,
            "write_mbps": 50
        }
    }
    
    # Calculate runtime
    runtime = scheduler._calculate_runtime_on_node(job_req, node_data)
    
    # Verify it's reasonable
    assert runtime >= job_req["runtime_minutes"]
    assert runtime <= job_req["runtime_minutes"] * 1.5  # Should add small I/O time

# Test with missing properties
def test_calculate_runtime_without_properties(mock_workflow, mock_executor_plugin):
    scheduler = MILPJobScheduler(mock_workflow, mock_executor_plugin)
    
    # Create job requirements without properties
    job_req = {
        "jobid": "test_job_1",
        "rule_name": "test_rule",
        "runtime_minutes": 30
    }
    
    # Create node data without properties
    node_data = {
        "resources": {
            "cores": 4,
            "memory_mb": 8192
        }
    }
    
    # Calculate runtime
    runtime = scheduler._calculate_runtime_on_node(job_req, node_data)
    
    # Verify it returns the base runtime
    assert runtime == job_req["runtime_minutes"]

if __name__ == "__main__":
    pytest.main(["-v", __file__])
