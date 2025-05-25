# MILP Snakemake Scheduler Plugin

[![PyPI version](https://img.shields.io/pypi/v/milp-snakemake-scheduler)](https://pypi.org/project/milp-snakemake-scheduler/) [![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

A **pip-installable** scheduler plugin for [Snakemake](https://snakemake.readthedocs.io/) that uses a Mixed-Integer Linear Programming (MILP) scheduler for Snakemake to optimally analyze the jobs for heterogeneous compute resources.

- **`milp`**: schedules only ready jobs via MILP.
- **`milp-ext`**: schedules across all pending jobs using MILP + critical-path analysis.

---
## Features

* **Resource-Aware Scheduling**: Considers CPU cores, memory, GPU, custom resources, and I/O constraints.
* **Critical-Path Analysis**: Optionally detect and prioritize the longest dependency chain.
* **Historical Estimation**: Learns from past executions to refine runtime and I/O size predictions.
* **Multi-Objective Optimization**: Balances makespan vs. energy consumption with configurable weights.
* **Graceful Fallbacks**: Falls back to greedy or ILP strategies if MILP fails or times out.
* **Plugin Auto-Discovery**: Integrates seamlessly via Snakemake entry points.

---

## Installation

### From PyPI

```bash
pip install milp-snakemake-scheduler
```

### From Source

```bash
git clone https://github.com/AasishKumarSharma/milp_snakemake_scheduler.git
cd milp_snakemake_scheduler
pip install .
```

For development (including test dependencies):

```bash
pip install -e .[dev]
```

---

## Quick Start

1. **Create a System Profile**: `system_profile.json` describes your clusters/nodes.

2. **(Optional) Customize**: `scheduler_config.yaml` to adjust estimation, objectives, and fallbacks.

3. **Run Snakemake**:

   3.1.  Run in **MILP-only** mode for ready jobs:

   ```bash
   snakemake --scheduler milp \
             --cores 8 \
             --scheduler-config scheduler_config.yaml \
             --system-profile system_profile.json \
             --jobs 16
   ```

   3.2.  Run in **Extended MILP** mode with critical-path for a ready jobs in a workflow:
   
   ```bash
   snakemake --scheduler milp-ext \
             --cores 8 \
             --scheduler-config scheduler_config.yaml \
             --system-profile system_profile.json \
             --jobs 16
   ```


4. **Dry-Run**:

   ```bash
   snakemake --scheduler milp -n --cores 8 --scheduler-config scheduler_config.yaml --system-profile system_profile.json
   ```

---

## Configuration

### system\_profile.json

Defines compute clusters and nodes. Example:

```json
{
  "clusters": {
    "local": {
      "nodes": {
        "default": {
          "resources": {"cores": 8, "memory_mb": 32768, "gpu_count": 0},
          "features": ["cpu", "x86_64", "avx2"],
          "properties": {"cpu_flops": 1e11, "memory_bandwidth_mbps": 25600}
        }
      }
    },
    "gpu_cluster": { ... },
    "high_memory": { ... }
  }
}
```

### scheduler\_config.yaml

Controls scheduler behavior. Example:

```yaml
scheduler:
  type: milp
  paths:
    system_profile: system_profile.json
  estimation:
    auto_estimate_file_sizes: true
    history:
      enabled: true
      adaptation_weight: 0.7
  optimization:
    objective_weights:
      makespan: 0.8
      energy: 0.2
    time_limit_seconds: 30
    fallback: greedy
```

Key options:

* `estimation.auto_estimate_file_sizes`: infer I/O if not provided.
* `optimization.time_limit_seconds`: MILP solver cutoff.
* `optimization.fallback`: `greedy` or `ilp` on failure.

---

## Example Snakefile

See `example_snakefile.py` in this repo. A minimal snippet:

```python
configfile: "scheduler_config.yaml"

rule all:
    input: "results/output.txt"

rule demo:
    output: "results/output.txt"
    threads: 2
    resources:
        mem_mb=1024,
        runtime=5
    params:
        job_specification={
            "features": ["cpu"],
            "resources": {"input_size_mb": 10, "output_size_mb": 5},
            "properties": {"cpu_flops": 1e9}
        }
    shell:
        "echo Hello > {output}"
```

---

## Testing

Run tests with:

```bash
pytest tests/
```

Include additional tests in `tests/` following the existing patterns.

---

## Contributing

1. Fork and clone the repo
2. Create a branch: `git checkout -b feature/new-flag`
3. Implement changes and add tests
4. Ensure all tests pass: `pytest`
5. Submit a pull request

Please follow [Conventional Commits](https://www.conventionalcommits.org/) for commit messages.

---

## License

This project is licensed under the [MIT License](LICENSE).

