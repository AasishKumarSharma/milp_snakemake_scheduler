"""Scheduler extension module for Snakemake."""

from snakemake.scheduler_plugins.milp_scheduler.scheduler import MILPJobScheduler

__all__ = ["MILPJobScheduler"]
