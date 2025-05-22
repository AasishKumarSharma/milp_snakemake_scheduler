from setuptools import setup, find_packages

setup(
    name="milp-snakemake-scheduler",
    version="0.1.1",
    author="Your Name",
    description="MILP‑based job scheduler plugin for Snakemake",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[
        "pulp>=2.0",
        "networkx>=2.5",
    ],
    entry_points={
        "snakemake.scheduler_plugins": [
            "milp = milp_scheduler.scheduler:milp_scheduler_factory",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
