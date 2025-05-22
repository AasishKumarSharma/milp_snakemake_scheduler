"""
Example Snakefile demonstrating how to use the MILP scheduler with job specifications.
"""

# Configuration for the workflow
configfile: "config.yaml"

# Define final target files
SAMPLES = ["sample1", "sample2", "sample3"]
CHROMOSOMES = list(range(1, 23)) + ["X", "Y"]

rule all:
    input:
        expand("results/{sample}/final_report.html", sample=SAMPLES)

# CPU-intensive preprocessing
rule preprocess:
    input:
        fastq="data/{sample}.fastq.gz"
    output:
        bam="results/{sample}/aligned.bam"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=45  # Expected runtime in minutes
    params:
        job_specification={
            "features": ["cpu", "avx2"],
            "resources": {
                "input_size_mb": 2000,
                "output_size_mb": 8000
            },
            "properties": {
                "cpu_flops": 8e11  # 800 GFLOPS required
            }
        }
    shell:
        """
        # Simulate CPU-intensive alignment
        echo "Processing {wildcards.sample}..."
        sleep 10
        mkdir -p $(dirname {output.bam})
        touch {output.bam}
        """

# GPU-accelerated variant calling
rule variant_calling:
    input:
        bam="results/{sample}/aligned.bam"
    output:
        vcf=expand("results/{{sample}}/variants_chr{chrom}.vcf", chrom=CHROMOSOMES)
    threads: 4
    resources:
        mem_mb=8000,
        runtime=30  # Expected runtime in minutes
    params:
        job_specification={
            "features": ["gpu", "cuda"],
            "resources": {
                "gpu_memory_mb": 4000,
                "input_size_mb": 8000,
                "output_size_mb": 500
            },
            "properties": {
                "gpu_flops": 5e12  # 5 TFLOPS required
            }
        }
    shell:
        """
        # Simulate GPU-accelerated variant calling
        echo "Calling variants for {wildcards.sample}..."
        sleep 5
        for chrom in {CHROMOSOMES}; do
            touch results/{wildcards.sample}/variants_chr$chrom.vcf
        done
        """

# Memory-intensive analysis
rule analyze_variants:
    input:
        vcfs=expand("results/{{sample}}/variants_chr{chrom}.vcf", chrom=CHROMOSOMES)
    output:
        report="results/{sample}/final_report.html"
    threads: 2
    resources:
        mem_mb=32000,
        runtime=15  # Expected runtime in minutes
    params:
        job_specification={
            "features": ["high_memory"],
            "resources": {
                "input_size_mb": 1000,
                "output_size_mb": 100
            },
            "properties": {
                "memory_bandwidth_mbps": 20000
            }
        }
    shell:
        """
        # Simulate memory-intensive analysis
        echo "Analyzing variants for {wildcards.sample}..."
        sleep 3
        mkdir -p $(dirname {output.report})
        echo "<html><body><h1>Results for {wildcards.sample}</h1></body></html>" > {output.report}
        """
