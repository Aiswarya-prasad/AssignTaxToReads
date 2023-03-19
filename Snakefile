"""
This pipeline is written by Aiswarya Prasad (aiswarya.prasad@unil.ch) and is
intended primarily for her own use in her PhD thesis project(s). It is being written
to work on snakemake v6.15.5 and run in a cluster using the slurm profile mentioned
 here (https://github.com/RomainFeron/snakemake-slurm) and minor modifications.

This was run in curnagl (in conda environment called snakemake_7.7) with
snakemake -p --use-conda --conda-prefix /work/FAC/FBM/DMF/pengel/spirit/aprasad/snakemake-conda-envs --conda-frontend mamba --profile slurm --restart-times 0 -r --cluster-cancel scancel --keep-going --rerun-incomplete -n
"""

import os
import itertools

configfile: "config.yaml"

SAMPLES = config["raw_files"].keys()

def convertToMb(string):
    """
    This function can convert text in the form
    xxG to mb
    If it does not end with G, it returns the string
    It does not handle other cases of invalid input
    """
    if string.endswith("G"):
        number = int(string.split("G")[0])
        return(number*1000)
    else:
        return(string)

def convertToSec(string):
    """
    This function can convert text in the form
    D-hh:mm:ss to seconds
    D - # days
    hh # hours
    mm # mins
    ss # secs
    """
    days = string.split("-")[0]
    hrs = string.split("-")[1].split(":")[0]
    min = string.split("-")[1].split(":")[1]
    sec = string.split("-")[1].split(":")[2]
    total = int(sec)
    total = total + 60*int(min)
    total = total + 60*60*int(hrs)
    total = total + 24*60*60*int(days)
    return(total)

rule all:
    input:
        expand("03_assign_reads_to_strain/{sample}_strain_counts.csv", sample=SAMPLES)

rule rename_reads:
    input:
        file = lambda wildcards: [path for path in config["raw_files"][wildcards.sample]]
    output:
        reads_renamed = "01_ReadsRenamed/{sample}_reads.fastq.gz"
    conda: 
        "envs/pacbio-ampli-env.yaml"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-0:10:00")
    threads: 2
    resources:
        mem_mb = 8000
    shell:
        """
        cp -n {input.file} {output.reads_renamed}
        """

rule fastqc:
    input:
        reads = rules.rename_reads.output.reads_renamed
    output:
        html="02_FastQCBeforeTrimming/{sample}_fastqc.html",
        zip="02_FastQCBeforeTrimming/{sample}_fastqc.zip"
    threads: 2
    conda: 
        "envs/pacbio-ampli-env.yaml"
    params:
        outdir = "02_FastQCBeforeTrimming",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-1:00:00")
    resources:
        mem_mb = 8000
    shell:
        """
        fastqc -t {threads} {input.reads} -o {params.outdir}
        """

rule identify_write_positions:
    input:
        fasta_sequences = config["database"]["all_sequences"]
    output:
        pickle_file = "database/16S_sequences/unique_position_info.pickle",
        alignment = "database/16S_sequences/16S_aligned.fasta"
    threads: 2
    conda: 
        "envs/pacbio-ampli-env.yaml"
    params:
        list_strains = [strain for strain in config["database"]["subset_genomes"]],
        database_path = config["database"]["database_path"], # alignment and pickle file will be found inside by script
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-1:10:00")
    resources:
        mem_mb = 16000
    log:
        "database/16S_sequences/identify_write_positions.log"
    benchmark:
        "database/16S_sequences/identify_write_positions.benchmark"
    shell:
        """
        python3 scripts/identify_unique_position_set.py --database_path {params.database_path} \
                                        --input_fasta {input.fasta_sequences} \
                                        --output_file {output.pickle_file} \
                                        --subset true
                                        --list_strains {parama.list_strains} | tee {log}
        """

rule assign_reads_to_strains:
    input:
        reads_file = rules.rename_reads.output.reads_renamed,
        pickle_file = "database/16S_sequences/unique_position_info.pickle",
        alignment = "database/16S_sequences/16S_aligned.fasta"
    output:
        csv = "03_assign_reads_to_strain/{sample}_strain_counts.csv",
        summary = "03_assign_reads_to_strain/{sample}_summary.txt"
    threads: 2
    conda: 
        "envs/pacbio-ampli-env.yaml"
    params:
        database_path = config["database"]["database_path"], # alignment and pickle file will be found inside by script
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-1:10:00")
    resources:
        mem_mb = 16000
    log:
        "database/16S_sequences/{sample}_assign_reads_to_strains.log"
    benchmark:
        "database/16S_sequences/{sample}_assign_reads_to_strains.benchmark"
    shell:
        """
        python3 scripts/assign_reads_to_strains.py --database_path {params.database_path} \
                                           --reads_file {input.reads_file} \
                                           --summary_file_path {output.summary} \
                                           --output_path {output.csv} \
                                           --sample_name {wildcards.sample}
        """

