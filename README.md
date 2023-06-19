# Assign reads to strain from databse

The fastq files from the pacbio run are first renamed and stored in the appropriate directory. The rest of the analysis (including trimming and filtering) is done within the R chunks written in Analysis_and_report.Rmd.

## Old approach

<This is the older approach used but it was not worth the computational time as I find that the dada2 approach works well>

This repository is the code used to do the initial processing of the pacbio reads obtained from pacbio sequencing of a defined community of about 10 strains per individual. Their genomes are sequenced and 16S region sequence known. The snakemake pipeline reads the raw reads, renamed them as specified in the config file and then reads through them and assigned each read to a 16S sequence in the database fasta file after subsetting it to the strains used in the experiment as specified by the config file. Final summary of counts can be found per sample in the directory called 03_assign_reads_to_strain.

To run snakemake in curnagl, use snakemake -p --use-conda --conda-prefix /scratch/aprasad/built-envs/ --conda-frontend mamba --profile slurm --restart-times 0 -r --cluster-cancel scancel --keep-going --rerun-incomplete --rerun-triggers mtime -n (snakemake v7.8.5)

These files can then be processed using the PCR data and metadata specified in the Rmd file used to make the final report. The path to nas on the Rmd file needs to be changed according to the computing system that the file is run or knitted in. If this computing system does not have nas_recherche mounted, the script will not work. If running interactively, make sure to do so the conda environment specified in the yaml file. In curnagl, this is at /work/FAC/FBM/DMF/pengel/spirit/aprasad/Miniconda3/envs/pacbio_ampli_env

The following commands are useful for monitoring progress. cd 03_assign_reads_to_strain while true; do for file in *.progress; do echo $file; cat $file | tr -d "\n"; echo " "; done; sleep 1; clear; done for file in *summary*; do echo $file; cat $file; echo; echo; done