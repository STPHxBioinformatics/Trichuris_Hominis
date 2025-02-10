#!/bin/bash

#SBATCH --job-name=InterproTrichurisAnnotationMax                   #This is the name of your job
#SBATCH --cpus-per-task=16                #This is the number of cores reserved
#SBATCH --mem-per-cpu=8G              #This is the memory reserved per core.
#Total memory reserved: 128GB

#SBATCH --time=1-00:00:00        #This is the time that your task will run
#SBATCH --qos=1day           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/scicore/home/schpie00/baer0006/ONT/de_novo_annotation_pipeline/interpro/stdoutinterpro1    #These are the STDOUT and STDERR files
#SBATCH --error=/scicore/home/schpie00/baer0006/ONT/de_novo_annotation_pipeline/interpro/stderrinterpro1
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=max.baer@swisstph.ch        #You will be notified via email when your task ends or fails

#This job runs from the current working directory


#Remember:
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $SLURM_JOBID stores the ID number of your job.


#load your required modules below
ml load InterProScan/5.63-95.0-foss-2021a;
ml load InterProScan_data/5.63-95.0-foss-2021a

#export your required environment variables below
#################################################


#add your command lines below
interproscan.sh -i Trichuris_cote_divoire_freeze_genes_no_special.fasta -t n -dp -iprlookup -goterms

