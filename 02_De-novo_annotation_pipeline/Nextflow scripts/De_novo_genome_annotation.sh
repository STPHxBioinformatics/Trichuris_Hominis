#!/bin/bash

#SBATCH --job-name=TrichurisAnnotationMax                   #This is the name of your job
#SBATCH --cpus-per-task=8                #This is the number of cores reserved
#SBATCH --mem-per-cpu=4G              #This is the memory reserved per core.
#Total memory reserved: 64GB

#SBATCH --time=06:00:00        #This is the time that your task will run
#SBATCH --qos=6hours           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/scicore/home/schpie00/baer0006/ONT/de_novo_annotation_pipeline/stdout2    #These are the STDOUT and STDERR files
#SBATCH --error=/scicore/home/schpie00/baer0006/ONT/de_novo_annotation_pipeline/stderr2
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=max.baer@swisstph.ch        #You will be notified via email when your task ends or fails

#This job runs from the current working directory


#Remember:
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $SLURM_JOBID stores the ID number of your job.


#load your required modules below
ml load Java/11.0.3_7

#export your required environment variables below
#################################################


#add your command lines below
nextflow run De_novo_genome_annotation.nf -resume
