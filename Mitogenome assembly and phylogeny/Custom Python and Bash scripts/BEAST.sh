#!/bin/bash

#SBATCH --job-name=TrichurisBeastMax                  #This is the name of your job
#SBATCH --cpus-per-task=32                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=8G              #This is the memory reserved per core.
#Total memory reserved: 156GB

#SBATCH --time=1-00:00:00        #This is the time that your task will run
#SBATCH --qos=1day          #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/scicore/home/schpie00/baer0006/Trichuris_Trichiura_CI/mitochondrial_baiting_and_analysis_pipeline/beast_runs/run_9_JC69_concatenated/stdout_7    #These are the STDOUT and STDERR files
#SBATCH --error=/scicore/home/schpie00/baer0006/Trichuris_Trichiura_CI/mitochondrial_baiting_and_analysis_pipeline/beast_runs/run_9_JC69_concatenated/stderr_7
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=max.baer@swisstph.ch        #You will be notified via email when your task ends or fails

#This job runs from the current working directory


#Remember:
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $SLURM_JOBID stores the ID number of your job.


#load your required modules below
ml load beagle-lib/3.1.2-GCCcore-10.3.0  
#export your required environment variables below
#################################################


#add your command lines below
/scicore/home/schpie00/baer0006/Trichuris_Trichiura_CI/Nextflow/beast/beast/bin/beast -overwrite -threads 32 concatenated_JC69_model_for_rate.xml