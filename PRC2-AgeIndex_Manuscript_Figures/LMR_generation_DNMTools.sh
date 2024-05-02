#!/bin/bash

#SBATCH --job-name=MergeLMRs_ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=11
#SBATCH --mem=200GB
#SBATCH --time=48:00:00

#Daniel J. Simpson 


module load dnmtools/1.2.2
module load samtools/1.9

mkdir Merged_sym


#NEO only HMRs/LMRs

dnmtools merge -o Merged_sym/NEO2_allPs.meth symCpGs/NEO2_P2_sym.meth symCpGs/NEO2_P5_sym.meth symCpGs/NEO2_P8_sym.meth

mkdir Merged_LMRs

#Finding HMRs/LMRs
dnmtools hmr -v -o Merged_LMRs/NEO2_allPs.hmr Merged_sym/NEO2_allPs.meth


#OLD only HMRs/LMRs

dnmtools merge -o Merged_sym/OLD3_allPs.meth symCpGs/OLD3_P2_sym.meth symCpGs/OLD3_P5_sym.meth symCpGs/OLD3_P8_sym.meth

mkdir Merged_LMRs

#Finding HMRs/LMRs
dnmtools hmr -v -o Merged_LMRs/OLD3_allPs.hmr Merged_sym/OLD3_allPs.meth


#LMRs/HMRs on ALL Fibroblast samples

dnmtools merge -o Merged_sym/ALL_fibropass.meth symCpGs/NEO2_P2_sym.meth symCpGs/NEO2_P5_sym.meth symCpGs/NEO2_P8_sym.meth symCpGs/OLD3_P2_sym.meth symCpGs/OLD3_P5_sym.meth symCpGs/OLD3_P8_sym.meth

mkdir Merged_LMRs

#Finding HMRs/LMRs
dnmtools hmr -v -o Merged_LMRs/ALL_fibropass.hmr Merged_sym/ALL_fibropass.meth


