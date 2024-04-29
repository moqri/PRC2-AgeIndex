#!/bin/bash

#SBATCH --job-name=WGBS_bis-methpip
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=400GB
#SBATCH --time=48:00:00
#SBATCH --array=1-6%6

# Daniel J. Simpson

#unified clean script so that I don't have to rewrite this everytime (hopefully)

#list of samples processed via array
index=$1
fileID=`sed "${SLURM_ARRAY_TASK_ID}p;d" ${index} | awk '{print $1}'`

#bismark formatted genome index
bisIndex=$2

#path to genome fasta file
Index=$3

#working directory
workdir=$4


module load trim-galore/0.6.7
#running local, more recent version of Bismark
export PATH="/home/dsimps93/Bismark-0.24.0:$PATH"
module load cutadapt/3.4
module load fastqc/0.11.9
module load bowtie2/2.5.0
module load samtools/1.9
module load dnmtools/1.2.2
module load ucsc_tools/309


mkdir $workdir/Trimmed

cd $TMPDIR

echo '************************************'
echo Start ${fileID} Trimming:
echo '************************************'
date


##added gzip so that output should be gzipped
trim_galore --paired --gzip --cores 8  $workdir/Fastas/${fileID}*1.fastq $workdir/Fastas/${fileID}*2.fastq --output_dir $workdir/Trimmed


###Bismark Alignment


echo '************************************'
echo Start ${fileID} bismark alignment:
echo '************************************'
date

mkdir $workdir/Bis_Mapped

# run alignment
bismark --multicore 8 --phred33-quals -N 0 -L 20 ${bisIndex} -1 $workdir/Trimmed/${fileID}*1*.fq.gz -2 $workdir/Trimmed/${fileID}*2*.fq.gz  -o $workdir/Bis_Mapped



echo '************************************'
echo Start ${fileID} conversion:
echo '************************************'
date

mkdir $workdir/DNMTools_SAMs

#Converting bismark bam file to dnmtools format
dnmtools format -v -f bismark -o $workdir/DNMTools_SAMs/${fileID}_formatted.sam  $workdir/Bis_Mapped/${fileID}_*.bam

mkdir $workdir/sortedSams


echo '************************************'
echo Start ${fileID} Sorting:
echo '************************************'
date

#Sorting Sam files
samtools sort -O sam -o $workdir/sortedSams/${fileID}_input-sorted.sam $workdir/DNMTools_SAMs/${fileID}_formatted.sam


mkdir $workdir/uniqueSams


echo '************************************'
echo Start ${fileID} Uniqing:
echo '************************************'
date

#Uniqing
dnmtools uniq -v $workdir/sortedSams/${fileID}_input-sorted.sam $workdir/uniqueSams/${fileID}_out-sorted.sam

#Calculate conversion rate, will comment out now to save time
#dnmtools bsrate [OPTIONS] -c <chroms> <input.sam>

echo '************************************'
echo Start ${fileID} methcounts:
echo '************************************'
date

mkdir $workdir/methcounts

#Calculating counts. Doing CpG Context only to save time
dnmtools counts -v -cpg-only -c $Index -o $workdir/methcounts/${fileID}.meth $workdir/uniqueSams/${fileID}_out-sorted.sam


echo '************************************'
echo Start ${fileID} Sym:
echo '************************************'
date

mkdir $workdir/symCpGs

#Output to symmetryic CpGs
dnmtools sym -o $workdir/symCpGs/${fileID}_sym.meth $workdir/methcounts/${fileID}.meth


echo '************************************'
echo Start ${fileID} Bigwigs:
echo '************************************'
date

mkdir $workdir/meth_bigwigs


awk -v OFS="\t" '{print $1, $2, $2+1, $4":"$6, $5, $3}'  $workdir/symCpGs/${fileID}_sym.meth >  $workdir/symCpGs/${fileID}_symmeth.bed

cut -f 1-3,5 $meths/${fileID}_symmeth.bed | wigToBigWig /dev/stdin /labs/vsebast/DJS/genomes/HomoSap/Primary_ver/hg38_ensembEd.chrom.sizes.txt $workdir/meth_bigwigs/${fileID}_meth.bw


