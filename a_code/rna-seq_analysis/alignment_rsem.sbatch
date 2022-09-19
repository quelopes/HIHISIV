#!/bin/bash

#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH -p cpu_dev
#SBATCH -J vitradet

module load bowtie2/2.3
module load perl/5.20
module load rsem/1.3
module load gnu-parallel/20190422

DIR_IN=/../data/hiv_sra
DIR_OUT=/../hihisiv_results
DIR_GENOME=/../hg19/WholeGenome

NTASKS=4
NTHREADS=24

cd $DIR_IN

parallel --jobs $NTASKS srun -N 1 -n 1 -c 24 rsem-calculate-expression --paired-end --bowtie2 \
 -p $NTHREADS --output-genome-bam --keep-intermediate-files --append-names --estimate-rspd \
 --time {}"_1.fastq" {}"_2.fastq" $DIR_GENOME/rsem_ref_hg19 $DIR_OUT/{} ::: \
 ls *_1.fastq | awk -F "_2.fastq" '{print $1}' | uniq
