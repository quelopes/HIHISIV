#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH -p cpu_dev
#SBATCH -J vitradet

module load bowtie2/2.3
module load rsem/1.3


DIR_IN=/../hihisiv_results/*genes.results


rsem-generate-data-matrix $DIR_IN > geneMat.txt
rsem-run-ebseq geneMat.txt 25,25 GeneMat.results
rsem-control-fdr GeneMat.results 0.05 GeneMat.de.txt









