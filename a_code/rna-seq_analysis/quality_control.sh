#!/bin/bash

DIR_IN=/../hihisiv_data_rna-seq/fastq

cd $DIR_IN

for i in $(ls $DIR_IN | grep fastq | awk -F "_[1:2].fastq" '{print $1}' | uniq)
do
echo "====="
echo $i
  fastqc $i"_1.fastq"
  fastqc $i"_2.fastq"
done
