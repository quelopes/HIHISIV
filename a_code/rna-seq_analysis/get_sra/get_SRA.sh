#!/bin/bash

URLFILE=$1 
NUM=8
for URL in `cat $URLFILE`
do aria2c --split=$NUM $URL
#do axel --num-connections=$NUM $URL
done
