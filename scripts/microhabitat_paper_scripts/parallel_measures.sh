#!/bin/sh

#  parallel_measures.sh
#  
#
#  Created by Katharine Prata on 7/4/22.
# TODO: use env_dist argument variable

# Argument variables
SAMPLE_FILE=$1

#cat ${SAMPLE_FILE} | time parallel -j+0 --eta 'python get_measures.py {}_env_0.5 annotations_all.csv'
cat ${SAMPLE_FILE} | time parallel -j4 --eta 'python get_measures.py {}_env_0.5 annotations_all.csv'
