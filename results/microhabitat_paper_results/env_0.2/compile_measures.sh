#!/bin/sh

#  compile_measures.sh
#
#
#  Created by Katharine Prata on 14/6/22.
#

folder=$1

FILENAMES=`ls ${folder}/`

for file in $FILENAMES
do cat ${folder}/${file} | tail -n 1 >> ${folder}.txt
done
