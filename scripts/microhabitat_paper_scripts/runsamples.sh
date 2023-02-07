#!/usr/bin/env bash
# Author: Katharine Prata
# Date created: 07/04/22

# Argument variables
SAMPLE_FILE=$1
LINES=$(cat $SAMPLE_FILE)

for LINE in ${LINES}
do
    python get_measures.py "${LINE}_env_0.1" annotations_all.csv
done

