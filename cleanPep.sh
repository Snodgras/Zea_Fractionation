#!/bin/bash
file=$1
if [ ${file} != Sb313.pep.faa ]
then
	cut -f 1 -d "_" ${file} > temp
	mv temp ${file}
fi
