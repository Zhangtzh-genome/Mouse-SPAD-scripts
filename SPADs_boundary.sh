#!/bin/bash
#### find the left boundary point and right boundary point of a SPAD.

usage="sh <script.sh> <path> <file> "

if [ $# -ne 2 ]; then
        echo "$usage"
        exit
fi


dataplcae=$1
File=$2


#cd ${dataplcae}/${stage}
cd ${dataplcae}
awk -v OFS="\t" '{print $1,$2-1,$2+1,"LB","1","+"}' ${File} > ${File}.LB.bed
awk -v OFS="\t" '{print $1,$3-1,$3+1,"RB","2","-"}' ${File} > ${File}.RB.bed
cat ${File}.LB.bed ${File}.RB.bed > ${File}.boundary.bed
bedtools sort -i ${File}.boundary.bed > ${File}.boundary.sort.bed
rm ${File}.LB.bed ${File}.RB.bed ${File}.boundary.bed






