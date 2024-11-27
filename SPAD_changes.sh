#!/bin/bash

#A script for defining 7 types of SPADs.

usage="sh <script.sh> <stage1_later> <stage2_earlier>"
if [ $# -ne 2 ]; then
        echo "$usage"
        exit
fi

s1=$1
s2=$2
## later stage lost SPADs
bedtools intersect -a ${s2}.final_bed3.bed -b ${s1}.final_bed3.bed -wao > ${s2}_${s1}_wao.bed
awk -v OFS="\t" '{if($7==0) print $1,$2,$3}' ${s2}_${s1}_wao.bed > ${s1}_${s2}_lost.bed


bedtools intersect -a ${s1}.final_bed3.bed -b ${s2}.final_bed3.bed -wao > ${s1}_${s2}_wao.bed
awk -v OFS="\t" '{print $0,$7/($3-$2),$7/($6-$5)}' ${s1}_${s2}_wao.bed > ${s1}_${s2}_wao_percent.bed
awk -v OFS="\t" '{if($8>0.75 && $9>0.75) print $1,$2,$3}' ${s1}_${s2}_wao_percent.bed > ${s1}_${s2}_stable.bed
awk -v OFS="\t" '{if($7==0) print $1,$2,$3}' ${s1}_${s2}_wao_percent.bed > ${s1}_${s2}_gained.bed

## First get newly formed SPADs
awk -v OFS="\t" '{if($7==0&&$8==0) print $1,$2,$3}' ${s1}_${s2}_wao_percent.bed > ${s1}_${s2}_new.bed
bedtools intersect -a ${s1}_${s2}_wao_percent.bed -b ${s1}_${s2}_new.bed -v > ${s1}_${s2}_wao_percent_no_new.bed

## Then get several SPADs merged
bedtools groupby -i ${s1}_${s2}_wao_percent_no_new.bed -g 1-3 -c 4,5,6,7,8,9 -o collapse,collapse,collapse,collapse,sum,sum > ${s1}_${s2}_groupby123sum.bed
awk -v OFS="\t" '{if($9>=1.6) print $1,$2,$3}' ${s1}_${s2}_groupby123sum.bed > ${s1}_${s2}_merge.bed

## Next get splits from earlier stage SPAds
bedtools groupby -i ${s1}_${s2}_groupby123sum.bed -g 4-6 -c 1,2,3,7,8,9 -o collapse,collapse,collapse,collapse,sum,sum > ${s1}_${s2}_split.tmp.bed
awk -v OFS="\t" '{if($8>=1.6)print $1,$2,$3}' ${s1}_${s2}_split.tmp.bed > ${s2}_split.tmp.bed
bedtools intersect -a ${s1}_${s2}_groupby123sum.bed -b  ${s2}_split.tmp.bed -wa > ${s1}_${s2}_split.tmp.bed
awk -v OFS="\t" '{print $1,$2,$3}' ${s1}_${s2}_split.tmp.bed > ${s1}_${s2}_split.bed

## Then get stable SPADs
bedtools intersect -a ${s1}_${s2}_groupby123sum.bed -b ${s1}_${s2}_merge.bed  ${s1}_${s2}_split.bed -v > ${s1}_${s2}_groupby123sum2.bed
awk -v OFS="\t" '{if($8>=0.75 && $8<=1 && $9>=0.75 && $9<=1) print $1,$2,$3}' ${s1}_${s2}_groupby123sum2.bed > ${s1}_${s2}_stable.bed

## Then get expand SPADs
awk -v OFS="\t" '{if($8<0.75 && $9>=0.75 && $9<=1) print $1,$2,$3}' ${s1}_${s2}_groupby123sum2.bed > ${s1}_${s2}_expand.bed

## Then get shrink SPADs
awk -v OFS="\t" '{if($8>=0.75 && $8<=1 && $9<0.75) print $1,$2,$3}' ${s1}_${s2}_groupby123sum2.bed > ${s1}_${s2}_shrink.bed

## Fianlly, what left are rearranged SPADs
bedtools intersect -a ${s1}.final_bed3.bed -b ${s1}_${s2}_new.bed ${s1}_${s2}_merge.bed ${s1}_${s2}_split.bed ${s1}_${s2}_stable.bed ${s1}_${s2}_expand.bed ${s1}_${s2}_shrink.bed -v > ${s1}_${s2}_rearrange.bed
