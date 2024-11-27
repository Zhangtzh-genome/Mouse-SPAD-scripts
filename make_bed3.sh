#!/bin/bash
### After p-value checking, make 3-column SPAD region files

for name in m8C
	#mOoGV mZygote mE2C mL2C m4C m8C mMorula mBlastocyst mESC mESC_new
do
	awk -v OFS="\t" '{print $1,$2,$3}' ${name}.p-value.final.noXYchr.bed > ${name}.final_bed3.bed
done
