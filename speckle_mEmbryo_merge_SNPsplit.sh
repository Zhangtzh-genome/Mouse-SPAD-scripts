#!/bin/sh
### This pipeline is for merging differents data from same stage.
### And get maternal and paternal data. 


species=mm10
celltype=${stage}
chrtype=SPADs
binsize=100000
smooth=300000
workplace=/your/path/to/data
sof=/your/path/to/Software
ref=/your/path/to/reference/${species}/${species}.fa
split_ref=/your/path/to/reference/mm10_SNPSplit/C57BL_6NJ_PWK_PhJ_dual_hybrid.based_on_GRCm38_N-masked
split_file=/your/path/to/reference/mm10_SNPSplit/all_PWK_PhJ_SNPs_C57BL_6NJ_reference.chr.based_on_GRCm38.txt
trim=/your/path/to/Software/Trimmomatic-0.39
trimadp=/your/path/to/Software/Trimmomatic-0.39/adapters
picard=${sof}/picard-tools-2.5.0/picard.jar
chrhmm=${sof}/ChromHMM/ChromHMM.jar
genomesize=${sof}/ChromHMM/CHROMSIZES/${species}.txt
mkdir -p /your/path/to/data/${celltype}_merge/${celltype}_merge_SNPsplit
output=/your/path/to/data/${celltype}_merge/${celltype}_merge_SNPsplit

cd ${output}

cat > ${celltype}.Maternal.txt << EOF
${celltype}	${chrtype}	${celltype}_merge_split.Maternal.sort.bam
EOF
cat > ${celltype}.Paternal.txt << EOF
${celltype}	${chrtype}	${celltype}_merge_split.Paternal.sort.bam      
EOF
cat > ${celltype}.merge.txt    << EOF
${celltype}	${chrtype}	${celltype}_merge_split.picard.bam
EOF

cat > ${celltype}_merge_SNPsplit.sh << EOF

cd ${output}
cat ${output}/${celltype}_merge_list
rm ${output}/${celltype}_merge_list

for i in ${Data1} ${Data2} ${Data3} ${Data4}
do
	echo ${i} >> ${output}/${celltype}_merge_list
done

##merge _1.fq.gz
cat \
${workplace}/${Data1}/${Data1}_1.fq.gz \
${workplace}/${Data2}/${Data2}_1.fq.gz \
${workplace}/${Data3}/${Data3}_1.fq.gz \
> ${celltype}_merge_1.fq.gz


##merge _2.fq.gz
cat \
${workplace}/${Data1}/${Data1}_2.fq.gz \
${workplace}/${Data2}/${Data2}_2.fq.gz \
${workplace}/${Data3}/${Data3}_2.fq.gz \
> ${celltype}_merge_2.fq.gz


##Trim for QC##
java -jar ${trim}/trimmomatic-0.39.jar  PE -threads 5 -phred33 ${output}/${celltype}_merge_1.fq.gz ${output}/${celltype}_merge_2.fq.gz -baseout ${output}/${celltype}_merge_fqfile.fq.gz ILLUMINACLIP:${trimadp}/NexteraPE-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:36

rm ${celltype}_merge_1.fq.gz 
rm ${celltype}_merge_2.fq.gz

##Map reads tp SNPsplit refgenome##
bowtie2 -p ${ppn} -q --phred33 --very-sensitive --end-to-end --no-unal --no-mixed --no-discordant -N 1 -X 3000 -x ${split_ref}/mm10_SNPsplit.N-masked.all.index -1 ${output}/${celltype}_merge_fqfile_1P.fq.gz -2 ${output}/${celltype}_merge_fqfile_2P.fq.gz -S ${output}/${celltype}_merge_split.sam

rm ${celltype}_merge_fqfile_1P.fq.gz
rm ${celltype}_merge_fqfile_2P.fq.gz
rm ${celltype}_merge_fqfile_1U.fq.gz
rm ${celltype}_merge_fqfile_2U.fq.gz

samtools view -q 10 -h ${output}/${celltype}_merge_split.sam > ${output}/${celltype}_merge_split.bam

rm ${output}/${celltype}_merge_split.sam

samtools sort ${output}/${celltype}_merge_split.bam -o ${output}/${celltype}_merge_split.sort.bam

rm ${output}/${celltype}_merge_split.bam

#dulplicate PCR#
##normal##
java -jar ${picard}  MarkDuplicates I=${output}/${celltype}_merge_split.sort.bam O=${output}/${celltype}_merge_split.picard.bam METRICS_FILE=${output}/${celltype}_merge_split.picard.metrics REMOVE_DUPLICATES=TRUE

samtools index ${output}/${celltype}_merge_split.picard.bam
bamCoverage -bs ${binsize} --smoothLength ${smooth} --normalizeUsing RPKM --ignoreDuplicates -b ${output}/${celltype}_merge_split.picard.bam -o ${celltype}_merge.bw
bamCoverage -bs ${binsize} --normalizeUsing RPKM --ignoreDuplicates -b ${output}/${celltype}_merge_split.picard.bam -o ${celltype}_merge.2.bw

java -jar ${chrhmm} BinarizeBam -paired -b ${binsize} ${genomesize} ${output} ${output}/${celltype}.merge.txt ${output}/${celltype}_merge_chrhmm
java -jar ${chrhmm} LearnModel -b ${binsize} -color 0,0,255 ${output}/${celltype}_merge_chrhmm ${output}/${celltype}_merge_chrhmm 2 ${species}

cat ${output}/${celltype}_merge_chrhmm/${celltype}_2_segments.bed|grep -E 'E2' > ${output}/${celltype}_merge_E2_${chrtype}_bs${binsize}.bed
sed -i 's/E2/SPAD/g' ${output}/${celltype}_merge_E2_${chrtype}_bs${binsize}.bed

##SNPsplit for Maternal&Paternal##
SNPsplit --paired --snp_file ${split_file} ${output}/${celltype}_merge_split.picard.bam

rm ${output}/${celltype}_merge_split.sort.bam

mv ${celltype}_merge_split.picard.genome1.bam ${celltype}_merge_split.Maternal.bam
mv ${celltype}_merge_split.picard.genome2.bam ${celltype}_merge_split.Paternal.bam

samtools sort ${celltype}_merge_split.Maternal.bam -o ${celltype}_merge_split.Maternal.sort.bam
samtools index ${celltype}_merge_split.Maternal.sort.bam
samtools sort ${celltype}_merge_split.Paternal.bam -o ${celltype}_merge_split.Paternal.sort.bam
samtools index ${celltype}_merge_split.Paternal.sort.bam

rm ${celltype}_merge_split.Maternal.bam
rm ${celltype}_merge_split.Paternal.bam

bamCoverage -bs ${binsize} --smoothLength ${smooth} --normalizeUsing RPKM --ignoreDuplicates -b ${celltype}_merge_split.Maternal.sort.bam -o ${celltype}_merge_split.Maternal.bw
bamCoverage -bs ${binsize} --normalizeUsing RPKM --ignoreDuplicates -b ${celltype}_merge_split.Maternal.sort.bam -o ${celltype}_merge_split.Maternal.2.bw

java -jar ${chrhmm} BinarizeBam -paired -b ${binsize} ${genomesize} ${output} ${output}/${celltype}.Maternal.txt ${output}/${celltype}_Maternal_chrhmm
java -jar ${chrhmm} LearnModel -b ${binsize} -color 0,0,255 ${output}/${celltype}_Maternal_chrhmm ${output}/${celltype}_Maternal_chrhmm 2 ${species}

cat ${output}/${celltype}_Maternal_chrhmm/${celltype}_2_segments.bed|grep -E 'E2' > ${output}/${celltype}_Maternal_E2_${chrtype}_bs${binsize}.bed
sed -i 's/E2/SPAD/g' ${output}/${celltype}_Maternal_E2_${chrtype}_bs${binsize}.bed

bamCoverage -bs ${binsize} --smoothLength ${smooth} --normalizeUsing RPKM --ignoreDuplicates -b ${celltype}_merge_split.Paternal.sort.bam -o ${celltype}_merge_split.Paternal.bw
bamCoverage -bs ${binsize} --normalizeUsing RPKM --ignoreDuplicates -b ${celltype}_merge_split.Paternal.sort.bam -o ${celltype}_merge_split.Paternal.2.bw

java -jar ${chrhmm} BinarizeBam -paired -b ${binsize} ${genomesize} ${output} ${output}/${celltype}.Paternal.txt ${output}/${celltype}_Paternal_chrhmm
java -jar ${chrhmm} LearnModel -b ${binsize} -color 0,0,255 ${output}/${celltype}_Paternal_chrhmm ${output}/${celltype}_Paternal_chrhmm 2 ${species}

cat ${output}/${celltype}_Paternal_chrhmm/${celltype}_2_segments.bed|grep -E 'E2' > ${output}/${celltype}_Paternal_E2_${chrtype}_bs${binsize}.bed
sed -i 's/E2/SPAD/g' ${output}/${celltype}_Paternal_E2_${chrtype}_bs${binsize}.bed
for k in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19
do
	awk '{if ($1 == "chr'$k'") print $0}' ${output}/${celltype}_Paternal_E2_${chrtype}_bs${binsize}.bed >> ${output}/${celltype}_Paternal_E2_${chrtype}_bs${binsize}.final.bed
done
EOF

mkdir -p  ${output}/${celltype}_Maternal_chrhmm
mkdir -p  ${output}/${celltype}_Paternal_chrhmm
mkdir -p  ${output}/${celltype}_merge_chrhmm
sh  ${celltype}_merge_SNPsplit.sh


