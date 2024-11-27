#!/bin/sh

for sample in ${Data}
do

species=mm10
celltype=${stage}
chrtype=SPADs
binsize=100000
smooth=300000
workplace=/your/path/to/data
sof=/your/path/to/Software
ref=/your/path/to/reference/mm10_SNPSplit/C57BL_6NJ_PWK_PhJ_dual_hybrid.based_on_GRCm38_N-masked
split_file=/your/path/to/reference/mm10_SNPSplit/all_PWK_PhJ_SNPs_C57BL_6NJ_reference.chr.based_on_GRCm38.txt
out=/your/path/to/data/${sample}
trim=/your/path/to/Software/Trimmomatic-0.39
trimadp=/your/path/to/Software/Trimmomatic-0.39/adapters
picard=${sof}/picard-tools-2.5.0/picard.jar
chrhmm=${sof}/ChromHMM/ChromHMM.jar
genomesize=${sof}/ChromHMM/CHROMSIZES/${species}.txt
mkdir -p /your/path/to/data/${sample}
cd ${workplace}/${sample}

cat > ${sample}.txt << EOF
${celltype}	${chrtype}	${sample}.picard.bam
EOF

cat > ${sample}.Maternal.txt <<EOF
${celltype}	${chrtype}	${sample}_split.Maternal.sort.bam
EOF

cat > ${sample}.Paternal.txt <<EOF
${celltype}	${chrtype}	${sample}_split.Paternal.sort.bam
EOF

cat > ${sample}_cuttag.sh << EOF


cd ${workplace}/${sample}


##fastqc rawdata##
fastqc  -o ${out}  -t 6 ${sample}_1.fq.gz ${sample}_2.fq.gz

##Trim for QC##
java -jar ${trim}/trimmomatic-0.39.jar  PE -threads 5 -phred33 ${sample}_1.fq.gz ${sample}_2.fq.gz -baseout ${sample}.fq.gz ILLUMINACLIP:${trimadp}/NexteraPE-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:36

###build reference index & mapping in genome
#gunzip ${ref}/${species}.fa.gz
#bowtie2-build  ${ref}/${species}.fa ${ref}/${species}_index
bowtie2 -p ${ppn} -q --phred33 --very-sensitive --end-to-end --no-unal --no-mixed --no-discordant -N 1 -X 3000 -x  ${ref}/mm10_SNPsplit.N-masked.all.index -1 ${workplace}/${sample}/${sample}_1P.fq.gz -2 ${workplace}/${sample}/${sample}_2P.fq.gz -S $sample.sam --un-conc ${workplace}/${sample}/${sample}.unmapped.sam

rm $sample\_1P.fq.gz
rm $sample\_1U.fq.gz
rm $sample\_2P.fq.gz
rm $sample\_2U.fq.gz

samtools view -q 10 -h  ${workplace}/${sample}/$sample.sam > ${workplace}/${sample}/$sample.bam



rm ${workplace}/${sample}/$sample.sam
samtools sort ${workplace}/${sample}/$sample.bam -o ${workplace}/${sample}/$sample.sort.bam

rm ${workplace}/${sample}/$sample.bam

#dulplicate PCR#
##normal##
java -jar ${picard}  MarkDuplicates I=${workplace}/${sample}/$sample.sort.bam O=${workplace}/${sample}/${sample}.picard.bam METRICS_FILE=$sample.picard.bam.metrics REMOVE_DUPLICATES=TRUE
rm ${workplace}/${sample}/$sample.sort.bam
samtools index -b ${workplace}/${sample}/${sample}.picard.bam
bamCoverage -bs ${binsize} --smoothLength ${smooth} --ignoreDuplicates --normalizeUsing RPKM -b ${workplace}/${sample}/${sample}.picard.bam -o ${sample}_bs${binsize}.bw
bamCoverage -bs ${binsize} --ignoreDuplicates --normalizeUsing RPKM -b ${workplace}/${sample}/${sample}.picard.bam -o ${sample}_bs${binsize}.2.bw

java -jar ${chrhmm} BinarizeBam -paired -b ${binsize} ${genomesize} ${workplace}/${sample} ${sample}.txt ${workplace}/${sample}/${sample}_chrhmm
java -jar ${chrhmm} LearnModel -b ${binsize} -color 0,0,255 ${sample}_chrhmm ${workplace}/${sample}/${sample}_chrhmm 2 ${species}

cat ${sample}_chrhmm/${celltype}_2_segments.bed|grep -E 'E2' > ${celltype}_E2_${chrtype}_bs${binsize}.pair.bed
cat ${sample}_chrhmm/${celltype}_2_segments.bed|grep -E 'E1' > ${celltype}_E1_${chrtype}_bs${binsize}.pair.bed
sed -i 's/E2/SPAD/g' ${celltype}_E2_${chrtype}_bs${binsize}.pair.bed
sed -i 's/E1/noSPAD/g' ${celltype}_E1_${chrtype}_bs${binsize}.pair.bed

rm ${sample}_1.fq.gz ${sample}_2.fq.gz

EOF

mkdir -p  ${sample}_chrhmm
mkdir -p  ${sample}_Maternal_chrhmm
mkdir -p  ${sample}_Paternal_chrhmm
sh ${sample}_cuttag.sh

done
