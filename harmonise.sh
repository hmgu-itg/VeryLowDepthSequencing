#!/bin/bash

case=$1
refpanel=$2

tabix -f -p vcf $case

basecase=$(echo $case | sed 's/.vcf.gz//')
baserefp=$(echo $refpanel | sed 's/.vcf.gz//')

bn=$(basename $baserefp)

echo "Extracting positions from reference panel..."
zgrep -v '#' $refpanel | awk 'BEGIN{OFS="\t";}{print $1,$2-1,$2}' > toalign.bed
echo "Querying reference..."
/nfs/team144/software/bedtools2/bin/fastaFromBed -fi /lustre/scratch113/teams/zeggini/ReferenceData/ref_human/hs37d5.fa -bed toalign.bed -fo toalign.fromb37.fa
echo "Matching against reference alleles in panel..."
paste <(grep -v '>' toalign.fromb37.fa) <(zgrep -v '#' $refpanel | awk '{print $4}') > discordance
echo "Aligning..."
lastfile=${baserefp}.aligned.vcf.gz
~/alignVCF2.pl $refpanel | bgzip > $lastfile
tabix -f $lastfile
echo "Norming (collapsing indels)..."
/nfs/team144/software/bcftools/bcftools norm -f /lustre/scratch113/teams/zeggini/ReferenceData/ref_human/hs37d5.fa -c w -m +any -D -Oz -o ${baserefp}.aligned.normed.vcf.gz $lastfile
tabix -f -p vcf ${baserefp}.aligned.normed.vcf.gz
echo "Merging..."
/nfs/team144/software/bcftools/bcftools merge -m all ${baserefp}.aligned.normed.vcf.gz $case -Oz -o ${basecase}.${bn} 
tabix -f -p vcf ${basecase}.${bn} 
echo "Demerging+filtering reference..."
zcat $refpanel | head -n 5000 | grep -m 1 '#CHROM' | cut -f10- | sed 's/\t/\n/g'> ref.samples
/nfs/team144/software/bcftools/bcftools view -S ref.samples ${basecase}.${bn} | grep -v './.' | bgzip > ${baserefp}.aligned.normed.demerged.filtered.vcf.gz
# tabix -f -p vcf ${baserefp}.aligned.normed.demerged.filtered.vcf.gz
echo "Demerging+filtering case..."
zcat $case | head -n 5000 | grep -m 1 '#CHROM' | cut -f10- | sed 's/\t/\n/g'> case.samples
/nfs/team144/software/bcftools/bcftools view -S case.samples ${basecase}.${bn} | grep -e PL -e '#' | bgzip > ${basecase}.demerged.filtered.vcf.gz 
tabix -f -p vcf ${basecase}.demerged.filtered.vcf.gz
echo "Removing null Pl..."
perl ~/xg1/nomorenullpl3.pl ${basecase}.demerged.filtered.vcf.gz ${basecase}.nonullPL.demerged.filtered.vcf.gz
tabix -f -p vcf ${basecase}.nonullPL.demerged.filtered.vcf.gz
cat <(zcat ${baserefp}.aligned.normed.demerged.filtered.vcf.gz | head -5000 | grep '^#'| head -10) <(zcat ${basecase}.nonullPL.demerged.filtered.vcf.gz | head -5000 | grep '^##'| grep -e INFO -e FORMAT) <(zcat ${baserefp}.aligned.normed.demerged.filtered.vcf.gz | head -5000 | grep '^#'| tail -n+11) > newheader
tabix -r newheader ${baserefp}.aligned.normed.demerged.filtered.vcf.gz > ${baserefp}.aligned.normed.demerged.filtered.reheaded.vcf.gz
tabix -f -p vcf ${baserefp}.aligned.normed.demerged.filtered.reheaded.vcf.gz
echo "FINISHED!"
