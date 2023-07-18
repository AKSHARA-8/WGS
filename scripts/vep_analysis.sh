#! /usr/bin/bash 

input_vcf=$1
output_path=$2

source deactivate
vep -i ${input_vcf} --cache -o ${output_path}_vep.txt --offline --tab --force_overwrite --af_1kg --af --af_gnomadg --pubmed --sift b --canonical --hgvs --shift_hgvs 1
filter_vep -i ${output_path}_vep.txt -o ${output_path}_filtered.txt --filter "(CANONICAL is YES) and (AF < 0.01 or not AF)" --force_overwrite
grep -v "##" ${output_path}_filtered.txt > ${output_path}_vep_delheaders.txt
source activate new_base
