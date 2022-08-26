#!/usr/bin/env bash

##########################################################################################
# From joint-called VCF files generated by DELLY or Manta, generate per-individual files #
##########################################################################################

family_list_file=$1
metadata_file=$2

module load bcftools/1.10.2
mkdir -p individual/DELLY individual/Manta TSV/DELLY TSV/Manta
for family in $(cat $family_list_file); do
        print_status "Now working on family $family...\n"
        samples_in_family=$(awk -F $'\t' -v family=$family '$1 == family && $8 != "Complete Genomics"' $metadata_file | cut -f 2)

        for sample in $samples_in_family; do
                for caller in DELLY Manta; do
                        bcftools view --min-ac 1 --samples $sample --output-type z --output-file individual/$caller/$sample.vcf.gz calls/$caller/$family.$caller.vcf.gz # The option --min-ac 1 discards reference calls
                        tabix individual/$caller/$sample.vcf.gz
                        format_vcf_to_tab.$caller.cyvcf2.py --modified-parse individual/$caller/$sample.vcf.gz TSV/$caller/$sample
                        head -n 1 TSV/$caller/${sample}_ALL.tab > TSV/$caller/${sample}.tsv
                        tail -n +2 TSV/$caller/${sample}_ALL.tab | sort -k 2,2V -k 3,3n -k 4,4n >> TSV/$caller/${sample}.tsv
                        rm -f TSV/$caller/${sample}_ALL.tab
                        bgzip TSV/$caller/${sample}.tsv
                        tabix --sequence 2 --begin 3 --end 4 TSV/$caller/${sample}.tsv.gz
                done
        done
done
