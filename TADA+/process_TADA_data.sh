#!/usr/bin/env bash

########################################################################################################################
# Wrapper script to perform all of the sample comparisons and exclusions, make the tables, etc. for the TADA+ analysis #
########################################################################################################################

# Identify excluded samples
get_excluded_samples_from_spreadsheet.R MSSNG_SSC_SPARK-WES1-WES2-PILOT_duplicates.xlsx

# Exclude samples from SPARK pilot data
exclude.py SPARK_WES_pilot/SPARK_WES_pilot.variant_data.tmp SPARK_pilot_samples_to_exclude.txt 1 > SPARK_WES_pilot/SPARK_WES_pilot.variant_data.tsv

# Concatenate all variant files together
concat_files_with_headers.sh 1 ASC/ASC.variant_data.hg38.tsv {MSSNG_ILMN,MSSNG_CG,SPARK_WES_pilot,SPARK_WES_1,SPARK_WES_2}/*.variant_data.tsv | perl -pe 's/^chr//g' > MSSNG+SPARK+ASC.variant_data.tmp1

# Exclude duplicate samples in MSSNG/SPARK (identified to be the same using PLINK)
exclude.py MSSNG+SPARK+ASC.variant_data.tmp1 MSSNG+SSC+SPARK_samples_to_exclude.txt 1 > MSSNG+SPARK+ASC.variant_data.tmp2

# Summarize the overlap in de novo variants among the samples
summarize_denovo_overlap.py MSSNG+SPARK+ASC.variant_data.tmp2 MSSNG+SPARK+ASC.variant_data.tmp2.overlap_summary MSSNG+SPARK+ASC.variant_data.tmp2.variants_with_duplicates MSSNG+SPARK+ASC.variant_data.tmp2.shared_denovo_counts MSSNG+SPARK+ASC.variant_data.tmp2.sample_variant_list

# Generate a list of samples in MSSNG and SPARK that share a de novo variant with a sample in ASC
# Only exclude samples where the same de novo is found AND the samples have the same sex
identify_samples_to_exclude_based_on_denovos.py MSSNG+SPARK+ASC.variant_data.tmp2.variants_with_duplicates $ALL_PLUS_METADATA_FILE ASC/ASC.sample_data.tsv > MSSNG_SPARK_overlap_with_ASC.details.txt
awk -F $'\t' '$8 == "yes"' MSSNG_SPARK_overlap_with_ASC.details.txt | cut -f 6 | sort | uniq > MSSNG_SPARK_overlap_with_ASC.txt
exclude.py MSSNG+SPARK+ASC.variant_data.tmp2 MSSNG_SPARK_overlap_with_ASC.txt 1 > MSSNG+SPARK+ASC.variant_data.tmp3

# Remove de novo variants with high recurrence counts
get_excluded_variants_from_spreadsheet.R recurrent_denovos_contributing_to_ASD_genes.annotated.xlsx
exclude.py MSSNG+SPARK+ASC.variant_data.tmp3 excluded_recurrent_variants.txt 0 > MSSNG+SPARK+ASC.variant_data.tsv

# Generate a list of PTVs that overlap between cases in ASC case-control data and de novo PTVs in MSSNG/SPARK.
tag_ASC_denovo_case_control_variants.py MSSNG+SPARK+ASC.variant_data.tsv ASC_variant_results.hg38.tsv > MSSNG+SPARK+ASC.variant_data.tagged.tsv
# Then extract variants that are 1) PTVs; 2) in MSSNG/SPARK; 3) Overlap with ASC case-control cohorts (DBS/SWE); 4) have at least one PTV in a case. Then take the unique variants, and extract the genes affected, the cohort (DBS or SWE), along with the number of cases in each
awk -F $'\t' '$5 == 2 && ($8 == "stop_gained" || $8 == "splice_acceptor_variant" || $8 == "splice_donor_variant" || $8 == "frameshift_variant") && ($11 == "DBS" || $11 == "SWE" || $11 == "DBS|SWE") && ($10 == "MSSNG_ILMN" || $10 == "MSSNG_CG" || $10 == "SPARK_WES_pilot" || $10 == "SPARK_WES_1" || $10 == "SPARK_WES_2") && $12 > 0' MSSNG+SPARK+ASC.variant_data.tagged.tsv | col_uniq.py 1 | cut -f 6,11,12 > cases_to_subtract.tsv
awk -F $'\t' '($5 == 2 && ($8 == "stop_gained" || $8 == "splice_acceptor_variant" || $8 == "splice_donor_variant" || $8 == "frameshift_variant") && ($11 == "DBS" || $11 == "SWE" || $11 == "DBS|SWE") && ($10 == "MSSNG_ILMN" || $10 == "MSSNG_CG" || $10 == "SPARK_WES_pilot" || $10 == "SPARK_WES_1" || $10 == "SPARK_WES_2") && $12 > 0) || NR == 1' MSSNG+SPARK+ASC.variant_data.tagged.tsv | cut -f 1,2,3,4,6,8,10 > supplementary_table_case-control_overlapping_variants.tsv


# Generate TADA counts table
generate_TADA_counts_table.py ASC/ASC.gene_data.tsv MSSNG+SPARK+ASC.variant_data.tsv cases_to_subtract.tsv MSSNG+SPARK+ASC.gene_data.tsv MSSNG+SPARK+ASC.gene_data.log

# Figure out sample size
ASC_sample_size="6430"
cat $SPARK_WES_PILOT_METADATA_FILE | awk -F $'\t' '$6 == "2" && $11 != "-" && $12 != "-" && $21 != "yes" && $23 != "yes"' | cut -f 2 | sort > SPARK_WES_pilot.txt
exclude.py SPARK_WES_pilot.txt <(cat SPARK_pilot_samples_to_exclude.txt MSSNG_SPARK_overlap_with_ASC.txt | sort | uniq) 0 > SPARK_WES_pilot.filtered.txt

cat $MSSNG_METADATA_FILE $SPARK_WES_1_METADATA_FILE $SPARK_WES_2_METADATA_FILE | awk -F $'\t' '$6 == "2" && $11 != "-" && $12 != "-" && $21 != "yes" && $23 != "yes"' | cut -f 2 | sort > MSSNG+SPARK.txt
exclude.py MSSNG+SPARK.txt <(cat MSSNG+SSC+SPARK_samples_to_exclude.txt MSSNG_SPARK_overlap_with_ASC.txt | sort | uniq) 0 > MSSNG+SPARK.filtered.txt

SPARK_WES_pilot_sample_size=$(cat SPARK_WES_pilot.filtered.txt | wc -l | tr -d "\n")
MSSNG_SPARK_sample_size=$(cat MSSNG+SPARK.filtered.txt | wc -l | tr -d "\n")
total_sample_size=$(($ASC_sample_size + $SPARK_WES_pilot_sample_size + $MSSNG_SPARK_sample_size))
original_ASC_cases="4811"
num_ASC_cases_to_subtract=$(cat cases_to_subtract.tsv | wc -l | tr -d " \n")
num_ASC_cases=$(($original_ASC_cases - $num_ASC_cases_to_subtract))
sample_size_file="sample_sizes.txt"
echo "ASC sample size is $ASC_sample_size" > $sample_size_file
echo "SPARK_WES_pilot sample size is $SPARK_WES_pilot_sample_size" >> $sample_size_file
echo "MSSNG+SPARK sample size is $MSSNG_SPARK_sample_size" >> $sample_size_file
echo "Total sample size is $total_sample_size" >> $sample_size_file
echo "Number of ASC cases is $num_ASC_cases" >> $sample_size_file


# Rearrange counts table and run TADA analysis
cat MSSNG+SPARK+ASC.gene_data.tsv | rearrange_columns.py "2,3,4,1,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27" > MSSNG+SPARK+ASC.gene_data.rearr.tsv
do_TADA.R MSSNG+SPARK+ASC.variant_data.tsv MSSNG+SPARK+ASC.gene_data.rearr.tsv $total_sample_size $num_ASC_cases > MSSNG+SPARK+ASC.TADA_results.tsv

# Generate list of genes with FDR < 0.1
tail -n +2 MSSNG+SPARK+ASC.TADA_results.tsv | awk -F $'\t' '$10 < 0.1' | cut -f 1 | sort > TADA+.MSSNG+SPARK+ASC.txt
ASC_gene_list="/hpf/largeprojects/tcagstor/users/btrost/papers/MSSNG/ASD_gene_lists/Satterstrom2020.txt"
comm -1 -3 TADA+.MSSNG+SPARK+ASC.txt $ASC_gene_list > TADA+.ASC_only.txt
comm -2 -3 TADA+.MSSNG+SPARK+ASC.txt $ASC_gene_list > TADA+.MSSNG+SPARK+ASC_only.txt
comm -1 -2 TADA+.MSSNG+SPARK+ASC.txt $ASC_gene_list > TADA+.both.txt


# Get only variants that contributed to those genes
include.py MSSNG+SPARK+ASC.gene_data.log TADA+.MSSNG+SPARK+ASC.txt 3 > variants_contributing_to_ASD_genes.tsv
tail -n +2 variants_contributing_to_ASD_genes.tsv | cut -f 1 | sort | uniq -c | niceuniq-c.sh | awk -F $'\t' '$2 > 1' | sort -rnk 2 > recurrent_denovos_contributing_to_ASD_genes.txt

# Get only genes from the counts file that had FDR < 0.1
include.py MSSNG+SPARK+ASC.gene_data.tsv TADA+.MSSNG+SPARK+ASC.txt 0 | cut -f 1,15-17,22-23 | awk -F $'\t' '{OFS="\t"; if (NR == 1) { cc_diff = "cc_diff.ptv"; total = "total" } else { cc_diff = $5 - $6; if (cc_diff < 0) { cc_diff = 0; } total = $2 + $3 + $4 + cc_diff}; print $1,$2,$3,$4,$5,$6,cc_diff,total}' > MSSNG+SPARK+ASC.gene_data.FDRlt0.1.tsv

mkdir -p recurrent_denovo_batch_files
for recurrent_denovo in $(cut -f 1 recurrent_denovos_contributing_to_ASD_genes.txt); do
	IFS=$'\n'
	denovo_dash=$(echo $recurrent_denovo | tr : "-")
	mkdir -p recurrent_denovo_batch_files/$denovo_dash
	for line in $(awk -F $'\t' -v recurrent_denovo=$recurrent_denovo '$1 == recurrent_denovo' MSSNG+SPARK+ASC.variant_data.tsv); do
		child=$(printf $line | cut -f 2)
		mom=$(printf $line | cut -f 3)
		dad=$(printf $line | cut -f 4)
		dataset=$(printf $line | cut -f 10)

		if [[ $dataset == "MSSNG_ILMN" ]]; then dataset="MSSNG"; fi

		chr=$(echo $recurrent_denovo | cut -f 1 -d :)
		pos=$(echo $recurrent_denovo | cut -f 2 -d :)

		printf "new\ngenome hg38\ngoto $chr:$pos\nload https://bounce2.tcag.ca/tcag/mirrored/$dataset/CRAM/$mom.cram\nload https://bounce2.tcag.ca/tcag/mirrored/$dataset/CRAM/$dad.cram\nload https://bounce2.tcag.ca/tcag/mirrored/$dataset/CRAM/$child.cram\n" > recurrent_denovo_batch_files/$denovo_dash/$child.bat
	done
done

# Annotate the TADA genes according to their EAGLE curations
TADA_gene_annotation.py /hpf/largeprojects/tcagstor/users/btrost/papers/MSSNG/ASD_gene_lists/docs/_ASDGeneCurationList_updated_2.18.2021.tsv TADA+.both.txt TADA+.ASC_only.txt TADA+.MSSNG+SPARK+ASC_only.txt > ASD_gene_annotation.tsv

SFARI_file="/hpf/largeprojects/tcagstor/users/btrost/papers/MSSNG/ASD_gene_lists/SFARI-Gene_genes_01-13-2021release_02-19-2021export.tsv"
ASC102="/hpf/largeprojects/tcagstor/users/btrost/papers/MSSNG/ASD_gene_lists/Satterstrom2020.txt"
module_map="module_map.data.tsv"
generate_ASD_gene_list_table.py MSSNG+SPARK+ASC.TADA_results.tsv MSSNG+SPARK+ASC.gene_data.rearr.tsv $SFARI_file $ASC102 EAGLE_genes.txt $module_map ../ASC+_no_MSSNG/only_with_MSSNG.txt > TADA+_ASD_gene_list.all.tsv # Generate file containing details for all genes
awk -F $'\t' 'NR == 1 || $2 < 0.1' TADA+_ASD_gene_list.all.tsv > TADA+_ASD_gene_list.tsv # Just get those that are FDR < 0.1
for gene in $(cat TADA+.ASC_only.txt); do cat ../ASC/ASC.TADA_results.tsv | awk -F $'\t' -v gene=$gene '$1 == gene'; done > TADA+.ASC_only.full.tsv # Just get those that dropped out of Satterstrom's list
generate_ASC_only_supplementary_table.py TADA+.ASC_only.txt Satterstrom_102_details.data.tsv TADA+_ASD_gene_list.all.tsv > ASC_only_supplementary_table.tsv # Generate the supplementary table containing the genes that dropped out of Satterstrom's list
