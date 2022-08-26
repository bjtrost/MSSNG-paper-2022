module load R/R-4.0.2
while read i;do Rscript comppund_het_merger.R -n $i -p /gpfs/projects/tmedicine/KFakhroLAB/MSSNG/MSSNG_CG_exonic+splicing/ -o /gpfs/projects/tmedicine/KFakhroLAB/MSSNG/compound_het_experiment/ -m /gpfs/projects/tmedicine/KFakhroLAB/MSSNG/MSSNG_metadata.tsv&;done < proband.list


while read i;do echo "Rscript fam_merger_2.0_unaffected.R -n $i -p /gpfs/projects/tmedicine/KFakhroLAB/MSSNG/MSSNG_ILMN_filtered/ -o /gpfs/projects/tmedicine/KFakhroLAB/MSSNG/families_unaffected_ILMNandCG/ -m /gpfs/projects/tmedicine/KFakhroLAB/MSSNG/MSSNG_metadata_rev.tsv" | bsub -n 1 -e err.txt -o out.txt -P 0032;done < unaffected_sibling.list


while read i;do echo "Rscript fam_merger_2.0.R -n $i -p /gpfs/projects/tmedicine/KFakhroLAB/MSSNG/SSC_filtered/ -o /gpfs/projects/tmedicine/KFakhroLAB/MSSNG/families/ -m /gpfs/projects/tmedicine/KFakhroLAB/MSSNG/SSC_metadata.tsv" | bsub -n 1 -e err.txt -o out.txt -P 0032;done < ssc_proband.list

while read i;do echo "Rscript fam_merger_2.0.R -n $i -p /gpfs/projects/tmedicine/KFakhroLAB/MSSNG/MSSNG_ILMN_filtered/ -o /gpfs/projects/tmedicine/KFakhroLAB/MSSNG/families/ -m /gpfs/projects/tmedicine/KFakhroLAB/MSSNG/MSSNG_metadata_rev.tsv" | bsub -n 1 -e err.txt -o out.txt -P 0032;done < proband_illumina.list

while read i;do echo "Rscript fam_merger_2.0.R -n $i -p /gpfs/projects/tmedicine/KFakhroLAB/MSSNG/MSSNG_CG_exonic+splicing/ -o /gpfs/projects/tmedicine/KFakhroLAB/MSSNG/families/ -m /gpfs/projects/tmedicine/KFakhroLAB/MSSNG/MSSNG_metadata_rev.tsv" | bsub -n 1 -e err.txt -o out.txt -P 0032;done < proband.list


while read i;do echo "Rscript fam_merger_compound_het.R -n $i -p /gpfs/projects/tmedicine/KFakhroLAB/MSSNG/SSC_filtered/ -o /gpfs/projects/tmedicine/KFakhroLAB/MSSNG/families_compound_het/ -m /gpfs/projects/tmedicine/KFakhroLAB/MSSNG/SSC_metadata.tsv" | bsub -n 1 -e err.txt -o out.txt -P 0032;done < ssc_proband.list

while read i;do echo "Rscript fam_merger_compound_het.R -n $i -p /gpfs/projects/tmedicine/KFakhroLAB/MSSNG/MSSNG_ILMN_filtered/ -o /gpfs/projects/tmedicine/KFakhroLAB/MSSNG/families_compound_het/ -m /gpfs/projects/tmedicine/KFakhroLAB/MSSNG/MSSNG_metadata_rev.tsv" | bsub -n 1 -e err.txt -o out.txt -P 0032;done < proband_illumina.list

while read i;do echo "Rscript fam_merger_compound_het.R -n $i -p /gpfs/projects/tmedicine/KFakhroLAB/MSSNG/MSSNG_CG_exonic+splicing/ -o /gpfs/projects/tmedicine/KFakhroLAB/MSSNG/families_compound_het/ -m /gpfs/projects/tmedicine/KFakhroLAB/MSSNG/MSSNG_metadata_rev.tsv" | bsub -n 1 -e err.txt -o out.txt -P 0032;done < proband.list



