rsync -avh /Users/sfoltz/Desktop/lab/zzz_other/resources/DEPO_final_20170206 DEPO_final_20170206.txt

rsync -avh 'sfoltz@linus6.gsc.wustl.edu:/gscmnt/gc2737/ding/Analysis/RNA-seq/fusion_paper/00_filtering_annotation/00_combined_fusion_file/Total_Fusions.tsv \
  /gscmnt/gc2737/ding/Analysis/RNA-seq/fusion_paper/00_filtering_annotation/00_combined_fusion_file/Hard_Filtered_Fusions.tsv \
  /gscmnt/gc2737/ding/Analysis/RNA-seq/fusion_paper/00_filtering_annotation/00_combined_fusion_file/Filtered_Fusions.tsv \
  /gscmnt/gc2737/ding/Analysis/RNA-seq/fusion_paper/00_filtering_annotation/00_combined_fusion_file/Fusions_with_many_partners.tsv \
  /gscmnt/gc2737/ding/Analysis/RNA-seq/fusion_paper/00_filtering_annotation/00_combined_fusion_file/Fusions_within_300kb.tsv \
  /gscmnt/gc2737/ding/Analysis/RNA-seq/fusion_paper/00_filtering_annotation/00_combined_fusion_file/Fusions_EFI.tsv \
  /gscmnt/gc2737/ding/Analysis/RNA-seq/fusion_paper/00_filtering_annotation/00_combined_fusion_file/Fusions_low_count.tsv \
  /gscmnt/gc2737/ding/Analysis/RNA-seq/fusion_paper/00_filtering_annotation/01_sample_set/sample_list.806.txt \
  /gscmnt/gc2737/ding/Analysis/RNA-seq/fusion_paper/00_filtering_annotation/01_sample_set/sample_list.with_file_names.txt \
  /gscmnt/gc2737/ding/Analysis/RNA-seq/fusion_paper/00_filtering_annotation/02_sv_annotation/discordant_reads/fusion_evidence_discordant_reads.100000.txt \
  /gscmnt/gc2737/ding/Analysis/RNA-seq/fusion_paper/00_filtering_annotation/02_sv_annotation/Filtered_Fusions_10000_delly_manta_20180817.txt \
  /gscmnt/gc2737/ding/Analysis/RNA-seq/fusion_paper/00_filtering_annotation/02_sv_annotation/Filtered_Fusions_100000_delly_manta_20180817.txt \
  /gscmnt/gc2737/ding/Analysis/RNA-seq/fusion_paper/00_filtering_annotation/03_cnv_annotation/all_output/combined_cnv_results.txt \
  /gscmnt/gc2737/ding/Analysis/RNA-seq/Gene_Expression/kallisto_tpm/out4/mmy_gene_tpm_table.tsv \
  /gscmnt/gc2737/ding/sample_info/clinical_data/clinical_data.20180813/Clinical_data.20180813.csv \
  /gscmnt/gc2737/ding/sample_info/SeqFISH.20180104.csv \
  /gscmnt/gc2737/ding/Analysis/RNA-seq/fusion_paper/00_filtering_annotation/00_combined_fusion_file/FilterDatabase/oncogene.tsv \
  /gscmnt/gc2737/ding/Analysis/RNA-seq/fusion_paper/00_filtering_annotation/00_combined_fusion_file/FilterDatabase/tsg.tsv \
  /gscmnt/gc2737/ding/Analysis/RNA-seq/fusion_paper/00_filtering_annotation/00_combined_fusion_file/FilterDatabase/mmy_known.tsv \
  /gscmnt/gc2737/ding/Analysis/RNA-seq/fusion_paper/00_filtering_annotation/00_combined_fusion_file/FilterDatabase/driver.tsv \
  /gscmnt/gc2737/ding/Analysis/RNA-seq/fusion_paper/00_filtering_annotation/00_combined_fusion_file/FilterDatabase/kinase.tsv ' .

head -n 55766 mmy_gene_tpm_table.tsv | cut -f3,4 > ensg_gene_list.tsv
