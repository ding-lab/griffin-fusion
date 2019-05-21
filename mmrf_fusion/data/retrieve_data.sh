rsync -avh /Users/sfoltz/Desktop/lab/zzz_other/resources/DEPO_final_20170206 DEPO_final_20170206.txt

rsync -avh 'sfoltz@virtual-workstation3.gsc.wustl.edu:/gscmnt/gc2737/ding/Analysis/RNA-seq/fusion_paper/00_filtering_annotation/00_combined_fusion_file/Total_Fusions.tsv \
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
  /gscmnt/gc2737/ding/Analysis/RNA-seq/fusion_paper/00_filtering_annotation/00_combined_fusion_file/FilterDatabase/kinase.tsv \
  /gscmnt/gc2737/ding/Analysis/RNA-seq/fusion_paper/03_kinase/Kinase_fusion_info.txt' .

head -n 55766 mmy_gene_tpm_table.tsv | cut -f3,4 > ensg_gene_list.tsv

wget https://www.cell.com/cms/10.1016/j.celrep.2018.03.050/attachment/41b9158b-96f6-457e-8837-551b91c40e67/mmc2.xlsx
mv mmc2.xlsx tcga_pancancer_fusions.xlsx

rsync -avh 'sfoltz@katmai.wusm.wustl.edu:/diskmnt/Projects/Users/sfoltz/griffin-fusion/mmrf_fusion/analysis/06_single_cell/results/27522_1/5000.discordant_reads.tsv' scRNA.discordant_reads.27522_1.tsv
rsync -avh 'sfoltz@katmai.wusm.wustl.edu:/diskmnt/Projects/Users/sfoltz/griffin-fusion/mmrf_fusion/analysis/06_single_cell/results/27522_1/chr4chr14.discovered_discordant_reads.tsv ' scRNA.discordant_reads.discover.27522_1.tsv
rsync -avh 'sfoltz@katmai.wusm.wustl.edu:/diskmnt/Projects/Users/sfoltz/griffin-fusion/mmrf_fusion/analysis/06_single_cell/results/56203_1/chr8chr14.discovered_discordant_reads.tsv ' scRNA.discordant_reads.discover.56203_1.tsv
rsync -avh 'sfoltz@katmai.wusm.wustl.edu:/diskmnt/Projects/Users/lyao/Tools/GRCh38_gencode_v29_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir/ref_annot.gtf.gene_spans' .
rsync -avh 'sfoltz@denali.wusm.wustl.edu:/diskmnt/Projects/Users/qgao/Priority/Gao_scRNA_cell_type/2.Prefiltered_celltype/backup_object_cell_type_in_sample_27522_1.rds ' scRNA.seurat_object.27522_1.rds
rsync -avh 'sfoltz@denali.wusm.wustl.edu:/diskmnt/Projects/Users/qgao/Priority/Gao_scRNA_cell_type/2.Prefiltered_celltype/backup_object_cell_type_in_sample_56203_1.rds ' scRNA.seurat_object.56203_1.rds

rsync -avh 'sfoltz@katmai.wusm.wustl.edu:/diskmnt/Projects/Users/sfoltz/griffin-fusion/mmrf_fusion/analysis/06_single_cell/scRNA.27522_1.bulk.Chimeric.out.junction' .
rsync -avh 'sfoltz@katmai.wusm.wustl.edu:/diskmnt/Projects/Users/lyao/MMY_RNAseq/STAR-Fusion/27522_1.ANNOTATS/star-fusion.fusion_predictions.tsv' scRNA.27522_1.bulk.star-fusion.fusion_predictions.tsv


#wxs_bm_data.withmutect.merged.maf.rc.caller.renamed.Bone_Marrow.tsv is a filtered,
#edited version of MGI:/gscmnt/gc2737/ding/Analysis/MMRF_MTP/SomaticWrapper/wxs_bm_data_merged/wxs_bm_data.withmutect.merged.maf.rc.caller
#see slack message with Hua Sun April 29, 2019
