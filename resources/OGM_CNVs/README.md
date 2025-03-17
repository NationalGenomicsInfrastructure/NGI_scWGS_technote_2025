# OGM CNV data from MM1S cell line

Main processing was done using the OGM software at Clinical genetics in Solna

Software: Bionano Solve 3.8 CNV annotation pipeline (Version: 810f67b223fdfa625256fe5ba2d9829164e52c51)

Files:
- `240222-01_-_Rare_Variant_Analysis_4_23_2024_10_41_34_Annotated_CNV.txt`: called CNVs

Header
```
# CNV pipeline calls
# Version: 810f67b223fdfa625256fe5ba2d9829164e52c51; installed on 2023-07-11 16:50:06
# annotation_source=Bionano Solve 3.8 CNV annotation pipeline
# command=run_cnv_annotation.py /home/bionano/access/local/jobs/5297/output/vap_params.txt
# case_cnv=/home/bionano/access/local/jobs/5297/output/data/alignmolvref/copynumber/cnv_calls_exp.txt case_cnv_stats=/home/bionano/access/local/jobs/5297/output/data/alignmolvref/copynumber/cnv_chr_stats.txt overlap_database=/home/bionano/access/share/jobs/5297/hg38.genes.bed cnv_dgv_database=/home/bionano/tools/pipeline/1.0/VariantAnnotation/1.0/config/data/homo_sapiens/dgv_hg38_gainLoss.txt cnv_genes_overlap_threshold=0 cnv_genes_overlap_method=OR cnv_genes_position_window=1000 cnv_dgv_overlap_threshold=30 cnv_dgv_overlap_method=OR cnv_dgv_position_window=1000 species=human_hg38 
```

- `cnv_rcmap_exp.txt`: Estimated copy number across genome

Header 
```
# CNV pipeline per-label coverage information
# Version: 810f67b223fdfa625256fe5ba2d9829164e52c51; installed on 2023-07-11 16:50:06
````

Renamed `240222-01_-_Rare_Variant_Analysis_4_23_2024_10_41_34_Annotated_CNV.txt` -> `cnv_calls.txt`