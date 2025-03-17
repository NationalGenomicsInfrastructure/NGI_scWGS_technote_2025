# WGS on single cells with BioScryb and Qiagen kits


BioScryb kit: ResolveDNA Library Preparation Kit
Qiagen kit: QIAseq® Single Cell DNA Library Kit with Unique Dual Indexes

Samples
- Qiagen:
  - 10 st single cell (MM1S)
  - 1 st mini-bulk, 10 cells (MM1S)
  - 1 st positive control DNA (mouse?)
- BioScryb:
  - 10 st single cell (MM1S)
  - 1 st mini-bulk, 10 cells (MM1S)
  - 1 st positive control DNA (mouse?)

## Directory organization

- `analysis`: Analysis results and run files.
  - `copykit`: CNV analysis with copykit
  - `downsample_40M`: Analysis on samples downsampled to 40M reads.
    - `contaminats`: Contamination analysis
    - `downstream`: Downstream QC
    - `sarek_GRCh38`: Genome alignment and QC using `nf-core/sarek`
  - `downsample_qc`: Analysis on samples downsampled to 40M reads.
    - `contaminats`: Contamination analysis
    - `downstream_qc`: Downstream QC
    - `sarek_GRCh38`: Genome alignment and QC using `nf-core/sarek`
  - `qiagen_repli_g_full`: Analysis on Qiagen samples without downsampling.
    - `downstream`: Downstream QC
    - `sarek`: Genome alignment and QC using `nf-core/sarek`
- `data`: Sequencing data generated for project.
  - `P33410`: Raw FASTQS
  - `P33410_downsampled` Downsampled FASTQs and run files for this generation
- `reports`: HTML reports related to different steps of analysis
- `resources`: Resource files related to analysis
  - `OGM_CNVs`: CNV calls and coverage files for Optical genome mapping (OGM) data
  - `PacBioHiFi_pr_023_003_MM1S_hificnv`: CNV calls and coverage files from PacBio HiFi data
- `scripts`: Python scripts related to analysis

## Links

Genstat project: [G.Applications_24_11, P33410](https://genomics-status.scilifelab.se/project/P33410)
Google drive: [GA_24_01](https://drive.google.com/drive/folders/1YCGlXWzlOfdedT8PfMR7-UReDPoCYJzM?usp=drive_link)
