Multiomics Network Based Integration App

-Links genes, proteins, metablites to disease and pathways

-Creates a network and final association table as output

-Provides enrichment analysis of involved pathways metabolites and diseases

INPUT:
- WGS varinant annotated csv file- containing columns Gene and CADD
- Transcriptomics differentially expressed genes csv file with columns Gene, LogFC, P value and Padj value
- Proteomics intensity data obtained after protein abundunce analysis using maxquant with colums- Gene, Protein/ID and Intensity

Provided Input:
- cvd_genomics- WGS variant annotated csv- PRJNA264546
- cvd_transcriptomics- Bulk rna seq differential expression results- PRJNA394884
- cvd_proteomics- Protein intensity- 
