# Matsumoto et al.2026_RNA-seq-analysis

## Description
This repository contains scripts used for RNA-seq analysis in Matsumoto et al. (2026).


## Dual RNA sequencing
RNA-seq data
>Sus scrofa (pig) swine intestinal epithelial cells: NCBI BioProject database (accession ID: PRJDB35408) approximately 4 Gb
>
Library preparation
>>NEBNext Poly(A) mRNA Magnetic Isolation Module (NEB E7490)(New England Biolabs, Ipswich, MA, USA)
>>NEBNext Ultra RNA Library Prep Kit for Illumina (E7530) (New England Biolabs, Ipswich, MA, USA)
>>
>Lactiplantibacillus plantarum JCM 1149T: NCBI BioProject database (accession ID: PRJDB35409) approximately 1 Gb
>
Library preparation
>>Ribo-Zero Plus rRNA Depletion Kit (Illumina, San Diego, CA, USA)
>>NEBNext Ultra II Directional RNA Library Prep Kit (New England Biolabs)
>

Sequencing 
>Illumina NovaSeq 6000 (paired-end 150-bp reads)


Reference genomes
>Sscrofa11.1 (GCA_000003025.6)
>>https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000003025.6/
>>
>L.plantarum subsp. plantarum ATCC 14917 = JCM 1149 = CGMCC 1.2437 (GCA_000143745.1)
>>https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000143745.1/


### Workflow
1. Adapter trimming and quality filtering (Python 3.9, Trimmomatic v0.38)
2. Read mapping (Python 3.9, HISAT2 v2.1.0)
3. Read counting (Python 3.9, featureCounts v1.6.3)
4. Differential expression analysis (R v4.2.3, DESeq2 v1.36.0)
5. Post-processing and visualization (Python 3.10, scikit-learn v1.0.2, seaborn v0.13.2, pandas v2.2.3, matplotlib v3.10.0, numpy v1.26.4)

### Statistical criteria
- Adjusted p-value < 0.05 (Benjamini–Hochberg)
