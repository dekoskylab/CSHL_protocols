# Clonal Variant Analysis of Antibody Engineering Libraries  

These scripts are helpful for the analysis of antibody clonal variants. By using techniques such as site-saturation mutagenesis, it is possible to generate all single amino acid mutations in an antibody variable region, which can be expressed in the surface of yeast, phage or other organisms. 

These libraries can be screened using antigens of interest. By enriching the library of particles that express antibody fragments with improved binding, it is possible to improve antibody function and evolve variants with enchance performance.

These scripts take as an input the output from IgBlast. First, the aligned nucleotide sequence is extracted and translated to amino acid. With this information, single or multiple mutants can be aligned to the antibody template used. Finally, all data from a single run can be compiled and antibody variant frequencies across selection rounds can be used to estimate antibody enrichment. 

## Steps

1. Run QC controls on FASTQ files. In our manuscript, we recommend using `FastQC` to profile the base call scores from Illumina MiSeq. Next, pair-end reads can be joined using `FLASh` and low quality reads can be removed using `fastq_quality_filter`
2. Align sequence data to germline sequences using `IgBlast`. This step is performed to obtain accurately the begining and the end of the antibody sequence
3. Run `Miseq_analysis_v3.py`.
4. Run `sort_blast` scripts
5. Run `aa_comparison_v3.py`
