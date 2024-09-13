# CSHL protocols: 
This folder contains the scripts and example data associated with the following publications:

**"Beyond Single Clones: High-Throughput Sequencing in Antibody Discovery"**

Ahmed S. Fahad1*, Matias F. Gutierrez-Gonzalez1*, Bharat Madan1, Brandon J. DeKosky1,2**

**"Clonal Variant Analysis of Antibody Engineering Libraries"**

Ahmed S. Fahad1*, Matias F. Gutierrez-Gonzalez1*, Bharat Madan1, Brandon J. DeKosky1,2**

**"Antibody Data Analysis from Diverse Immune Libraries"**

Ahmed S. Fahad1*, Matias F. Gutierrez-Gonzalez1*, Bharat Madan1, Brandon J. DeKosky1,2**

**"Clonal Lineage and Gene Diversity Analysis of Paired Antibody Heavy and Light Chains"**

Ahmed S. Fahad1*, Matias F. Gutierrez-Gonzalez1*, Bharat Madan1, Brandon J. DeKosky1,2**

1-The Ragon Institute of MGH, Harvard, and MIT, Cambridge, MA, USA 2-Department of Chemical Engineering, Massachusetts Institute of Technology, Cambridge, MA, USA

*Equal contribution **Corresponding author. E-mail: dekosky@mit.edu  

Requirements:

- A UNIX-like enviroment
- SLURM workload manager
- IgBlast
- Fastx_toolkit
- FLASH
- usearch v5
- perl5
- Python3
- Custom scripts in this repository
- Sample data provided in the manuscripts, or original HTS data in fastq format

# Protocol description

-**Protocol 1: Clonal Variant Analysis of Antibody Engineering Libraries**
  This protocol described a general method for tracking antibody clonal variants in display systems. In the example provided, the antibody vFP16.02 was diversified using site-saturation mutagenesis. 
  This methods generates all possible one amino acid mutations across the variable region on the heavy chain and light chain. The sample data provided consists on heavy chain variants sorted against two HIV Envelope proteins with different fusion peptide sequences. After quality control, sequences are aligned against the template antibody and the frequency of each variant is calculated. A custom script analyzed the frequency of each varaint across sorting rounds and calculates an enrichment ratio.

-**Protocol 2: Antibody Data Analysis from Diverse Immune Libraries** 
  This protocol is works in a similar way as Protocol 1. However, the input are immune libraries, which need to be aligned to germline V(D)J sequences.Next, clonal frequencies are calculated and the prevalence of each clone is calculates across sorting rounds. The prevalence is used to calculate an enrichment ratio. The sample data consist on immune libraries from donors recovering from ZIKV infection. Yeast libraries were sorted against ZIKV VLPs.

-**Protocol 3: Clonal Lineage and Gene Diversity Analysis of Paired Antibody Heavy and Light Chains**
This protocol provides the tools for paired VH:VL analysis of immune repertoires. With these analysis, the users can obtain information on gene usage, antibody isotype, and clonal lineage analysis
