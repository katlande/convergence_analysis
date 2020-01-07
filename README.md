# Convergent Gene Analysis in *Helianthus* species
### (1) Scripts for converting GEA and GWAS hit windows into gene lists
### (2) Scripts for annotating gene lists, extracting gene ontology data from gene lists, and testing tissue enrichment of gene lists

This analysis begins with Null W tables of significant windows identified in a GWAS or GEA analysis.
The pipeline for a single Null W table will be outlined here, and is a summary of the [Convergence_Analysis_Full.R](https://github.com/katlande/convergence_analysis/blob/master/Convergence_Analysis_Full.R) script.

#### Set-up
Set up for this script requires only two inputs:
* The name of the test that produced the Null W tables (e.g. Baypass)
* The name of the variable that the initial analysis was looking at (e.g. Soil)

#### __________________

#### Section (1)
###### Required Files: Null W tables

**This section of the script** simply condenses Null W tables for each subvariable into a single Null_W table. 
> E.g., for each species comparison, the separate Null W tables for different soil variables would be combined into a single Null W table containing *every* significant window for soil variables.

#### __________________

#### Section (2)
###### Required Files: window rpkm library, gene rpkm library, permutation script, condensed Null W tables, *optional: rpkm library python script*

**This section of the script:**
1. Matches hit windows to genes with known expression in 9 tissues
2. Checks for expression enrichment in each tissue
3. Corrects output data to an FDR of 0.05 and identifies tissues with enriched gene expression

In *Helianthus*, the baseline expression of each gene in 9 tissues is available on sunflowergenome.org (LINK!!!).
Using this data, an rpkm value for each genome window could be produced. In this analysis, windows were consistent 
so the rpkm library was generated a single time with the python2.7 script (LINK!!!).

The first half of section 2's R code pairs each window to an averaged RPKM value for 9 tissues, and removes all hit 
windows that do not contain genes. These windows with RPKM values are then fed into the python2.7 script (LINK!!!),
which compares the hit windows against all gene-containing windows in the *Helianthus* genome (LINK!!!) and checks 
for tissue enrichment using permutations. Lastly, the output files of the python script are fed back into R, where 
the second half of section 2's R code runs a Benjamini-Hochberg FDR correction to 0.05. The resulting output files
can be used to determine if convergent genes in a species comparison are enriched or depleted in any tissues.

#### __________________

#### Section (3)
###### Required Files: Condensed Null W tables, windows to genes text file, arabidopsis homolog file, *optional: windows_to_412.py*

**This section of the script** converts significant windows into *Helianthus* genes, pairs them with *Arabidopsis* homologs, and 
finds gene annotations for each significant gene.

This section uses the windows_to_genes.txt (LINK) file to convert each window in a Null W table into a gene (or genes). 
As windows were consistent across analyses, this file was only generated a single time using the python2.7 script 
windows_to_412.py (LINK!!!). Once significant windows are converted to genes, the script pairs each gene to its closest 
*Arabidopsis* homolog using the arabidopsis homolog file (LINK!!!) previously generated through BLAST. Annotations can 
then be obtained from the TAIR bulk annotation tool (LINK!!!).

#### __________________

#### Section (4)
###### Required Files: Annotated gene files

**This section of the script** runs a gene ontology analysis on the significant genes using TopGo.







