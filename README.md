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

Example Null W tables can be found [here](https://github.com/katlande/convergence_analysis/tree/master/Example_Data).

#### __________________

#### Section (2)
###### Required Files: window_rpkm_lib.txt, gene_window_library.txt, Tissue_Enrichment.py, condensed Null W tables, *optional: windows_to_genes.py*

**This section of the script:**
1. Matches hit windows to genes with known expression in 9 tissues
2. Checks for expression enrichment in each tissue
3. Corrects output data to an FDR of 0.05 and identifies tissues with enriched gene expression

In *Helianthus*, the baseline expression of each gene in 9 tissues is available on [sunflowergenome.org](https://sunflowergenome.org/jbrowse_current/?data=extdata%2Fbronze&loc=Ha1%3A116882..163932&tracks=Genes%2CTranscript%2CVariants%2CSeed%2CLigule%2COvary&highlight=).
Using this data, an rpkm value for each genome window could be produced. In this analysis, windows were consistent 
so the rpkm library was generated a single time with the python2.7 script [windows_to_genes.py](https://github.com/katlande/convergence_analysis/blob/master/Python2.7_Scripts/windows_to_genes.py).

The first half of section 2's R code pairs each window to an averaged RPKM value for 9 tissues using the [window_rpkm_lib.txt](https://github.com/katlande/convergence_analysis/blob/master/Genome_Data_Files/window_rpkm_lib.txt) file, and removes all hit 
windows that do not contain genes. These windows with RPKM values are then fed into the python2.7 script [Tissue_Enrichment.py](https://github.com/katlande/convergence_analysis/blob/master/Python2.7_Scripts/Tissue_Enrichment.py),
which compares the hit windows against all gene-containing windows in the *Helianthus* genome using the [gene_window_library.txt](https://github.com/katlande/convergence_analysis/blob/master/Genome_Data_Files/gene_window_library.txt) file and checks 
for tissue enrichment using permutations. Lastly, the output files of the python script are fed back into R, where 
the second half of section 2's R code runs a Benjamini-Hochberg FDR correction to 0.05. The resulting output files
can be used to determine if convergent genes in a species comparison are enriched or depleted in any tissues.

#### __________________

#### Section (3)
###### Required Files: Condensed Null W tables,windows_to_genes.txt, arabidopsis homolog file, *optional: windows_to_412.py*

**This section of the script** converts significant windows into *Helianthus* genes, pairs them with *Arabidopsis* homologs, and 
finds gene annotations for each significant gene.

This section uses the [windows_to_genes.txt](https://github.com/katlande/convergence_analysis/blob/master/Genome_Data_Files/windows_to_genes.txt) file to convert each window in a Null W table into a gene (or genes). 
As windows were consistent across analyses, this file was only generated a single time using the python2.7 script 
[windows_to_412.py](https://github.com/katlande/convergence_analysis/blob/master/Python2.7_Scripts/windows_to_412.py). Once significant windows are converted to genes, the script pairs each gene to its closest 
*Arabidopsis* homolog using the arabidopsis homolog file [Ha412_to_Arabidopsis.csv](https://github.com/katlande/convergence_analysis/blob/master/Genome_Data_Files/Ha412_to_Arabidopsis.csv) previously generated through BLAST. Annotations can 
then be obtained from the [TAIR bulk annotation tool](https://www.arabidopsis.org/tools/bulk/genes/index.jsp).

#### __________________

#### Section (4)
###### Required Files: Annotated gene files produced in section 3

**This section of the script** runs a gene ontology analysis on the significant genes using TopGo.
