#This script encompasses all the R analyses done for tissue enrichment, gene annotation, and gene ontology for the GEA and GWAS convergent genes.

#Section 1: Initial data set up from raw Null_W tables
#Section 2: Tissue Enrichment Analysis
#Section 3: Gene Annotation
#Section 4: Gene Ontology 

rm(list=ls())
library(readxl)
library(tidyverse)
library(rlist)
library(tidyr)
library(stringr)
library(reshape2)

#To start, specify which test and condition you're looking at:
Test <- "Spearman" #Spearman, Baypass, GWAS
Condition <- "Soil" #Climate, Soil, Phenotype

#########################################################################################
################################ Section 1:  Data Set Up ################################
#########################################################################################

#We start from the raw Null_W tables output by the analyses.
#I set my working directory automatically to work with my filing system. Other users may have to do this manually or alter the working directory names. 
raw_data_dir <- paste("~/Desktop/Convergence_Paper/",Test,"/Raw_Data/",Condition,sep="")
setwd(raw_data_dir)
raw_data <- list.files(path=raw_data_dir)#whenever using the list.file() command, make sure the path contains ONLY the fils/directories you want included in the list it create (this command makes a list of all the species comparisons, as the path contains 6 folders names for each species comparison).

#The first step is to take all the hits for individual variables and condense them into a single file. If  you already have all the hits for a comparison condensed into a single file, you can skip to Section 1.5.

#We will do this in a loop for each file:
r = 1
while(r<length(raw_data)+1){
  temp_dir <- paste(raw_data_dir,"/",raw_data[r],sep="")#each comparison is a subdirectory within the raw data directory. Here I am changing the directory to each subdirectory to pull out files. Other users may have to alter this part of the code or change their filing system to make it run.
  setwd(temp_dir)
  temp_file_list <- list.files(path=temp_dir)#pull out all the files in the current subdirectory an make them into a list
  temp_file_list_noNA <- c("")
  for (file in temp_file_list){#check all the files in the temporary list, some may be empty and this will cause errors later in the script
    if (file.size(file) == 0) next#if any files are empty, they will be ignored
    temp_file_list_noNA <- list.append(temp_file_list_noNA, file)#while files with content are added to the updated file list
  }
  temp_file_list_noNA <- temp_file_list_noNA[c(2:length(temp_file_list_noNA))]#remove the first index from the true file list as it was just an empty placeholder
  t = 1
  while(t < length(temp_file_list_noNA)+1){#we'll us a nested loop her to bind all the files in the subdirectory into a single dataframe
    temp <- read.table(temp_file_list_noNA[t], header = F)#read each file as a temporary data frame
    temp <- temp[c(1,3,4)]#select the window, the p value, and the variable columns only
    colnams <- c("Variable", "Window", "p_Val")
    colnames(temp) <- colnams
    temp$Test <- Test#add a test column
    temp$Comparison <- raw_data[r]#add a comparison column. This part only works because the name of each of my subdirectories is the name of the comparison. If this is not the case for other users, you may have to alter your filing system or alter the code here.
    if(t == 1){
      all_data <- temp#the first dataframe is the main dataframe
    } else {
      all_data <- rbind(all_data,temp)#and we'll row bind all the following dataframes onto it
    }
    t = t + 1#repeat the nested loop for all files in the subdirectory
  }
  all_data <- all_data %>% drop_na()#remove NA values 
  temp_output_name <- paste(Test, "_", Condition, "_", raw_data[r], ".txt", sep="")#make a temporary output name for the file
  setwd(raw_data_dir)#set to main directory for file output
  write.table(all_data, temp_output_name, row.names = F, quote = F, sep = "\t")#write a text file of the full data for each comparison. We will use it in the next steps.
  r = r + 1#repeat the main loop for all subdirectories in the raw_data directory
}#You should now have one file per comparison in you main raw_data directory.

#### SECTION 1.5 ####
#This section is ONLY USED if you are starting with Null W tables that contain every variable, rather than individual Null W tables for each variable. If you ran all of section 1.0, skip to section 2.

#raw_data_dir <- paste("~/Desktop/Convergence_Paper/",Test,"/Raw_Data/",Condition,sep="")
#setwd(raw_data_dir)
#files <- list.files(path=raw_data_dir, pattern = ".txt")
#The next two lines change depending on your file names. I am pull the name of the comparison out of the file names. My input files are in the following format:
#"Test_Condition_Comparison.txt"
#remove1 and remove2 are the strings I have to remove from my file name to be left with ONLY the comparison name:
#remove1 <- paste(Test,"_",Condition,"_",sep="")
#remove2 <- paste(".txt")

#t = 1
#while(t < length(files)+1){#we'll us a nested loop her to bind all the files in the subdirectory into a single dataframe
#  temp <- read.table(files[t], header = T)#read each file as a temporary data frame
#  temp <- temp[c(1,3,4)]#select the window, the p value, and the variable columns only
#  colnams <- c("Variable", "Window", "p_Val")
#  colnames(temp) <- colnams
#  temp$Test <- Test#add a test column
#  temp$Comparison <- files[t]#the comparison column is the file name.
#  temp$Comparison<-gsub(remove1,"",as.character(temp$Comparison),n)#remove the strings remove1 and remove2 from the cells in the comparison column
#  temp$Comparison<-gsub(remove2,"",as.character(temp$Comparison),n)
#  t = t + 1#repeat the nested loop for all files in the subdirectory
#  out_name <- paste(Test, "_", Condition, "_", temp$Comparison[1], ".txt", sep="")#chose an output name
#  write.table(temp, out_name, row.names = F, quote = F, sep = "\t")#write a table to use in section 2 for each species comparison
#}


#########################################################################################
############################# Section 2:  Tissue Enrichment #############################
#########################################################################################

#This section will set up the data for the Tissue Enrichment analysis, instruct on how to run the enrichment python script, and do the back-end analysis of the output data.

#Set up: For a file input into the python script, we just need a list of all the hit windows that CONTAIN GENES. Using a different python script (windows_to_genes.py), I previously generated a library of all the windows used in the GWAS and GEA analyses. The library file (window_rpkm_lib.txt) contains rpkm data for each window (mined from sunflowergenome.org) and an identifier as to whether or not each window contains gene region(s).

genome_files <- "~/Desktop/Convergence_Paper/Genome_Data"#Set the name of the directory that contains all my genome data to an object for later use.
setwd(genome_files)
windowlibrary <- read.table("window_rpkm_lib.txt", header = TRUE)#load in the library
gene_lib <- subset(windowlibrary, is_gene. == "['YES']" )#pull out only the windows that contain gene regions
gene_windows <- gene_lib[c(1)]#Just a list of all the windows with genes

#raw_data_dir <- paste("~/Desktop/Convergence_Paper/",Test,"/Raw_Data/",Condition,sep="")#If you skipped section 1, run this line. Otherwise, ignore.
setwd(raw_data_dir)
enrichment_input <- list.files(path=raw_data_dir, pattern = ".txt")#input ONLY the text files in your raw data directory.
gene_windows_dir <- paste("~/Desktop/Convergence_Paper/",Test,"/Significant_Gene_Windows/",Condition,sep="")#make a working directory for your output files

#The input for the python script requires a text file for each comparison that is JUST the gene windows. We will make them in a loop:
for (file in enrichment_input){#for each comparison, make a significant hit window file.
  setwd(raw_data_dir)
  temp_sample <- read.table(file, header = T)
  colnames(temp_sample)[2] <- "window"
  temp_sample <- temp_sample[2]#isolate the significant windows
  significant_genes <- merge(gene_windows,temp_sample, by = "window", all = F)#merge the gene windows and the hit windows files to keep only hit windows that contain genes
  filename <- str_remove(file,".txt")
  out_name <- paste(filename,"gene_windows.txt",sep="")
  setwd(gene_windows_dir)#set wd to the path you want your output files in
  write.table(significant_genes, out_name, col.names = F, row.names = F, quote = F)
}

#Now, we can take these files and input them into the enrichment python script, Tissue_Enrichment.py. The script must be run once for EACH input file.

#This script runs with PYTHON 2!!!!!! Using Python 3 will NOT WORK!!!!!! Additionally, the NumPy package is required for this script to run.

#Run Tissue_Enrichment.py. This script runs permutations tests on the tissue rpkm values, using full genome data as the null distribution and your sample hits as the sample distribution.
#(1) The code will first ask if you're running on a server or not. If you type y, the only difference is that it will prompt you to enter file names manually rather than by giving you a dialog box. If you do not have the Tkinter package installed on your local device, use server mode.
#(2) The code will ask if you want to run low or high memory mode. If you choose yes, the code will take twice as many random samples to create a null mean during each permutation. For large libraries and many permutations, this can take a considerable amount of computational time. Essentially, it will shrink the standard deviation of your null distribution and make fine enrichment more visible. From analyzing the low memory output, you can decide if a high memory run is necessary or not.
#(3) You will be asked for a gene library. This file contains ALL the windows in the genome that contain gene regions, and the average rpkm of all the genes in each window in nine tissues. It is the "window_rpkm_lib.txt" file used in the script above, but with the non-gene windows removed.
#(4) You will then be asked for you list of significant genes. This is the list of gene windows we just generated. Each time you run the script, this is the only input you will change.
#(5) You will be asked to enter an output file name. To make the rest of the script run smoothly, I recommend: "Comparison_Test_Condition_10k.txt"
#(6) You will be asked for the number of permutations you want to run. Between 10,000 and 100,000 is the broadly recommeneded range for permutation tests you want to publish. We used 10,000 in this analysis.

#Run the script for each input file, then once all have been run we can do an FDR correction.

enrichment_output_dir <- paste("~/Desktop/Convergence_Paper/",Test,"/Tissue_Enrichment/",Condition,sep="")#new working directory
setwd(enrichment_output_dir)
enrichment_files <- list.files(path=enrichment_output_dir, include.dirs = F, pattern = ".txt")#new file list

r = 1
while(r < length(enrichment_files)+1){#read in all the enrichment files
  temp <- read.table(enrichment_files[r], header = T)
  c <- enrichment_files[r]
  remove_temp <- paste("_",Test,"_",Condition,"_10k.txt",sep="")#alter this based on what you chose for your output name, unless you used the recommeneded output format.
  c <- str_remove(c, remove_temp)
  temp$comparison <- c#set the comparison name as a column
  if(r == 1){enrichment_data <- temp} else {enrichment_data <- rbind(enrichment_data, temp)}#rbind it all together
  r = r + 1
}
enrichment_data$condition <- Condition
enrichment_data$test <- Test#set columns for condition and test as well

#Now that the data is all together, we'll do a Benjamini-Hochberg FDR correction. We'll correct to 0.1, as well as 0.05.

t = 1
while (t < (nrow(enrichment_data)+1)) {#for each row of the dataframe, we'll check if the sample gene expression is enriched or depleted. If the enrichment_p is more significant than the depleteion_p, it's enriched, and vice versa. The enrichment or depletion p (whichever was more significant) is then set as the "true_p."
  if (enrichment_data$Enrichment_p[t] < enrichment_data$Depletion_p[t]){enrichment_data$true_p[t] <- enrichment_data$Enrichment_p[t]
  enrichment_data$enriched_or_depleted[t] <- "ENRICHED"}
  if (enrichment_data$Enrichment_p[t] > enrichment_data$Depletion_p[t]){enrichment_data$true_p[t] <- enrichment_data$Depletion_p[t]
  enrichment_data$enriched_or_depleted[t] <- "DEPLETED"}
  t = t+1
}

BH_test_10k <- enrichment_data
BH_test_10k <- BH_test_10k[order(BH_test_10k$true_p),]#order the dataframe by p-value

t = 1
while (t < (nrow(BH_test_10k)+1)) {#add a column that ranks all the results by true p-value, from most to least significant
  BH_test_10k$rank[t] <- t
  t = t+1
}

t = 1#FDR correction to 0.1 FDR:
while (t < (nrow(BH_test_10k)+1)){BH_test_10k$p_0.1[t] <- ((BH_test_10k$rank[t]/nrow(BH_test_10k))*0.1)
if (BH_test_10k$p_0.1[t] > BH_test_10k$true_p[t]){BH_test_10k$FDR_0.1[t] <- "PASS"}
if (BH_test_10k$p_0.1[t] < BH_test_10k$true_p[t]){BH_test_10k$FDR_0.1[t] <- "FAIL"}
t= t+1 }

t = 1#FDR correction to 0.05 FDR:
while (t < (nrow(BH_test_10k)+1)){BH_test_10k$p_0.05[t] <- ((BH_test_10k$rank[t]/nrow(BH_test_10k))*0.05)
if (BH_test_10k$p_0.05[t] > BH_test_10k$true_p[t]){BH_test_10k$FDR_0.05[t] <- "PASS"}
if (BH_test_10k$p_0.05[t] < BH_test_10k$true_p[t]){BH_test_10k$FDR_0.05[t] <- "FAIL"}
t= t+1 }

BH_test_10k <- BH_test_10k[c(-2,-3,-13)]#Remove columns that are unnecessary
BH_test_10k$Organ <- str_remove(BH_test_10k$Organ, ".rpkm")#remove .rpkm from the organ names
enriched_df <- subset(BH_test_10k, BH_test_10k$enriched_or_depleted == "ENRICHED")#create two files -- enriched tissues and depleted tissues
depleted_df <- subset(BH_test_10k, BH_test_10k$enriched_or_depleted == "DEPLETED")

enrichment_results_dir <- paste("~/Desktop/Convergence_Paper/",Test,"/Tissue_Enrichment/Results",sep="")
setwd(enrichment_results_dir)

#make file names for the outputs:
Enrichment_file_Output <- paste(Test,"_",Condition,"_Tissue_Enrichment.txt",sep="")
Depletion_file_Output <- paste(Test,"_",Condition,"_Tissue_Depletion.txt",sep="")

#write output files:
write.table(enriched_df, Enrichment_file_Output, row.names = F, col.names = T, sep = "\t", quote = F)
write.table(depleted_df, Depletion_file_Output, row.names = F, col.names = T, sep = "\t", quote = F)

#Tissue enrichment analysis is now done! Yay!

#########################################################################################
############################## Section 3:  Gene Annotation ##############################
#########################################################################################

#To get gene annotations, as well as gene ontologies, we are going to use Arabidopsis homologs for our genes. First, we have to figure out which sunflower genes were actually significant. To do this, we'll use the windows_to_genes.txt file. This file was generated from the windows_to_412.py python script, it works similarly to the windows_to_genes.py script, but pulls out gene names rather than window RPKMs.

setwd(genome_files)
windows_to_genes <- read.table("windows_to_genes.txt", header = T, fill = T)#you will note that this file contains a column for EACH gene found within a window. We will reformat this file to make it easier to use in a moment.

#If you skipped part 2, run this. Otherwise, you can ignore it:
#raw_data_dir <- paste("~/Desktop/Convergence_Paper/",Test,"/Raw_Data/",Condition,sep="")
#enrichment_input <- list.files(path=raw_data_dir, pattern = ".txt")

#We'll process each of the comparisons separately in a loop. Here we'll pull out hit genes, their p-values, and which variables they are significant for:
y = 1
while(y < length(enrichment_input)+1){
  remove <- paste(Test,"_",Condition,"_",sep="")
  comparison <- str_remove(enrichment_input[y],remove)
  comparison <- str_remove(comparison,".txt")#get the comparison name from the file name.
  
  setwd(raw_data_dir)
  hits <- read.table(enrichment_input[y], header = T)#read each file
  hits_imp <- hits[c(1:3)]#remove the variable, the window name, and the p-value
  colnams <- c("phenotype","window", "p_val")
  colnames(hits_imp) <- colnams
  
  hit_windows <- merge(hits_imp,windows_to_genes, by = "window", all = F)#merge the full window to gene library with the hits
  
  hit_windows <- hit_windows[order(hit_windows$is_gene., decreasing = T),]#order the dataframe so that the highest number of genes per window in the file is in the first row
  gene_num<-hit_windows$is_gene.[1]#save the max number of genes per window in the hits file as an object
  
  
  #Below: reformat the file so that gene is it's own column, rather than multiple columns. The window column will no longer contain unique values:
  q = 6
  while( q < (ncol(hit_windows)+1)){#all gene columns are just called "gene"
    colnames(hit_windows)[q] <- "gene"
    q = q+1
  }
  
  q = 6
  r = 1
  while( r < (gene_num+1)){
    temp<-hit_windows[c(1:3,q)]#make a separate temporary object for each gene column
    temp<-na.omit(temp)
    if(r == 1){gene_list <- temp} else {gene_list <- rbind(gene_list, temp)}#then combine them all
    r = r + 1
    q = q + 1
  }
  
  #You should now have an object with a single gene column.
  gene_list2 <- gene_list[c(3,4)]#remove non-numerical data
  gene_list2 <- aggregate(gene_list2, by = gene_list2[c(2)], FUN = mean)#The gene column may not be unique. To make it unique, aggregate the gene values together by taking the mean of the p-value for non-unique genes. You'll probably get some warnings because the gene column is non-numeric. Don't worry, it still worked, R just created an extra column full of NAs that we can delete:
  gene_list2 <- gene_list2[c(-3)]
  gene_phenotpes <- gene_list[c(2,4)]#pull out the gene names and the significant variable for those genes
  
  count.duplicates <- function(DF){#Define a function that counts the number of times an object in a column is duplicated
    x <- do.call('paste', c(DF, sep = '\r'))
    ox <- order(x)
    rl <- rle(x[ox])
    cbind(DF[ox[cumsum(rl$lengths)],,drop=FALSE],count = rl$lengths)
  }
  
  temp<-gene_phenotpes$gene#pull out the genes from the phenotype table
  temp<-as.data.frame(temp)#make it a temporary dataframe
  duplicate_table<- count.duplicates(temp)#count the number of variables each gene is significant for
  duplicate_table <- duplicate_table[order(duplicate_table$count, decreasing = T),]#order the table so the highest value is in the top row
  pheno_list <- unique((gene_phenotpes$phenotype))#make a list of all the variables that had hits in this analysis
  pheno_list <- as.list(levels(pheno_list))
  
  #Make columns for all possible phenotype hits:
  q = 1
  while( q < (length(pheno_list)+1)){#for each object in the list of variables
    gene_phenotpes$temp <- 0#make a new column and fill it with 0s
    colnames(gene_phenotpes)[2+q] <- pheno_list[q]#and make the object in the variable list the name of the column
    q=q+1#repeat for all variables
  }
  
  q = 1#Now run a loop that replaces the 0 with a 1 if a gene has a hit for that variable:
  while(q < nrow(gene_phenotpes)+1){
    col <- match(gene_phenotpes$phenotype[q],pheno_list)
    temp_row <- gene_phenotpes[c(q),]
    gene_phenotpes <- gene_phenotpes[c(-q),]
    tru_col <- (col+2)
    temp_row[tru_col] <- 1
    gene_phenotpes <- rbind(temp_row,gene_phenotpes)
    q = q+1
  }
  
  gene_phenotpes <- gene_phenotpes[c(-1)]#remove the phenotype column, it's no longer necessary
  pheno_final <- aggregate(gene_phenotpes, by = gene_phenotpes[c(1)], FUN = mean)#Now we aggregate the data frame together so that the gene column is only unique values. You'll get warnings again, which you can ignore for the same reason. For each gene, every phenotype column with a non-0 number in it had a hit for that phenotype!
  pheno_final <- pheno_final[c(-2)]
  
  pheno_final <- pheno_final %>% column_to_rownames("gene")
  pheno_final[pheno_final > 0] <- "HIT"#Replace all non-0s with a "HIT" classifier
  pheno_final[pheno_final == 0] <- "NO HIT"#Replace all 0s with a "NO HIT" classifier
  pheno_final <- pheno_final %>% as.data.frame()
  pheno_final <- pheno_final %>% rownames_to_column("gene")
  
  output_file <- cbind(gene_list2, pheno_final)#column bind the gene pvalue and phenotype hit information into one dataframe
  output_file <- output_file[c(-3)]#remove extraneous columns
  
  Annotation_output <- paste("~/Desktop/Convergence_Paper/",Test,"/Annotations",sep="")
  setwd(Annotation_output)
  out_name <- paste(comparison,"_",Test,"_",Condition,"_hits.txt",sep="")
  write.table(output_file, out_name, row.names = F, col.names = T, sep = "\t", quote = F)
  #gene_list2 is now a 2 column df of gene name && p value
  y = y +1
}

#Now we actually have to annotate this data. The first step is to find Arabidopsis homologs. Thankfully, a list of TAIR homologs for Ha412 genes already exists! We'll read it in and find matches:

setwd(genome_files)
Imp_Arab <- read.csv("Ha412_to_Arabidopsis.csv") #arabidopsis homolgs to 412 genes
colnames(Imp_Arab)[2] <- "gene"
Imp_Arab <- Imp_Arab[c(2:4)]

Annotation_dir <- paste("~/Desktop/Convergence_Paper/",Test,"/Annotations",sep="")
setwd(Annotation_dir)
annot_files <- list.files(path=Annotation_dir, include.dirs = F, pattern = Condition)

#NOTE: This part of the script is unfortunately very manual. Make a directory that is easy to acces such as your desktop the "easy directory": 
easy_dir <- "~/Desktop"

#Make a list of hit gene homologs for each comparison:
for(file in annot_files){
  remove <- paste("_",Test,"_",Condition,"_hits.txt",sep="")
  Comparison <- str_remove(file,remove)
  input_name <- file
  
  #Match Ha412 genes to their TAIR homolog:
  setwd(Annotation_dir)
  Genes <- read.table(input_name, sep = "\t", header = T) #list of significant genes
  All_412 <- merge(Imp_Arab,Genes) #label significant genes with homolog names
  if(file == annot_files[1]){annot_file_list <- list(All_412)#save each of these files with a unique name, for later
  } else {annot_file_list <- list.append(annot_file_list, All_412)}
  colnames(All_412)[3] <- "homology_p" #rename/clean file
  colnames(All_412)[4] <- paste(Test,"_p") #rename/clean file
  
  #Output the homologs only as a list:
  homologlist <- All_412$homolog
  setwd(easy_dir)
  homo_name <- paste(Comparison,"_convergence_TAIR.txt",sep="")
  write.table(homologlist, homo_name, row.names = F, col.names = F, quote = F)
}

#To annotate, I use an online tool from TAIR. There should be one homolog list for each comparison now in your easy-to-access directory. Take these files, and input them into the TAIR bulk annotation service (done at https://www.arabidopsis.org/tools/bulk/genes/index.jsp).
#Take the output and save it as "comparison_input.csv"
#I find the best way to get this to work is to copy the annotations into a text file, then from the text file copy them into excel and save as a csv in your "easy_dir". This part is, unfortunately, a bit of a pain.

#Input the csv's back into R to produce final annotation tables:
setwd(easy_dir)
annotated_files <- list.files(path=easy_dir, include.dirs = F, pattern = "_input.csv")
annotated_dir <- paste(Annotation_dir,"/Annotated/",Condition,sep="")


f = 1
while(f < length(annotated_files)+1){#make annotation files in a loop:
  filename <- str_remove(annotated_files[f],"_input.csv")
  setwd(easy_dir)
  homologs <- read.csv(annotated_files[f]) #put annotated homolgs in
  names(homologs)[names(homologs) == "Locus.Identifier"] <- "homolog" #rename homolog column 
  
  complete <- merge(homologs,annot_file_list[f])
  setwd(annotated_dir)
  output_name <- paste(filename,"_",Test,"_",Condition,"_Annotations.txt",sep="")
  write.table(complete, output_name, row.names = F, quote = F, sep = "\t")
  f = f +1
}

#If the annotation files are in your annotated directory, you can delete all the files from your desktop. Annotations are complete! Huzzah!


#########################################################################################
############################### Section 4:  Gene Ontology ###############################
#########################################################################################

library(topGO)

#Set names for your GO working directories. Input data is the annotation files, we already have a directory for this from section 3.
GO_out_dir <- paste("~/Desktop/Convergence_Paper/",Test,"/GO/",Condition,sep="")#output directory for GO files
GO_PDF_dir <- paste("~/Desktop/Convergence_Paper/",Test,"/GO/GO_PDFS",sep="")#output directory for node diagrams (optional, makes looking at GO results easier)
#Input files for GO:
GO_inputs <- list.files(path=annotated_dir)

#Now we'll run GO biological process ontology in a loop:
for(file in GO_inputs){
  setwd(annotated_dir)#Working in the directory of your annotations
  data <- read.delim(file, sep = "\t", header = T, fill = T)#input each file
  
  remove <- paste("_",Test,"_",Condition,"_Annotations.txt",sep="")
  Comparison <- str_remove(file,remove)
  
  require(topGO)
  require(org.At.tair.db)#load the TAIR gene ontology data
  
  data <- data[c(1,9)]#from the annotation data, we only need the name of the arabidopsis homolog and it's p-value from the test
  data <- data %>% as.data.frame()
  
  #Create a genelist object:
  data_geneList <- data$p_val#list of gene p-values
  names(data_geneList) <- data$homolog#homolog names
  
  #We will conduct two gene ontologies per comparison. One for biological process GO terms, and one for molecular function GO terms. We'll start with BP.
  data_Selection <- function(allScore){ return(allScore < 0.05)} # function that returns TRUE/FALSE for p-values<0.05. All the genes are already significant anyways so this step doesn't actually do anything except put the data in a usable format.
  data_All <- annFUN.org(whichOnto="BP", feasibleGenes=NULL, mapping="org.At.tair.db", ID="GeneName")#specify your ontology. We want biological process mapped to the TAIR database.
  data_GOdata <- new("topGOdata",#make a topGO object using the inputs from above
                     ontology = "BP",
                     allGenes = data_geneList,
                     geneSelectionFun = data_Selection,
                     annot = annFUN.org, mapping = "org.At.tair.db")
  #We conduct two statistical tests for enrichment, fisher's exact test (ignores p-values and just looks at enrichment) and the Kolomogorov-Smirnov test (uses the p-value of genes to rank GO terms and then looks for rank deviations from the full population of TAIR genes).
  resultFisher <- runTest(data_GOdata, algorithm = "classic", statistic = "fisher")
  resultKS <- runTest(data_GOdata, algorithm = "classic", statistic = "ks")
  #Make an ontology table:
  if(length(data_GOdata@graph@nodes) > 59){#if there are more than 60 nodes, include the top 60 only.
    data_tab_BP <- GenTable(data_GOdata, classicFisher = resultFisher, classicKS = resultKS, orderBy = "classicKS", ranksOf = "classicFisher", topNodes = 60)
  } else {data_tab_BP <- GenTable(data_GOdata, classicFisher = resultFisher, classicKS = resultKS, orderBy = "classicKS", ranksOf = "classicFisher", topNodes = length(data_GOdata@graph@nodes))}#if there are fewer than 60 nodes, include all
  
  data_tab_BP <- data_tab_BP[c(1,2,3,6,8)]#extract useful rows and reformat
  colnames(data_tab_BP)[4] <- "Fisher_Rank"
  colnames(data_tab_BP)[5] <- "KS_p"
  
  #########################################################################################
  #If you want to make a node diagram, run the following. Otherwise, skip ahead:
  #########################################################################################
  
  pdf_temp <- paste(Comparison,"_",Test,"_",Condition,"_Convergence_BP.pdf",sep="")#Title for your PDF (file name)
  main_title <- paste("\n\n",Comparison," Convergence Gene Ontology - Biological Process\n",Test," || ",Condition,sep="")#pdf header
  
  setwd(GO_PDF_dir)#PDF output directory
  pdf(pdf_temp) 
  par(cex = 0.4)
  data_node <- showSigOfNodes(data_GOdata, score(resultKS), firstSigNodes = 4, useInfo = "def", useFullNames = TRUE, .NO.CHAR = 60)#Make a pdf of the top 4 most significant GO terms and their parent terms
  par(cex = 0.6)
  title(main=main_title)
  dev.off()#save PDF
  
  #########################################################################################
  #########################################################################################
  
  #Output gene ontology raw data:
  setwd(GO_out_dir)
  output_name <- paste(Comparison,"_",Test,"_",Condition,"_Convergence_BP.txt",sep="")
  write.table(data_tab_BP,output_name, row.names = F, sep = "\t", quote = F)
}

#We run the loop again, this time for molecular functions:
for(file in GO_inputs){
  setwd(annotated_dir)#Working in the directory of your annotations
  data <- read.delim(file, sep = "\t", header = T, fill = T)#input each file
  
  remove <- paste("_",Test,"_",Condition,"_Annotations.txt",sep="")
  Comparison <- str_remove(file,remove)
  
  require(topGO)
  require(org.At.tair.db)#load the TAIR gene ontology data
  
  data <- data[c(1,9)]#from the annotation data, we only need the name of the arabidopsis homolog and it's p-value from the test
  data <- data %>% as.data.frame()
  
  #Create a genelist object:
  data_geneList <- data$p_val#list of gene p-values
  names(data_geneList) <- data$homolog#homolog names
  
  #We will conduct two gene ontologies per comparison. One for molecular function GO terms, and one for molecular function GO terms. We'll start with MF.
  data_Selection <- function(allScore){ return(allScore < 0.05)} # function that returns TRUE/FALSE for p-values<0.05. All the genes are already significant anyways so this step doesn't actually do anything except put the data in a usable format.
  data_All <- annFUN.org(whichOnto="MF", feasibleGenes=NULL, mapping="org.At.tair.db", ID="GeneName")#specify your ontology. We want molecular function mapped to the TAIR database.
  data_GOdata <- new("topGOdata",#make a topGO object using the inputs from above
                     ontology = "MF",
                     allGenes = data_geneList,
                     geneSelectionFun = data_Selection,
                     annot = annFUN.org, mapping = "org.At.tair.db")
  #We conduct two statistical tests for enrichment, fisher's exact test (ignores p-values and just looks at enrichment) and the Kolomogorov-Smirnov test (uses the p-value of genes to rank GO terms and then looks for rank deviations from the full population of TAIR genes).
  resultFisher <- runTest(data_GOdata, algorithm = "classic", statistic = "fisher")
  resultKS <- runTest(data_GOdata, algorithm = "classic", statistic = "ks")
  #Make an ontology table:
  if(length(data_GOdata@graph@nodes) > 59){#if there are more than 60 nodes, include the top 60 only.
    data_tab_MF <- GenTable(data_GOdata, classicFisher = resultFisher, classicKS = resultKS, orderBy = "classicKS", ranksOf = "classicFisher", topNodes = 60)
  } else {data_tab_MF <- GenTable(data_GOdata, classicFisher = resultFisher, classicKS = resultKS, orderBy = "classicKS", ranksOf = "classicFisher", topNodes = length(data_GOdata@graph@nodes))}#if there are fewer than 60 nodes, include all
  
  
  data_tab_MF <- data_tab_MF[c(1,2,3,6,8)]#extract useful rows and reformat
  colnames(data_tab_MF)[4] <- "Fisher_Rank"
  colnames(data_tab_MF)[5] <- "KS_p"
  
  #########################################################################################
  #If you want to make a node diagram, run the following. Otherwise, hash it out.
  #########################################################################################
  
  pdf_temp <- paste(Comparison,"_",Test,"_",Condition,"_Convergence_MF.pdf",sep="")#Title for your PDF (file name)
  main_title <- paste("\n\n",Comparison," Convergence Gene Ontology - Molecular Function\n",Test," || ",Condition,sep="")#pdf header
  setwd(GO_PDF_dir)#PDF output directory
  pdf(pdf_temp) 
  par(cex = 0.4)
  data_node <- showSigOfNodes(data_GOdata, score(resultKS), firstSigNodes = 4, useInfo = "def", useFullNames = TRUE, .NO.CHAR = 60)#Make a pdf of the top 4 most significant GO terms and their parent terms
  par(cex = 0.6)
  title(main=main_title)
  dev.off()#save PDF
  
  #########################################################################################
  #########################################################################################
  
  #Output gene ontology raw data:
  setwd(GO_out_dir)
  output_name <- paste(Comparison,"_",Test,"_",Condition,"_Convergence_MF.txt",sep="")
  write.table(data_tab_MF,output_name, row.names = F, sep = "\t", quote = F)
}

#That's all folks!
 
