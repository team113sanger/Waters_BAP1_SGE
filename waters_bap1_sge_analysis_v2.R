#SGE BAP1 analysis code -- code to take variant counts from the QUANTS pipeline, extensively annotate, and compute SGE functional scores and classifications. 

# START 

#dependencies that may or may not be required
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(DEGreport)
library(pheatmap)
library(GGally)
library(RColorBrewer)
library(cowplot)
library(plyr)
library(ggpubr)
library(stringr)
library(pROC)
library(FSA)

#Data preparation
#makes lists of file locations to unzip, then run the output list of commands using termincal
gen.count.locations <- function(exon,sg) {
  paste0("cd /Users/aw28/Documents/sge_all_bap_0620/BAP1_hiseq_data/counts/","W",exon,"_",sg,"/results/pycroquet
gunzip *.query_to_library_counts.tsv.gz")
}
#Function to apply the above function to listed exons and designated sgRNA_A
a_locs<-lapply (c(1,2,3,5,6,7,8,9,10,14,15,16), function(y) {gen.count.locations(exon=y, sg="A")
})
#Function to apply the above function to listed exons and designated sgRNA_B
b_locs<-lapply (c(1,2,3,5,6,7,8,9,10,14,15,16), function(y) {gen.count.locations(exon=y, sg="B")
})
#make a file of all the files directories to gunzip
lapply(a_locs, write, "./a_locs.txt", append=TRUE)
lapply(b_locs, write, "./b_locs.txt", append=TRUE)

#Function to read in separate count.csv files, $PATH will need to be changed. 
read.counts <- function(exon,sg) {
  #directories 
  dir=paste0("/Users/aw28/Documents/BAP1_analysis/counts_sge/","W",exon,"_",sg,"/results/pycroquet/")
  dir_plasmid=paste0("/Users/aw28/Documents/BAP1_analysis/counts_plasmid/","W",exon,"_",sg,"/results/pycroquet/")
  dir_out=paste0("/Users/aw28/Documents/BAP1_analysis/post_review_updated_analysis/")
  out=paste0(dir_out, "E", exon, "_SG", sg, "_count_frame.csv")
  #makes a list of count file names held in the dir 
  temp = list.files(dir, pattern="*rep_1.query_to_library_counts.tsv|rep_2.query_to_library_counts.tsv|rep_3.query_to_library_counts.tsv|R1.query_to_library_counts.tsv|R2.query_to_library_counts.tsv|R3.query_to_library_counts.tsv")
  temp_plasmid = list.files(dir_plasmid, pattern="*.counts.tsv")
  #makes a list of full paths to each file
  file_all_count=paste0(dir,temp)
  file_all_count_plasmid=paste0(dir_plasmid,temp_plasmid)
  #reads the data from these paths into a list of dataframes
  myfiles = lapply(file_all_count, read.table, skip = 2, header =T, comment.char = "")
  myfiles_plasmid = lapply(file_all_count_plasmid, read.table, skip = 2, header =T, comment.char = "", fill=TRUE)
  myfiles_plasmid = lapply(myfiles_plasmid,"[",c(1,6))
  #makes a large dataframe with all counts - n.b. be careful of the ordering of the new column names make sure agrees to the order in the lists above - this version is ascending numerical
  df <- myfiles %>% purrr::reduce(full_join, by=c("X.id")) %>% dplyr::select(contains("X.id")|contains("rep")|contains("BAP1"))
  df_plasmid <- myfiles_plasmid %>% purrr::reduce(full_join, by=c("X.id")) %>% dplyr::select(contains("X.id")|contains("read"))
  df<- df %>% left_join(df_plasmid, by=c("X.id"))
  df <- df %>% `colnames<-` (c("id", "D10R1","D10R2","D10R3","D14R1","D14R2","D14R3","D21R1","D21R2","D21R3","D4R1","D4R2","D4R3","D7R1","D7R2","D7R3","PLASMID")) %>% replace(is.na(.), 0) 
  #saves dataframe to output location and re-imports it as a named df with names taken from definitions above
  write.csv(df, file=out, row.names = FALSE)
  x=paste0("count_table_", "E", exon, "_SG", sg)
  assign(x,read.csv(out), envir = .GlobalEnv)
}

###############-----------RUN FOR ALL LIBRARIES_STRART-------------###########################################################
#Function to apply the above function to listed exons and designated sgRNA_A
lapply (c(1,2,3,4,5,6,7,9,10,"11_1","11_2","12_1","12_2","13_1","13_2","13_3",14,15,16,"17_1","17_2"), function(y) {read.counts(exon=y, sg="A")
})
#Function to apply the above function to listed exons and designated sgRNA_B
lapply (c(1,2,4,5,6,7,8,9,10,"11_1","11_2","12_1","12_2","13_1","13_2","13_3",14,15,16,"17_1","17_2"), function(y) {read.counts(exon=y, sg="B")
})
###############-----------RUN FOR ALL LIBRARIES_END----------------###########################################################

###############----------------------------------------------------###########################################################
###############-----------COUNTS_END-------------------------------###########################################################
###############----------------------------------------------------###########################################################

###############----------------------------------------------------###########################################################
###############-----------VCF_CLEANING_FOR_VEP_INPUT_START---------###########################################################
###############----------------------------------------------------###########################################################

#Define directories where VaLiAnT VCFs are stored
dir_vcf=paste0("/Users/aw28/Desktop/BAP1_analysis/vep_formatted/")
vcf_out=paste0("/Users/aw28/Desktop/BAP1_analysis/vep_input/")
#makes a list of count file names held in the dir 
temp_vcf = list.files(dir_vcf, pattern="*.vcf")
#makes a list of full paths to each file
file_all_vcf=paste0(dir_vcf,temp_vcf)
#reads the data from these paths into a list of dataframes
my_vcf_files = sapply(file_all_vcf, read.table, comment.char = "#", header = F, simplify=FALSE)
#names columns
colnames_vcf <- c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")
my_vcf_files = lapply(my_vcf_files, setNames, colnames_vcf)
#gets list of directory names so one knows which targetons are which dataframe
vcf_names <-names(my_vcf_files)
#gets file names from full paths
vcf_names_base <-basename(vcf_names)
rm(vcf_names)

#########################################################################
#FUNCTIONS
#########################################################################

#Function to extract the oligo ID from the INFO field and place in the ID column
wrangle.vcf<-function(number) {
  tmp_vcf_data<-my_vcf_files[[number]]
  tmp_vcf_data[c('INFO1', 'INFO2', 'INFO3')] <- str_split_fixed(tmp_vcf_data$INFO, ';', 3)
  tmp_vcf_data$INFO2<-str_remove(tmp_vcf_data$INFO2, "SGE_OLIGO=")
  tmp_vcf_data[ ,c('INFO1', 'INFO3', 'ID')] <- list(NULL)
  tmp_vcf_data$ID = tmp_vcf_data$INFO2
  tmp_vcf_data$INFO2 = NULL
  tmp_vcf_data <- tmp_vcf_data %>% relocate(ID, .before = REF)
  colnames_vcf_out <- c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")
  setNames(tmp_vcf_data, colnames_vcf_out)
}

#run the above function over all VCFs
my_vcf_files = sapply (c(1:44), function(z) {wrangle.vcf(number=z)
}, simplify = FALSE, USE.NAMES=TRUE)
#name the dataframes based on the original VCF file 
names(my_vcf_files) <- vcf_names_base
#function to output each VCF as a separate dataframe
make.tables<-function(object_name) {
  table_output_vcf<-my_vcf_files[[object_name]]
  suffix<-names(my_vcf_files)[object_name]
  write.table(table_output_vcf, file=paste0(vcf_out,"VEP_INPUT_", suffix), sep = "\t", row.names = FALSE, quote = FALSE)
}
#run the above function over all VCFs to output appropriately named VEP INPUTS 
lapply (c(1:44), function(zy) {make.tables(object_name=zy)
})

###############----------------------------------------------------###########################################################
###############-----------VCF_CLEANING_FOR_VEP_INPUT_END-----------###########################################################
###############----------------------------------------------------###########################################################

###############----------------------------------------------------###########################################################
###############-----------COUNT_FILTERING_INPUT_START--------------###########################################################
###############----------------------------------------------------###########################################################

#directoy of metafiles
dir_meta=paste0("/Users/aw28/Documents/BAP1_analysis/meta/")
#makes a list of count file names held in the dir 
temp_meta = list.files(dir_meta, pattern="*.csv")
#makes a list of full paths to each file
file_all_meta=paste0(dir_meta,temp_meta)
#reads the data from these paths into a list of dataframes
my_meta_files = sapply(file_all_meta, read.csv, header =T, comment.char = "#", simplify=FALSE)
#meta names, full paths
meta_names <-names(my_meta_files)
#gets file names from full paths
meta_names_base <-basename(meta_names)
rm(meta_names)
#name the dataframes based on the original VCF file 
names(my_meta_files) <- meta_names_base
#function to output each VCF as a separate dataframe
make.tables<-function(object_name) {
  table_output_meta<-my_meta_files[[object_name]]
  identifier<-names(my_meta_files)[object_name]
  names(table_output_meta)[names(table_output_meta) == 'oligo_name'] <- 'id'
  assign(identifier, table_output_meta, envir = .GlobalEnv)
}
#run the above function over all metas to output appropriately named metas into the global environment
lapply (c(1:44), function(zy) {make.tables(object_name=zy)
})

#merge meta and counts_for_SGA
E1_SGA_master<-merge(count_table_E1_SGA, chr3_52409736_52409980_minus_sgRNA_w1_a_meta.csv, by="id", all.x=FALSE)
E2_SGA_master<-merge(count_table_E2_SGA, chr3_52409602_52409846_minus_sgRNA_w2_a_meta.csv, by="id", all.x=FALSE)
E3_SGA_master<-merge(count_table_E3_SGA, chr3_52409466_52409711_minus_sgRNA_w3_a_meta.csv, by="id", all.x=FALSE)
E4_SGA_master<-merge(count_table_E4_SGA, chr3_52408417_52408661_minus_sgRNA_w4_a_meta.csv, by="id", all.x=FALSE)
E5_SGA_master<-merge(count_table_E5_SGA, chr3_52407887_52408131_minus_sgRNA_w5_a_meta.csv, by="id", all.x=FALSE)
E6_SGA_master<-merge(count_table_E6_SGA, chr3_52407284_52407528_minus_sgRNA_w6_a_meta.csv, by="id", all.x=FALSE)
E7_SGA_master<-merge(count_table_E7_SGA, chr3_52407123_52407367_minus_sgRNA_w7_a_meta.csv, by="id", all.x=FALSE)
E9_SGA_master<-merge(count_table_E9_SGA, chr3_52406199_52406443_minus_sgRNA_w9_a_meta.csv, by="id", all.x=FALSE)
E10_SGA_master<-merge(count_table_E10_SGA, chr3_52405725_52405969_minus_sgRNA_w10_a_meta.csv, by="id", all.x=FALSE)
E11_1_SGA_master<-merge(count_table_E11_1_SGA, chr3_52405129_52405373_minus_sgRNA_w11_x_a_meta.csv, by="id", all.x=FALSE)
E11_2_SGA_master<-merge(count_table_E11_2_SGA, chr3_52405033_52405277_minus_sgRNA_w11_x_a_meta.csv, by="id", all.x=FALSE)
E12_1_SGA_master<-merge(count_table_E12_1_SGA, chr3_52404452_52404696_minus_sgRNA_w12_x_a_meta.csv, by="id", all.x=FALSE)
E12_2_SGA_master<-merge(count_table_E12_2_SGA, chr3_52404348_52404592_minus_sgRNA_w12_x_a_meta.csv, by="id", all.x=FALSE)
E13_1_SGA_master<-merge(count_table_E13_1_SGA, chr3_52403738_52403982_minus_sgRNA_w13_1_a_meta.csv, by="id", all.x=FALSE)
E13_2_SGA_master<-merge(count_table_E13_2_SGA, chr3_52403556_52403800_minus_sgRNA_w13_2_a_meta.csv, by="id", all.x=FALSE)
E13_3_SGA_master<-merge(count_table_E13_3_SGA, chr3_52403378_52403622_minus_sgRNA_w13_3_a_meta.csv, by="id", all.x=FALSE)
E14_SGA_master<-merge(count_table_E14_SGA, chr3_52403097_52403341_minus_sgRNA_w14_a_meta.csv, by="id", all.x=FALSE)
E15_SGA_master<-merge(count_table_E15_SGA, chr3_52402703_52402947_minus_sgRNA_w15_a_meta.csv, by="id", all.x=FALSE)
E16_SGA_master<-merge(count_table_E16_SGA, chr3_52402511_52402755_minus_sgRNA_w16_a_meta.csv, by="id", all.x=FALSE)
E17_1_SGA_master<-merge(count_table_E17_1_SGA, chr3_52402282_52402526_minus_sgRNA_w17_x_a_meta.csv, by="id", all.x=FALSE)
E17_2_SGA_master<-merge(count_table_E17_2_SGA, chr3_52402173_52402418_minus_sgRNA_w17_x_a_meta.csv, by="id", all.x=FALSE)


#merge meta and counts_for_SGB
E1_SGB_master<-merge(count_table_E1_SGB, chr3_52409736_52409980_minus_sgRNA_w1_b_meta.csv, by="id", all.x=FALSE)
E2_SGB_master<-merge(count_table_E2_SGB, chr3_52409602_52409846_minus_sgRNA_w2_b_meta.csv, by="id", all.x=FALSE)
E4_SGB_master<-merge(count_table_E4_SGB, chr3_52408417_52408661_minus_sgRNA_w4_b_meta.csv, by="id", all.x=FALSE)
E5_SGB_master<-merge(count_table_E5_SGB, chr3_52407887_52408131_minus_sgRNA_w5_b_meta.csv, by="id", all.x=FALSE)
E6_SGB_master<-merge(count_table_E6_SGB, chr3_52407284_52407528_minus_sgRNA_w6_b_meta.csv, by="id", all.x=FALSE)
E7_SGB_master<-merge(count_table_E7_SGB, chr3_52407123_52407367_minus_sgRNA_w7_b_meta.csv, by="id", all.x=FALSE)
E8_SGB_master<-merge(count_table_E8_SGB, chr3_52406736_52406980_minus_sgRNA_w8_b_meta.csv, by="id", all.x=FALSE)
E9_SGB_master<-merge(count_table_E9_SGB, chr3_52406199_52406443_minus_sgRNA_w9_b_meta.csv, by="id", all.x=FALSE)
E10_SGB_master<-merge(count_table_E10_SGB, chr3_52405725_52405969_minus_sgRNA_w10_b_meta.csv, by="id", all.x=FALSE)
E11_1_SGB_master<-merge(count_table_E11_1_SGB, chr3_52405129_52405373_minus_sgRNA_w11_x_b_meta.csv, by="id", all.x=FALSE)
E11_2_SGB_master<-merge(count_table_E11_2_SGB, chr3_52405033_52405277_minus_sgRNA_w11_x_b_meta.csv, by="id", all.x=FALSE)
E12_1_SGB_master<-merge(count_table_E12_1_SGB, chr3_52404452_52404696_minus_sgRNA_w12_x_b_meta.csv, by="id", all.x=FALSE)
E12_2_SGB_master<-merge(count_table_E12_2_SGB, chr3_52404348_52404592_minus_sgRNA_w12_x_b_meta.csv, by="id", all.x=FALSE)
E13_1_SGB_master<-merge(count_table_E13_1_SGB, chr3_52403738_52403982_minus_sgRNA_w13_1_b_meta.csv, by="id", all.x=FALSE)
E13_2_SGB_master<-merge(count_table_E13_2_SGB, chr3_52403556_52403800_minus_sgRNA_w13_2_b_meta.csv, by="id", all.x=FALSE)
E13_3_SGB_master<-merge(count_table_E13_3_SGB, chr3_52403378_52403622_minus_sgRNA_w13_3_b_meta.csv, by="id", all.x=FALSE)
E14_SGB_master<-merge(count_table_E14_SGB, chr3_52403097_52403341_minus_sgRNA_w14_b_meta.csv, by="id", all.x=FALSE)
E15_SGB_master<-merge(count_table_E15_SGB, chr3_52402703_52402947_minus_sgRNA_w15_b_meta.csv, by="id", all.x=FALSE)
E16_SGB_master<-merge(count_table_E16_SGB, chr3_52402511_52402755_minus_sgRNA_w16_b_meta.csv, by="id", all.x=FALSE)
E17_1_SGB_master<-merge(count_table_E17_1_SGB, chr3_52402282_52402526_minus_sgRNA_w17_x_b_meta.csv, by="id", all.x=FALSE)
E17_2_SGB_master<-merge(count_table_E17_2_SGB, chr3_52402173_52402418_minus_sgRNA_w17_x_b_meta.csv, by="id", all.x=FALSE)

#remove rows that have fewer than 10 counts total across all 15 samples - OPTIONAL REMOVE VARIANTS IN CONSTANT REGTIONS
filter.counts<-function(exon, sg){
  master_input<-paste0("E",exon,"_SG",sg,"_master")
  filtered_counts <- get(master_input)
  filtered_counts$rowsum <- rowSums(filtered_counts[, c("D4R1", "D4R2", "D4R3", "D7R1", "D7R2", "D7R3", "D10R1", "D10R2", "D10R3", "D14R1",  "D14R2",  "D14R3",  "D21R1",  "D21R2",  "D21R3")])
  filtered_counts<- filtered_counts %>% filter(rowsum>=10)
  #optional_remove high counts
  #filtered_counts<- filtered_counts %>% filter(rowsum<=10000)
  #OPTIONAL - REMOVE VARIANTS PAM MUTATION CODONS
  #filtered_counts<- filtered_counts %>% filter(pam_mut_sgrna_id %in% "" | is.na(pam_mut_sgrna_id))
  #OPTIONAL - REMOVE VARIANTS IN CONSTANT REGIONS - CAN REMOVE THIS STEP IF ADAPTERS STRIPPED SAME AS SEQUENCING PRIMERS
  filtered_counts<- filtered_counts %>% filter(vcf_var_in_const==0 | is.na(vcf_var_in_const))
  tablename<-paste0(master_input,"_filtered_counts")
  assign(tablename, filtered_counts, envir = .GlobalEnv)
}

#Function to apply the above function to listed exons and designated sgRNA_A
lapply (c(1,2,3,4,5,6,7,9,10,"11_1","11_2","12_1","12_2","13_1","13_2","13_3",14,15,16,"17_1","17_2"), function(y) {filter.counts(exon=y, sg="A")
})
#Function to apply the above function to listed exons and designated sgRNA_B
lapply (c(1,2,4,5,6,7,8,9,10,"11_1","11_2","12_1","12_2","13_1","13_2","13_3",14,15,16,"17_1","17_2"), function(y) {filter.counts(exon=y, sg="B")
})



###############----------------------------------------------------###########################################################
###############-----------VEP_START--------------------------------###########################################################
###############----------------------------------------------------###########################################################

#merging vep outputs to create vep + meta
#Define directories where VEP OUTPUTS are stored
dir_vep=paste0("/Users/aw28/Documents/BAP1_analysis/vep_output/")
#makes a list of count file names held in the dir 
temp_vep = list.files(dir_vep, pattern="*.tsv")
#makes a list of full paths to each file
file_all_vep=paste0(dir_vep,temp_vep)
#reads the data from these paths into a list of dataframes
my_vep_files = sapply(file_all_vep, read.table, skip = 92, header =T, comment.char = "", simplify=FALSE)
#meta names, full paths
vep_names <-names(my_vep_files)
#gets file names from full paths
vep_names_base <-basename(vep_names)
rm(vep_names)
#name the dataframes based on the original VEP file 
names(my_vep_files) <- vep_names_base
#function to output each VCF as a separate dataframe
make.vep.tables<-function(object_name) {
  table_output_vep<-my_vep_files[[object_name]]
  identifier<-names(my_vep_files)[object_name]
  names(table_output_vep)[names(table_output_vep) == 'X.Uploaded_variation'] <- 'id'
  table_output_vep <- table_output_vep %>% filter(Feature %in% "ENST00000460680.6")
  assign(identifier, table_output_vep, envir = .GlobalEnv)
}
#run the above function over all metas to output appropriately named VEPs into the global environment
lapply (c(1:44), function(zy) {make.vep.tables(object_name=zy)
})

#merge filtered master tables with VEP annotation tables and syn_filter
#merge filtered counts + meta with vep output to give filtered annotated dataframes SGA 
E1_SGA_master_annotated<-merge(E1_SGA_master_filtered_counts, VEP_OUTPUT_chr3_52409736_52409980_minus_sgRNA_w1_a_pam.tsv, by="id", all.x=TRUE)
E2_SGA_master_annotated<-merge(E2_SGA_master_filtered_counts, VEP_OUTPUT_chr3_52409602_52409846_minus_sgRNA_w2_a_pam.tsv, by="id", all.x=TRUE)
E3_SGA_master_annotated<-merge(E3_SGA_master_filtered_counts, VEP_OUTPUT_chr3_52409466_52409711_minus_sgRNA_w3_a_pam.tsv, by="id", all.x=TRUE)
E4_SGA_master_annotated<-merge(E4_SGA_master_filtered_counts, VEP_OUTPUT_chr3_52408417_52408661_minus_sgRNA_w4_a_pam.tsv, by="id", all.x=TRUE)
E5_SGA_master_annotated<-merge(E5_SGA_master_filtered_counts, VEP_OUTPUT_chr3_52407887_52408131_minus_sgRNA_w5_a_pam.tsv, by="id", all.x=TRUE)
E6_SGA_master_annotated<-merge(E6_SGA_master_filtered_counts, VEP_OUTPUT_chr3_52407284_52407528_minus_sgRNA_w6_a_pam.tsv, by="id", all.x=TRUE)
E7_SGA_master_annotated<-merge(E7_SGA_master_filtered_counts, VEP_OUTPUT_chr3_52407123_52407367_minus_sgRNA_w7_a_pam.tsv, by="id", all.x=TRUE)
E9_SGA_master_annotated<-merge(E9_SGA_master_filtered_counts, VEP_OUTPUT_chr3_52406199_52406443_minus_sgRNA_w9_a_pam.tsv, by="id", all.x=TRUE)
E10_SGA_master_annotated<-merge(E10_SGA_master_filtered_counts, VEP_OUTPUT_chr3_52405725_52405969_minus_sgRNA_w10_a_pam.tsv, by="id", all.x=TRUE)
E11_1_SGA_master_annotated<-merge(E11_1_SGA_master_filtered_counts, VEP_OUTPUT_chr3_52405129_52405373_minus_sgRNA_w11_x_a_pam.tsv, by="id", all.x=TRUE)
E11_2_SGA_master_annotated<-merge(E11_2_SGA_master_filtered_counts, VEP_OUTPUT_chr3_52405033_52405277_minus_sgRNA_w11_x_a_pam.tsv, by="id", all.x=TRUE)
E12_1_SGA_master_annotated<-merge(E12_1_SGA_master_filtered_counts, VEP_OUTPUT_chr3_52404452_52404696_minus_sgRNA_w12_x_a_pam.tsv, by="id", all.x=TRUE)
E12_2_SGA_master_annotated<-merge(E12_2_SGA_master_filtered_counts, VEP_OUTPUT_chr3_52404348_52404592_minus_sgRNA_w12_x_a_pam.tsv, by="id", all.x=TRUE)
E13_1_SGA_master_annotated<-merge(E13_1_SGA_master_filtered_counts, VEP_OUTPUT_chr3_52403738_52403982_minus_sgRNA_w13_1_a_pam.tsv, by="id", all.x=TRUE)
E13_2_SGA_master_annotated<-merge(E13_2_SGA_master_filtered_counts, VEP_OUTPUT_chr3_52403556_52403800_minus_sgRNA_w13_2_a_pam.tsv, by="id", all.x=TRUE)
E13_3_SGA_master_annotated<-merge(E13_3_SGA_master_filtered_counts, VEP_OUTPUT_chr3_52403378_52403622_minus_sgRNA_w13_3_a_pam.tsv, by="id", all.x=TRUE)
E14_SGA_master_annotated<-merge(E14_SGA_master_filtered_counts, VEP_OUTPUT_chr3_52403097_52403341_minus_sgRNA_w14_a_pam.tsv, by="id", all.x=TRUE)
E15_SGA_master_annotated<-merge(E15_SGA_master_filtered_counts, VEP_OUTPUT_chr3_52402703_52402947_minus_sgRNA_w15_a_pam.tsv, by="id", all.x=TRUE)
E16_SGA_master_annotated<-merge(E16_SGA_master_filtered_counts, VEP_OUTPUT_chr3_52402511_52402755_minus_sgRNA_w16_a_pam.tsv, by="id", all.x=TRUE)
E17_1_SGA_master_annotated<-merge(E17_1_SGA_master_filtered_counts, VEP_OUTPUT_chr3_52402282_52402526_minus_sgRNA_w17_x_a_pam.tsv, by="id", all.x=TRUE)
E17_2_SGA_master_annotated<-merge(E17_2_SGA_master_filtered_counts, VEP_OUTPUT_chr3_52402173_52402418_minus_sgRNA_w17_x_a_pam.tsv, by="id", all.x=TRUE)


#merge filtered counts + meta with vep output to give filtered annotated dataframes SGB
E1_SGB_master_annotated<-merge(E1_SGB_master_filtered_counts, VEP_OUTPUT_chr3_52409736_52409980_minus_sgRNA_w1_b_pam.tsv, by="id", all.x=TRUE)
E2_SGB_master_annotated<-merge(E2_SGB_master_filtered_counts, VEP_OUTPUT_chr3_52409602_52409846_minus_sgRNA_w2_b_pam.tsv, by="id", all.x=TRUE)
E4_SGB_master_annotated<-merge(E4_SGB_master_filtered_counts, VEP_OUTPUT_chr3_52408417_52408661_minus_sgRNA_w4_b_pam.tsv, by="id", all.x=TRUE)
E5_SGB_master_annotated<-merge(E5_SGB_master_filtered_counts, VEP_OUTPUT_chr3_52407887_52408131_minus_sgRNA_w5_b_pam.tsv, by="id", all.x=TRUE)
E6_SGB_master_annotated<-merge(E6_SGB_master_filtered_counts, VEP_OUTPUT_chr3_52407284_52407528_minus_sgRNA_w6_b_pam.tsv, by="id", all.x=TRUE)
E7_SGB_master_annotated<-merge(E7_SGB_master_filtered_counts, VEP_OUTPUT_chr3_52407123_52407367_minus_sgRNA_w7_b_pam.tsv, by="id", all.x=TRUE)
E8_SGB_master_annotated<-merge(E8_SGB_master_filtered_counts, VEP_OUTPUT_chr3_52406736_52406980_minus_sgRNA_w8_b_pam.tsv, by="id", all.x=TRUE)
E9_SGB_master_annotated<-merge(E9_SGB_master_filtered_counts, VEP_OUTPUT_chr3_52406199_52406443_minus_sgRNA_w9_b_pam.tsv, by="id", all.x=TRUE)
E10_SGB_master_annotated<-merge(E10_SGB_master_filtered_counts, VEP_OUTPUT_chr3_52405725_52405969_minus_sgRNA_w10_b_pam.tsv, by="id", all.x=TRUE)
E11_1_SGB_master_annotated<-merge(E11_1_SGB_master_filtered_counts, VEP_OUTPUT_chr3_52405129_52405373_minus_sgRNA_w11_x_b_pam.tsv, by="id", all.x=TRUE)
E11_2_SGB_master_annotated<-merge(E11_2_SGB_master_filtered_counts, VEP_OUTPUT_chr3_52405033_52405277_minus_sgRNA_w11_x_b_pam.tsv, by="id", all.x=TRUE)
E12_1_SGB_master_annotated<-merge(E12_1_SGB_master_filtered_counts, VEP_OUTPUT_chr3_52404452_52404696_minus_sgRNA_w12_x_b_pam.tsv, by="id", all.x=TRUE)
E12_2_SGB_master_annotated<-merge(E12_2_SGB_master_filtered_counts, VEP_OUTPUT_chr3_52404348_52404592_minus_sgRNA_w12_x_b_pam.tsv, by="id", all.x=TRUE)
E13_1_SGB_master_annotated<-merge(E13_1_SGB_master_filtered_counts, VEP_OUTPUT_chr3_52403738_52403982_minus_sgRNA_w13_1_b_pam.tsv, by="id", all.x=TRUE)
E13_2_SGB_master_annotated<-merge(E13_2_SGB_master_filtered_counts, VEP_OUTPUT_chr3_52403556_52403800_minus_sgRNA_w13_2_b_pam.tsv, by="id", all.x=TRUE)
E13_3_SGB_master_annotated<-merge(E13_3_SGB_master_filtered_counts, VEP_OUTPUT_chr3_52403378_52403622_minus_sgRNA_w13_3_b_pam.tsv, by="id", all.x=TRUE)
E14_SGB_master_annotated<-merge(E14_SGB_master_filtered_counts, VEP_OUTPUT_chr3_52403097_52403341_minus_sgRNA_w14_b_pam.tsv, by="id", all.x=TRUE)
E15_SGB_master_annotated<-merge(E15_SGB_master_filtered_counts, VEP_OUTPUT_chr3_52402703_52402947_minus_sgRNA_w15_b_pam.tsv, by="id", all.x=TRUE)
E16_SGB_master_annotated<-merge(E16_SGB_master_filtered_counts, VEP_OUTPUT_chr3_52402511_52402755_minus_sgRNA_w16_b_pam.tsv, by="id", all.x=TRUE)
E17_1_SGB_master_annotated<-merge(E17_1_SGB_master_filtered_counts, VEP_OUTPUT_chr3_52402282_52402526_minus_sgRNA_w17_x_b_pam.tsv, by="id", all.x=TRUE)
E17_2_SGB_master_annotated<-merge(E17_2_SGB_master_filtered_counts, VEP_OUTPUT_chr3_52402173_52402418_minus_sgRNA_w17_x_b_pam.tsv, by="id", all.x=TRUE)


#filter the above tables to get normalization tables
normalization.counts<-function(exon, sg){
  master_input<-paste0("E",exon,"_SG",sg,"_master_annotated")
  filtered_counts <- get(master_input)
  filtered_counts <-filtered_counts[!duplicated(filtered_counts$mseq), ]
  filtered_counts<- filtered_counts %>% filter(Consequence %in% "synonymous_variant" | Consequence %in% "intron_variant")
  tablename<-paste0("E",exon,"_SG",sg,"_normalization_counts")
  assign(tablename, filtered_counts, envir = .GlobalEnv)
}

#Function to apply the above function to listed exons and designated sgRNA_A
lapply (c(1,2,3,4,5,6,7,9,10,"11_1","11_2","12_1","12_2","13_1","13_2","13_3",14,15,16,"17_1","17_2"), function(y) {normalization.counts(exon=y, sg="A")
})
#Function to apply the above function to listed exons and designated sgRNA_B
lapply (c(1,2,4,5,6,7,8,9,10,"11_1","11_2","12_1","12_2","13_1","13_2","13_3",14,15,16,"17_1","17_2"), function(y) {normalization.counts(exon=y, sg="B")
})


#make normalization matricies
deseq.norm.input <- function(exon, sg){
  x_normalization_counts <- paste0("E",exon,"_SG",sg,"_normalization_counts")
  x_normalization_counts <-get(x_normalization_counts)
  deseq_input_noramlization_counts <- x_normalization_counts %>% dplyr::select("mseq", "D4R1", "D4R2", "D4R3", "D7R1", "D7R2", "D7R3", "D10R1", "D10R2", "D10R3", "D14R1",  "D14R2",  "D14R3",  "D21R1",  "D21R2",  "D21R3")
  deseq_input_noramlization_counts <- deseq_input_noramlization_counts %>% remove_rownames %>% column_to_rownames(var="mseq")
  deseq_input_noramlization_counts <- as.matrix.data.frame(deseq_input_noramlization_counts)
  normdfname<-paste0("E",exon,"_SG",sg,"_normalization_counts_MATRIX")
  assign(normdfname, deseq_input_noramlization_counts, envir = .GlobalEnv)
}


#Function to apply the above function to listed exons and designated sgRNA_A
lapply (c(1,2,3,4,5,6,7,9,10,"11_1","11_2","12_1","12_2","13_1","13_2","13_3",14,15,16,"17_1","17_2"), function(y) {deseq.norm.input(exon=y, sg="A")
})
#Function to apply the above function to listed exons and designated sgRNA_B
lapply (c(1,2,4,5,6,7,8,9,10,"11_1","11_2","12_1","12_2","13_1","13_2","13_3",14,15,16,"17_1","17_2"), function(y) {deseq.norm.input(exon=y, sg="B")
})



###############----------------------------------------------------###########################################################
###############-----------DESEQ2_START-----------------------------###########################################################
###############----------------------------------------------------###########################################################

#make count matricies
deseq.count.input <- function(exon,sg){
  x_filtered_counts <- paste0("E",exon,"_SG",sg,"_master_annotated")
  x_filtered_counts <- get(x_filtered_counts)
  #line below removes duplicates within the targeton file based on mseq - some mutators produce the same mseq but oligo library contains one instance
  x_filtered_counts <-x_filtered_counts[!duplicated(x_filtered_counts$mseq), ]
  deseq_input_filtered_counts <- x_filtered_counts %>% dplyr::select("mseq", "D4R1", "D4R2", "D4R3", "D7R1", "D7R2", "D7R3", "D10R1", "D10R2", "D10R3", "D14R1",  "D14R2",  "D14R3",  "D21R1",  "D21R2",  "D21R3")
  deseq_input_filtered_counts <- deseq_input_filtered_counts %>% remove_rownames %>% column_to_rownames(var="mseq")
  deseq_input_filtered_counts <- as.matrix.data.frame(deseq_input_filtered_counts)
  countdfname<-paste0("E",exon,"_SG",sg,"_master_annotated_filtered_counts_MATRIX")
  assign(countdfname, deseq_input_filtered_counts, envir = .GlobalEnv)
}


#Function to apply the above function to listed exons and designated sgRNA_A
lapply (c(1,2,3,4,5,6,7,9,10,"11_1","11_2","12_1","12_2","13_1","13_2","13_3",14,15,16,"17_1","17_2"), function(y) {deseq.count.input(exon=y, sg="A")
})
#Function to apply the above function to listed exons and designated sgRNA_B
lapply (c(1,2,4,5,6,7,8,9,10,"11_1","11_2","12_1","12_2","13_1","13_2","13_3",14,15,16,"17_1","17_2"), function(y) {deseq.count.input(exon=y, sg="B")
})


#read in experimental designs SGA
annot_a<-read.csv("/Users/aw28/Documents/BAP1_analysis/experimental_design/annot_a.csv")
annot_a<-as.matrix.data.frame(annot_a)
annot_a_continuous<-read.csv("/Users/aw28/Documents/BAP1_analysis/experimental_design/annot_a_continuous.csv")
annot_a_continuous$condition<-as.numeric(annot_a_continuous$condition)



#########################################################################
#FUNCTIONS_DESEQ_START
#########################################################################

#########################################################################
#OUTPUT FUNCTIONS START
#########################################################################

#PCA SCREE PLOT
plotPCA.hk <- function (object, intgroup = "condition", ntop = 500, pc_1 = 1, pc_2 = 2, returnData = FALSE, scree = FALSE)
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, pc_1], PC2 = pca$x[, pc_2], group = group,
                  intgroup.df, name = colData(object)[,1])
  if (returnData) {
    attr(d, "percentVar") <- percentVar[pc_1:pc_2]
    return(d)
  }
  if (scree) {
    xx<- barplot(round(percentVar, digits = 2)*100, names.arg=c(1:length(percentVar)),xlab="PC",ylab="% Variance",ylim=c(0,100), main="Scree Plot")
    text(x = xx, y = round(percentVar, digits = 4)*100, label = round(percentVar, digits = 4)*100, pos = 3, cex = 0.8, col = "black")
  }
  else {
    ggplot(data = d, aes_string(x = "PC1", y = "PC2", color="condition", shape="type", label = "sample")) + scale_shape_manual(values=seq(0,127)) + geom_point(size = 3) + xlab(paste0("PC",pc_1,": ", round(percentVar[pc_1] * 100, digits = 2), "% variance")) + ylab(paste0("PC",pc_2,": ", round(percentVar[pc_2] * 100, digits=2), "% variance")) + coord_fixed() + geom_text_repel(size=3)
  }
}

#needed for scatter pairs plots
ggpairs_ext <- function(data, mapping, pts=list(), smt=list(), ...){
  ggplot(data = data, mapping = mapping, ...) +
    do.call(geom_point, pts) +
    do.call(geom_smooth, smt)
}

# Function to estimate size factors (typically run for control oligos only)
estimate_control_size_factors <- function ( countData = countData, colData = colData, design = design, minRowSum = 10, ref = NULL ) {
  print( "Estimating control size factors...")
  dds <- DESeqDataSetFromMatrix( countData = countData, 
                                 colData = colData, 
                                 design = as.formula( design ) )
  
  dds <- dds[ rowSums( counts( dds ) ) > minRowSum, ]
  
  if ( ! is.null( ref ) ) {
    dds$condition <- relevel( dds$condition, ref = ref )
  }
  
  control_size_factors <- sizeFactors( estimateSizeFactors( dds ) )
  
  return( control_size_factors )
}

#function to add text to column names
appendDataFrameColumns<-function(df, prefix='', suffix='', sep='')
{
  colnames(df) <- paste(prefix, colnames(df), suffix, sep=sep)
  
  return(df)
}

#########################################################################
#OUTPUT FUNCTIONS END
#########################################################################

#the main deseq funtion - setwd() to where you want all the files to be outputted
run.deseq<-function(exon_assay, sg_assay){
  x_normalization_counts = paste0("E",exon_assay,"_","SG",sg_assay,"_normalization_counts_MATRIX")
  x_filtered_counts = paste0("E",exon_assay,"_","SG",sg_assay,"_master_annotated_filtered_counts_MATRIX")
  x_normalization_counts=get(x_normalization_counts)
  x_filtered_counts =get(x_filtered_counts )
  # Get control size factors
  SF <- estimate_control_size_factors(  countData = x_normalization_counts,
                                        colData = annot_a,
                                        ref = "D4",
                                        design = "~ condition")
  #Prep DESEQ
  dds <- DESeqDataSetFromMatrix(countData = x_filtered_counts, colData = annot_a, design = ~condition)
  dds$condition <- factor(dds$condition, levels=c("D4", "D7", "D10","D14","D21"))
  dds$condition <- relevel(dds$condition, ref = "D4")
  dds <- estimateSizeFactors(dds)
  sizeFactors(dds) <- SF
  #Run Deseq
  dds <- DESeq(dds)
  rld <- rlog(dds)
  res<-results(dds)
  res<-as.data.frame(res)
  z_score <- assay(rld) %>% as.matrix() %>% t() %>% scale() %>% t() %>% as.data.frame()
  colnames(z_score) <- paste0(colnames(z_score), "_z_score")
  z_score$mseq <- row.names(z_score)  
  row.names(z_score) <- NULL
  z_score <-  z_score %>% select(mseq, everything())
  #shirnkage type set to normal, apeglm would be better for RNASeq - but have selected no shrinkage is table
  table_wald <- degComps(dds, combs = "condition", contrast = list("condition_D7_vs_D4","condition_D10_vs_D4" , "condition_D14_vs_D4", "condition_D21_vs_D4"), alpha = 0.05, skip = FALSE, type = "normal", pairs = FALSE, fdr = "default")
  
  #Summary Table
  #production of summary table - RAW - change to shrunken if want LFC shrinkage - apeglm is most appropriate for RNASeq
  summary <- purrr::reduce(c(deg(table_wald[[1]], "raw") %>% appendDataFrameColumns(suffix="_D4_D7") %>% as.data.frame() %>% rownames_to_column(var="mseq") %>% list(),
                             deg(table_wald[[2]], "raw") %>% appendDataFrameColumns(suffix="_D4_D10") %>% as.data.frame() %>% rownames_to_column(var="mseq") %>% list(),
                             deg(table_wald[[3]], "raw") %>% appendDataFrameColumns(suffix="_D4_D14") %>% as.data.frame() %>% rownames_to_column(var="mseq") %>% list(),
                             deg(table_wald[[4]], "raw") %>% appendDataFrameColumns(suffix="_D4_D21") %>% as.data.frame() %>% rownames_to_column(var="mseq") %>% list()), left_join, by="mseq") %>% select(-c("baseMean_D4_D10","baseMean_D4_D14","baseMean_D4_D21")) %>% dplyr::rename(baseMean = "baseMean_D4_D7") %>% data.frame()
  #bind the z_scores to the summary table
  summary<-merge(summary, z_score, by="mseq", all.x=TRUE)
  
  #OUTPUT PLOTS AND TABLES 
  ident<-paste0("E",exon_assay,"_","SG",sg_assay)
  #SCREE
  pdf(paste0(ident,"_scree.pdf"))
  plotPCA.hk(rld,intgroup=c("condition", "type"), returnData=FALSE,pc_1=1, pc_2=2, scree=TRUE)
  dev.off()
  #DISPERSION
  pdf(paste0(ident,"_dispersion.pdf"))
  plotDispEsts(dds, ylim =c(1e-4,2e1))
  dev.off()
  #HEATMAP
  sampleDistMatrix <- as.matrix( dist( t( assay(rld) ) ) )
  pdf(paste0(ident,"_heatmap.pdf"))
  pheatmap(sampleDistMatrix, trace="none", col=colorRampPalette(rev(brewer.pal(9, "Blues")) )(255), adjRow = c(1,1))
  dev.off()
  remove(sampleDistMatrix)
  ###Scatter plot, checking replicate consistency
  pdf(paste0(ident,"_scatter_matrix.pdf"), width=8, height=8)
  print(ggpairs(as.data.frame(assay(rld)) %>% select(starts_with("D4")), lower = list(continuous = wrap(ggpairs_ext, pts=list(size=0.4,colour="black"), smt=list(method="lm", se=F, size=0.2, colour="blue"))), title = "Regularized-log Transformed Read Count") + theme(plot.title = element_text(hjust = 0.5, size=13, face="bold")))
  print(ggpairs(as.data.frame(assay(rld)) %>% select(starts_with("D7")), lower = list(continuous = wrap(ggpairs_ext, pts=list(size=0.4,colour="black"), smt=list(method="lm", se=F, size=0.2, colour="blue"))), title = "Regularized-log Transformed Read Count") + theme(plot.title = element_text(hjust = 0.5, size=13, face="bold")))
  print(ggpairs(as.data.frame(assay(rld)) %>% select(starts_with("D10")), lower = list(continuous = wrap(ggpairs_ext, pts=list(size=0.4,colour="black"), smt=list(method="lm", se=F, size=0.2, colour="blue"))), title = "Regularized-log Transformed Read Count" ) + theme(plot.title = element_text(hjust = 0.5, size=13, face="bold")))
  print(ggpairs(as.data.frame(assay(rld)) %>% select(starts_with("D14")), lower = list(continuous = wrap(ggpairs_ext, pts=list(size=0.4,colour="black"), smt=list(method="lm", se=F, size=0.2, colour="blue"))), title = "Regularized-log Transformed Read Count") + theme(plot.title = element_text(hjust = 0.5, size=13, face="bold")))
  print(ggpairs(as.data.frame(assay(rld)) %>% select(starts_with("D21")), lower = list(continuous = wrap(ggpairs_ext, pts=list(size=0.4,colour="black"), smt=list(method="lm", se=F, size=0.2, colour="blue"))), title = "Regularized-log Transformed Read Count") + theme(plot.title = element_text(hjust = 0.5, size=13, face="bold")))
  dev.off()
  
  #PCA PLOTS
  pcaData <- plotPCA(rld, intgroup=c("condition", "type"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  pca<-ggplot(pcaData, aes(PC1, PC2, color=condition, shape=type)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed()+
    theme_classic()
  ggsave(paste0(ident,"_PCA.pdf"), pca, height=6, width=8)
  
  #Export rld
  write.table(as.data.frame(assay(rld)), sep="\t",file=paste0(ident,"_rld.txt"), col.names=NA)
  #Export normalized read count
  write.table(counts(dds,normalized=TRUE), sep="\t",file=paste0(ident,"_norm_count.txt"), col.names=NA)
  #Export size factor
  write.table(dds@assays@data@listData %>% as.data.frame(),sep="\t",file=paste0(ident,"_normalization_table.txt"), col.names=NA)
  #Export full table
  write.table(dds@rowRanges@elementMetadata@listData %>% as.data.frame() ,sep="\t",file=paste0(ident,"_disper_table.txt"), col.names=NA)
  
  #continuous DESEQ to add SGE RATE 
  # Get control size factors
  SF <- estimate_control_size_factors(  countData = x_normalization_counts,
                                        colData = annot_a_continuous,
                                        design = "~ condition")
  #Prep DESEQ
  dds <- DESeqDataSetFromMatrix(countData = x_filtered_counts, colData = annot_a_continuous, design = ~condition)
  dds <- estimateSizeFactors(dds)
  sizeFactors(dds) <- SF
  #Run Deseq - SHRUNKEN LFC - CHANGE TO table_wald[[2]] if want shrunken LFC, type="apgelm" is appropriate for RNASeq perhaps not SGE
  dds <- DESeq(dds)
  table_wald <- degComps(dds, combs = "condition", alpha = 0.05, skip = FALSE, type = "normal", pairs = FALSE, fdr = "default")
  rate <- as.data.frame(table_wald[[1]])%>%rownames_to_column(var="mseq")%>%appendDataFrameColumns(suffix="_continuous")
  colnames(`rate`)[colnames(`rate`) == "mseq_continuous"] <- "mseq"
  summary<-merge(summary, rate, by="mseq", all.x=TRUE)
  
  
  #Substract the LFC median of the synonymous_variant and intron_variant
  adj<- as.data.frame(x_normalization_counts)
  adj<-adj %>% mutate(mseq=rownames(adj)) %>% select(mseq, everything()) %>% remove_rownames() %>% select(1)
  adj <- adj %>% left_join(summary %>%select("mseq","log2FoldChange_D4_D7","log2FoldChange_D4_D10","log2FoldChange_D4_D14","log2FoldChange_D4_D21", "log2FoldChange_continuous")) %>% summarise(median_D4_D7=median(log2FoldChange_D4_D7),median_D4_D10=median(log2FoldChange_D4_D10),median_D4_D14=median(log2FoldChange_D4_D14),median_D4_D21=median(log2FoldChange_D4_D21), median_continuous=median(log2FoldChange_continuous))
  median_scaled<- summary %>% mutate(median_D4_D7=adj$median_D4_D7) %>% mutate(median_D4_D10=adj$median_D4_D10)%>% mutate(median_D4_D14=adj$median_D4_D14)%>% mutate(median_D4_D21=adj$median_D4_D21) %>% mutate(median_continuous=adj$median_continuous) %>% mutate(adj_lfc_D4_D7=log2FoldChange_D4_D7-median_D4_D7) %>% mutate(adj_lfc_D4_D10=log2FoldChange_D4_D10-median_D4_D10)%>% mutate(adj_lfc_D4_D14=log2FoldChange_D4_D14-median_D4_D14) %>% mutate(adj_lfc_D4_D21=log2FoldChange_D4_D21-median_D4_D21) %>% mutate(adj_lfc_continuous=log2FoldChange_continuous-median_continuous) %>% mutate(adj_score_D4_D7=adj_lfc_D4_D7/lfcSE_D4_D7) %>% mutate(adj_score_D4_D10=adj_lfc_D4_D10/lfcSE_D4_D10)%>% mutate(adj_score_D4_D14=adj_lfc_D4_D14/lfcSE_D4_D14) %>% mutate(adj_score_D4_D21=adj_lfc_D4_D21/lfcSE_D4_D21) %>% mutate(adj_score_continuous=adj_lfc_continuous/lfcSE_continuous)                              
  median_scaled<-median_scaled %>% select("mseq","median_D4_D7","median_D4_D10","median_D4_D14", "median_D4_D21", "median_continuous","adj_lfc_D4_D7", "adj_lfc_D4_D10", "adj_lfc_D4_D14", "adj_lfc_D4_D21", "adj_lfc_continuous", "adj_score_D4_D7", "adj_score_D4_D10", "adj_score_D4_D14", "adj_score_D4_D21", "adj_score_continuous")
  
  
  #recalculate p values based on z score produced from median scaling, standard error does not need to be recalculated after the scaling as this is a simple translation (ie moving the y axis up and down)
  median_scaled<-median_scaled %>% mutate(uncombined_two_tailed_p_D4_D7= pnorm(abs(adj_score_D4_D7),lower.tail = FALSE) *2) %>% mutate(uncombined_BH_FDR_D4_D7 = p.adjust(uncombined_two_tailed_p_D4_D7, method = "BH")) %>% 
    mutate(uncombined_two_tailed_p_D4_D10= pnorm(abs(adj_score_D4_D10),lower.tail = FALSE) *2) %>% mutate(uncombined_BH_FDR_D4_D10 = p.adjust(uncombined_two_tailed_p_D4_D10, method = "BH")) %>%
    mutate(uncombined_two_tailed_p_D4_D14= pnorm(abs(adj_score_D4_D14),lower.tail = FALSE) *2) %>% mutate(uncombined_BH_FDR_D4_D14 = p.adjust(uncombined_two_tailed_p_D4_D14, method = "BH")) %>%
    mutate(uncombined_two_tailed_p_D4_D21= pnorm(abs(adj_score_D4_D21),lower.tail = FALSE) *2) %>% mutate(uncombined_BH_FDR_D4_D21 = p.adjust(uncombined_two_tailed_p_D4_D21, method = "BH")) %>%
    mutate(uncombined_two_tailed_p_continuous= pnorm(abs(adj_score_continuous),lower.tail = FALSE) *2) %>% mutate(uncombined_BH_FDR_continuous = p.adjust(uncombined_two_tailed_p_continuous, method = "BH"))
  summary<-merge(summary, median_scaled, by="mseq", all.x=TRUE)
  #merge summary to annotation file - this is the final output file
  annotation_file <-paste0(ident,"_master_annotated")
  annotation_file <-get(annotation_file)
  OUT<-summary
  OUT$SG<-sg_assay
  OUT$IDENT<-ident
  OUT$EXON_GROUP<-exon_assay
  OUT<-merge(OUT, annotation_file[ , c("id","mseq","pam_mut_sgrna_id")],  by="mseq", all.x=TRUE)
  OUT$combined_name<-paste0(OUT$id,"_",OUT$EXON_GROUP)
  OUT<-OUT %>% relocate(SG, .after = mseq) %>% relocate(IDENT, .after = SG) %>% relocate(EXON_GROUP, .after = IDENT) %>% relocate(combined_name, .after = mseq) %>% relocate(pam_mut_sgrna_id, .after = EXON_GROUP) %>% relocate(id, .after = mseq)
  #make multiple outputs - one full like before and one input for the combination
  OUT_full_anot<-merge(OUT, annotation_file,  by="id", all.x=TRUE)
  names(OUT_full_anot)[names(OUT_full_anot) == 'pam_mut_sgrna_id.x'] <- 'pam_mut_sgrna_id'
  names(OUT_full_anot)[names(OUT_full_anot) == 'mseq.x'] <- 'mseq'
  OUT_full_anot$pam_mut_sgrna_id.y<-NULL
  OUT_full_anot$mseq.y<-NULL
  write.csv(OUT, file=(paste0("./",ident,"_OUT.csv")), row.names = FALSE)
  write.csv(OUT_full_anot, file=(paste0("./",ident,"_OUT_full_anot.csv")), row.names = FALSE)
  #name and output files to global environment
  #deseqname<-deparse(substitute(x_filtered_counts))
  deseqname<-paste0(ident,"_OUT")
  deseqname_full<-paste0(ident,"_OUT_full_anot")
  assign(deseqname, OUT, envir = .GlobalEnv)
  assign(deseqname_full, OUT_full_anot, envir = .GlobalEnv)
}

#########################################################################
#FUNCTIONS_DESEQ_END
#########################################################################

#########################################################################
#RUN ALL OF THE ANALYSIS AND QC PLOTS OUTPUT START
#########################################################################

#Function to apply to listed exons and designated sgRNA_A
lapply (c(1,2,3,4,5,6,7,9,10,"11_1","11_2","12_1","12_2","13_1","13_2","13_3",14,15,16,"17_1","17_2"), function(xx) {run.deseq(exon_assay=xx,sg_assay="A")
})
#Function to apply to listed exons and designated sgRNA_B
lapply (c(1,2,4,5,6,7,8,9,10,"11_1","11_2","12_1","12_2","13_1","13_2","13_3",14,15,16,"17_1","17_2"), function(xx) {run.deseq(exon_assay=xx,sg_assay="B")
})
#########################################################################
#RUN ALL OF THE ANALYSIS AND QC PLOTS OUTPUT END
#########################################################################


#########################################################################
#RUNS THAT REQUIRE COMPLETE EXLCUSION- START 
#########################################################################

#BAP1: exclude following replicates:
#All of E2_SGB - POOR HDR RATE AND ALSO REPLICATE CONSISTENCY ESCPECIALLY IN: E2_SGB_D7R1; E2_SGB_D10R1; E2_SGB_D14R1; E2_SGB_D21R1
#All of E13_1_SGB - POOR HDR RATE, STONG POSITIONAL EFFECT AND ALSO REPLICATE CONSISTENCY ESCPECIALLY IN: E13_1_SGB_D7R3; E13_1_SGB_D10R3; E13_1_SGB_D14R3; E13_1_SGB_D21R3

#########################################################################
#RUNS THAT REQUIRE COMPLETE EXLCUSION- END
#########################################################################


#########################################################################
#RUNS THAT REQUIRE REPLICATE EXLCUSION- START 
#########################################################################

##############################################################################
#RUN REPEATS THAT REQUIRE REPLICATE EXCLUSION - START - EXON 12_1_SGA START --------------------------------------------------------------------------
##############################################################################

#E12_1_SGA_D7R3; E12_1_SGA_D10R3
#setwd to sub directory './exclude')
annot_a<-read.csv("/Users/aw28/Documents/BAP1_analysis/experimental_design/annot_a_exclude_12.csv")
annot_a_continuous<-read.csv("/Users/aw28/Documents/BAP1_analysis/experimental_design/annot_a_continuous_exclude_12.csv")
annot_a<-as.matrix.data.frame(annot_a)
annot_a_continuous$condition<-as.numeric(annot_a_continuous$condition)

#MAKE NEW NORMALIZATION MATRIX THAT HAS UNDESIRED REPLICATES EXLUDED
#make normalization matricies - EXCLUDED
deseq.norm.input <- function(exon, sg){
  x_normalization_counts <- paste0("E",exon,"_SG",sg,"_normalization_counts")
  x_normalization_counts <-get(x_normalization_counts)
  deseq_input_noramlization_counts <- x_normalization_counts %>% dplyr::select("mseq", "D4R1", "D4R2", "D4R3", "D7R1", "D7R2", "D10R1", "D10R2",  "D14R1", "D14R2",  "D14R3",  "D21R1", "D21R2",  "D21R3")
  deseq_input_noramlization_counts <- deseq_input_noramlization_counts %>% remove_rownames %>% column_to_rownames(var="mseq")
  deseq_input_noramlization_counts <- as.matrix.data.frame(deseq_input_noramlization_counts)
  normdfname<-paste0("E",exon,"_SG",sg,"_normalization_counts_MATRIX")
  assign(normdfname, deseq_input_noramlization_counts, envir = .GlobalEnv)
}
#Function to apply the above function to listed exons and designated sgRNA_A
lapply (c("12_1"), function(y) {deseq.norm.input(exon=y, sg="A")
})

#MAKE NEW COUNT MATRIX THAT HAS UNDESIRED REPLICATES EXLUDED
#make count matricies
deseq.count.input <- function(exon,sg){
  x_filtered_counts <- paste0("E",exon,"_SG",sg,"_master_annotated")
  x_filtered_counts <- get(x_filtered_counts)
  #line below removes duplicates within the targeton file based on mseq - some mutators produce the same mseq but oligo library contains one instance
  x_filtered_counts <-x_filtered_counts[!duplicated(x_filtered_counts$mseq), ]
  deseq_input_filtered_counts <- x_filtered_counts %>% dplyr::select("mseq", "D4R1", "D4R2", "D4R3", "D7R1", "D7R2", "D10R1", "D10R2", "D14R1", "D14R2",  "D14R3", "D21R1", "D21R2", "D21R3")
  deseq_input_filtered_counts <- deseq_input_filtered_counts %>% remove_rownames %>% column_to_rownames(var="mseq")
  deseq_input_filtered_counts <- as.matrix.data.frame(deseq_input_filtered_counts)
  countdfname<-paste0("E",exon,"_SG",sg,"_master_annotated_filtered_counts_MATRIX")
  assign(countdfname, deseq_input_filtered_counts, envir = .GlobalEnv)
}
#Function to apply the above function to listed exons and designated sgRNA_B
lapply (c("12_1"), function(y) {deseq.count.input(exon=y, sg="A")
})

#RUN DESEQ AGAIN WITH THE EXLCUDED COUNT AND NORMALIZATION MATRICIES AND NEW EXPERIMENTAL DESIGN
lapply (c("12_1"), function(xx) {run.deseq(exon_assay=xx,sg_assay="A")
})

##############################################################################
#RUN REPEATS THAT REQUIRE REPLICATE EXCLUSION - START - EXON 12_1_SGA END --------------------------------------------------------------------------
##############################################################################

#########################################################################
#RUN REPEATS THAT REQUIRE REPLICATE EXCLUSION - END 
#########################################################################


###############----------------------------------------------------###########################################################
###############-----------DESEQ2_END-------------------------------###########################################################
###############----------------------------------------------------###########################################################

#the above DESeq2 scripts use unique mseq per targeton - mseqs will be duplicated if not made unique as. This is because different 'id' are not unique for sequence 'mseq'.
#mutator function used by valiant is part of the 'id' nomenclature, so the same variant at the 'mseq' level is produced by varied means and will produce a duplicated mseq at the annotation stage
#the mseqs are expanded again at the final annotation stage in the deseq process and all ids will again be present in the below bindings

#This binds together all the outputs where there are two guide RNA libraries that have been processed successfully 
#First we bind A libraries together
total_combined_bap1_a <- do.call("rbind.fill", list(E1_SGA_OUT,
                                                    #E2_SGA_OUT,
                                                    #E3_SGA_OUT,
                                                    E4_SGA_OUT,
                                                    E5_SGA_OUT,
                                                    E6_SGA_OUT,
                                                    E7_SGA_OUT,
                                                    E9_SGA_OUT,
                                                    E10_SGA_OUT,
                                                    E11_1_SGA_OUT,
                                                    E11_2_SGA_OUT,
                                                    E12_1_SGA_OUT,
                                                    E12_2_SGA_OUT,
                                                    #E13_1_SGA_OUT,
                                                    E13_2_SGA_OUT,
                                                    E13_3_SGA_OUT,
                                                    E14_SGA_OUT,
                                                    E15_SGA_OUT,
                                                    E16_SGA_OUT,
                                                    E17_1_SGA_OUT,
                                                    E17_2_SGA_OUT))
#Now bind B libraries together
total_combined_bap1_b <- do.call("rbind.fill", list(E1_SGB_OUT,
                                                    #E2_SGB_OUT,
                                                    E4_SGB_OUT,
                                                    E5_SGB_OUT,
                                                    E6_SGB_OUT,
                                                    E7_SGB_OUT,
                                                    E9_SGB_OUT,
                                                    #E8_SGB_OUT,
                                                    E10_SGB_OUT,
                                                    E11_1_SGB_OUT,
                                                    E11_2_SGB_OUT,
                                                    E12_1_SGB_OUT,
                                                    E12_2_SGB_OUT,
                                                    #E13_1_SGB_OUT,
                                                    E13_2_SGB_OUT,
                                                    E13_3_SGB_OUT,
                                                    E14_SGB_OUT,
                                                    E15_SGB_OUT,
                                                    E16_SGB_OUT,
                                                    E17_1_SGB_OUT,
                                                    E17_2_SGB_OUT))
#can bind together libraries - the 'id' field produced by valiant in these dataframes are guide library agnostic (ie. the duplicate variant will will have the same id)
total_bap1 <- do.call("rbind.fill", list(total_combined_bap1_a,total_combined_bap1_b))
#This binds the targetons that do not have two guide RNA libraries
total_single_bap1 <- do.call("rbind.fill", list(E2_SGA_OUT,
                                                E3_SGA_OUT,
                                                E13_1_SGA_OUT,
                                                E8_SGB_OUT))
#can use the below to write out the single targeton data
#write.csv(total_single_bap1,"./single_targetons.csv", row.names = FALSE)

#this produces a data frame with all variants processed - duplicated
total_unique_id_count<-do.call("rbind.fill", list(total_combined_bap1_a,total_combined_bap1_b,total_single_bap1))
#this produces a data frame with all variants processed - not duplicated for 'id' but will be duplicated for mseq
total_unique_id_count<-total_unique_id_count[!duplicated(total_unique_id_count$id), ]
#this produces a data frame with all variants processed - not duplicated for 'mseq' - what the experimental flasks would have contained as collated unique libraries
#total_unique_id_count<-total_unique_id_count[!duplicated(total_unique_id_count$mseq), ]

total_bap1_pam_full<-total_bap1 %>%
  group_by(mseq) %>%
  dplyr::summarize(dup_mseq= n(), 
                   PAM_status_full_arrayed = paste(pam_mut_sgrna_id, collapse = ','))%>%
  mutate(PAM_status_simplified=case_when(str_detect(PAM_status_full_arrayed,"_b")~ "B", TRUE~ "-")) %>%
  mutate(PAM_status_simplified=case_when(str_detect(PAM_status_full_arrayed,"_a")~ "A", TRUE~ PAM_status_simplified)) %>%
  select(c(1,2,4))

total_bap1_pam_corrected<-merge(total_bap1,total_bap1_pam_full,by="mseq",all.x=TRUE)

total_bap1_pam_corrected<-total_bap1_pam_corrected %>% arrange(EXON_GROUP,id,SG)


total_single_bap1_pam_full<-total_single_bap1 %>%
  group_by(mseq) %>%
  dplyr::summarize(dup_mseq= n(), 
                   PAM_status_full_arrayed = paste(pam_mut_sgrna_id, collapse = ','))%>%
  mutate(PAM_status_simplified=case_when(str_detect(PAM_status_full_arrayed,"_b")~ "B", TRUE~ "-")) %>%
  mutate(PAM_status_simplified=case_when(str_detect(PAM_status_full_arrayed,"_a")~ "A", TRUE~ PAM_status_simplified)) %>%
  select(c(1,2,4))

total_single_bap1_pam_corrected<-merge(total_single_bap1,total_single_bap1_pam_full,by="mseq",all.x=TRUE)

total_single_bap1_pam_corrected<-total_single_bap1_pam_corrected  %>% arrange(EXON_GROUP,id,SG)


################# COMBINATION OF GUIDES ############# START
################# COMBINATION OF GUIDES ############# START
################# COMBINATION OF GUIDES ############# START
#in this step we take only the libraries where there is a possibility of library duplication (ie. the same variant assayed with guide library a and library b)
#we group by the exon group in which the variant is found - so there will be duplicated ids as a whole
#fields of interest for each duplicated variant are put into a new field and separated by comma then key fields for calculation are ungrouped by comma into new fields
################# COMBINATION OF GUIDES ############# START
################# COMBINATION OF GUIDES ############# START
################# COMBINATION OF GUIDES ############# START


sga_sgb_collapse <- total_bap1_pam_corrected  %>% 
  group_by(id,EXON_GROUP) %>%
  dplyr::summarize(Variant_duplication= n(), 
                   Variant_Sources = paste(SG, collapse = ','),
                   PAM_status = paste(pam_mut_sgrna_id, collapse = ','),
                   PAM_status_simplified_collated = paste(PAM_status_simplified, collapse = ','),
                   dup_mseq_collated = paste(dup_mseq, collapse = ','),
                   mseq_combined = paste(mseq, collapse = ','),
                   baseMean = paste(baseMean, collapse = ','), 
                   raw_LFC_D4_D7 = paste(log2FoldChange_D4_D7, collapse = ','), 
                   raw_LFC_D4_D10 = paste(log2FoldChange_D4_D10, collapse = ','), 
                   raw_LFC_D4_D14 = paste(log2FoldChange_D4_D14, collapse = ','), 
                   raw_LFC_D4_D21 = paste(log2FoldChange_D4_D21, collapse = ','),
                   raw_LFC_continuous = paste(log2FoldChange_continuous, collapse = ','),
                   
                   lfcSE_D4_D7 = paste(lfcSE_D4_D7, collapse = ','), 
                   lfcSE_D4_D10 = paste(lfcSE_D4_D10, collapse = ','), 
                   lfcSE_D4_D14 = paste(lfcSE_D4_D14, collapse = ','), 
                   lfcSE_D4_D21 = paste(lfcSE_D4_D21, collapse = ','), 
                   lfcSE_continuous = paste(lfcSE_continuous, collapse = ','), 
                   
                   adj_lfc_D4_D7 = paste(adj_lfc_D4_D7, collapse = ','),
                   adj_lfc_D4_D10 = paste(adj_lfc_D4_D10, collapse = ','),
                   adj_lfc_D4_D14 = paste(adj_lfc_D4_D14, collapse = ','),
                   adj_lfc_D4_D21 = paste(adj_lfc_D4_D21, collapse = ','),
                   adj_lfc_continuous = paste(adj_lfc_continuous, collapse = ',')) %>%
  
  ungroup() %>%
  separate(baseMean, sep=",", c("sga_baseMean","sgb_baseMean")) %>%
  separate(mseq_combined, sep=",", c("mseq_a","mseq_b")) %>%
  separate(adj_lfc_D4_D7, sep=",", c("sga_adj_lfc_D4_D7","sgb_adj_lfc_D4_D7")) %>%
  separate(adj_lfc_D4_D10, sep=",", c("sga_adj_lfc_D4_D10","sgb_adj_lfc_D4_D10")) %>%
  separate(adj_lfc_D4_D14, sep=",", c("sga_adj_lfc_D4_D14","sgb_adj_lfc_D4_D14")) %>%
  separate(adj_lfc_D4_D21, sep=",", c("sga_adj_lfc_D4_D21","sgb_adj_lfc_D4_D21")) %>%
  separate(adj_lfc_continuous, sep=",", c("sga_adj_lfc_continuous","sgb_adj_lfc_continuous")) %>%
  
  separate(raw_LFC_D4_D7, sep=",", c("sga_raw_LFC_D4_D7","sgb_raw_LFC_D4_D7")) %>%
  separate(raw_LFC_D4_D10, sep=",", c("sga_raw_LFC_D4_D10","sgb_raw_LFC_D4_D10")) %>%
  separate(raw_LFC_D4_D14, sep=",", c("sga_raw_LFC_D4_D14","sgb_raw_LFC_D4_D14")) %>%
  separate(raw_LFC_D4_D21, sep=",", c("sga_raw_LFC_D4_D21","sgb_raw_LFC_D4_D21")) %>%
  separate(raw_LFC_continuous, sep=",", c("sga_raw_LFC_continuous","sgb_raw_LFC_continuous")) %>%
  
  separate(lfcSE_D4_D7, sep=",", c("sga_lfcSE_D4_D7","sgb_lfcSE_D4_D7")) %>%
  separate(lfcSE_D4_D10, sep=",", c("sga_lfcSE_D4_D10","sgb_lfcSE_D4_D10"))%>%
  separate(lfcSE_D4_D14, sep=",", c("sga_lfcSE_D4_D14","sgb_lfcSE_D4_D14"))%>%
  separate(lfcSE_D4_D21, sep=",", c("sga_lfcSE_D4_D21","sgb_lfcSE_D4_D21")) %>%
  separate(lfcSE_continuous, sep=",", c("sga_lfcSE_continuous","sgb_lfcSE_continuous")) %>%
  
  #obtain values for instances where a variant is specific to one library - basemean
  mutate(sgb_baseMean=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ sga_baseMean, TRUE ~ sgb_baseMean)) %>%
  mutate(sga_baseMean=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ "NA", TRUE ~ sga_baseMean)) %>%
  #obtain values for instances where a variant is specific to one library - mseq (the oligo used to generate the variant)
  mutate(mseq_b=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ mseq_a, TRUE ~ mseq_b)) %>%
  mutate(mseq_a=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ "NA", TRUE ~ mseq_a)) %>%
  
  #obtain values for instances where a variant is specific to one library - median scaled LFCs (that is the median LFC of intron and synonymous variants for that targeton subtracted post-deseq2)
  #these values will not be used to calulate a weighted mean as there is only one variant instance, there will be no weight added due to error and the median scaled values will be progressed
  mutate(sgb_adj_lfc_D4_D7=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ sga_adj_lfc_D4_D7, TRUE ~ sgb_adj_lfc_D4_D7)) %>%
  mutate(sga_adj_lfc_D4_D7=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ "NA", TRUE ~ sga_adj_lfc_D4_D7)) %>%
  mutate(sgb_adj_lfc_D4_D10=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ sga_adj_lfc_D4_D10, TRUE ~ sgb_adj_lfc_D4_D10)) %>%
  mutate(sga_adj_lfc_D4_D10=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ "NA", TRUE ~ sga_adj_lfc_D4_D10)) %>%
  mutate(sgb_adj_lfc_D4_D14=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ sga_adj_lfc_D4_D14, TRUE ~ sgb_adj_lfc_D4_D14)) %>%
  mutate(sga_adj_lfc_D4_D14=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ "NA", TRUE ~ sga_adj_lfc_D4_D14)) %>%
  mutate(sgb_adj_lfc_D4_D21=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ sga_adj_lfc_D4_D21, TRUE ~ sgb_adj_lfc_D4_D21)) %>%
  mutate(sga_adj_lfc_D4_D21=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ "NA", TRUE ~ sga_adj_lfc_D4_D21)) %>%
  mutate(sgb_adj_lfc_continuous=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ sga_adj_lfc_continuous, TRUE ~ sgb_adj_lfc_continuous)) %>%
  mutate(sga_adj_lfc_continuous=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ "NA", TRUE ~ sga_adj_lfc_continuous)) %>%
  
  #obtain values for instances where a variant is specific to one library - standard error from DESeq2 
  mutate(sgb_lfcSE_D4_D7=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ sga_lfcSE_D4_D7, TRUE ~ sgb_lfcSE_D4_D7)) %>%
  mutate(sga_lfcSE_D4_D7=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ "NA", TRUE ~ sga_lfcSE_D4_D7)) %>%
  mutate(sgb_lfcSE_D4_D10=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ sga_lfcSE_D4_D10, TRUE ~ sgb_lfcSE_D4_D10)) %>%
  mutate(sga_lfcSE_D4_D10=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ "NA", TRUE ~ sga_lfcSE_D4_D10)) %>%
  mutate(sgb_lfcSE_D4_D14=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ sga_lfcSE_D4_D14, TRUE ~ sgb_lfcSE_D4_D14)) %>%
  mutate(sga_lfcSE_D4_D14=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ "NA", TRUE ~ sga_lfcSE_D4_D14)) %>%
  mutate(sgb_lfcSE_D4_D21=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ sga_lfcSE_D4_D21, TRUE ~ sgb_lfcSE_D4_D21)) %>%
  mutate(sga_lfcSE_D4_D21=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ "NA", TRUE ~ sga_lfcSE_D4_D21)) %>%
  mutate(sgb_lfcSE_continuous=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ sga_lfcSE_continuous, TRUE ~ sgb_lfcSE_continuous)) %>%
  mutate(sga_lfcSE_continuous=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ "NA", TRUE ~ sga_lfcSE_continuous)) %>%
  
  #obtain values for instances where a variant is specific to one library - raw log2fold changes (before median scaling), not used in calculation but to preserve provenance of value
  mutate(sgb_raw_LFC_D4_D7=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ sga_raw_LFC_D4_D7, TRUE ~ sgb_raw_LFC_D4_D7)) %>%
  mutate(sga_raw_LFC_D4_D7=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ "NA", TRUE ~ sga_raw_LFC_D4_D7)) %>%
  mutate(sgb_raw_LFC_D4_D10=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ sga_raw_LFC_D4_D10, TRUE ~ sgb_raw_LFC_D4_D10)) %>%
  mutate(sga_raw_LFC_D4_D10=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ "NA", TRUE ~ sga_raw_LFC_D4_D10)) %>%
  mutate(sgb_raw_LFC_D4_D14=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ sga_raw_LFC_D4_D14, TRUE ~ sgb_raw_LFC_D4_D14)) %>%
  mutate(sga_raw_LFC_D4_D14=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ "NA", TRUE ~ sga_raw_LFC_D4_D14)) %>%
  mutate(sgb_raw_LFC_D4_D21=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ sga_raw_LFC_D4_D21, TRUE ~ sgb_raw_LFC_D4_D21)) %>%
  mutate(sga_raw_LFC_D4_D21=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ "NA", TRUE ~ sga_raw_LFC_D4_D21)) %>%
  mutate(sgb_raw_LFC_continuous=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ sga_raw_LFC_continuous, TRUE ~ sgb_raw_LFC_continuous)) %>%
  mutate(sga_raw_LFC_continuous=case_when(Variant_duplication=="1" & Variant_Sources == "B" ~ "NA", TRUE ~ sga_raw_LFC_continuous)) #%>%


#make key field contents nurmeric for further calculations
sga_sgb_collapse[, c(10:41)] <- sapply(sga_sgb_collapse[, c(10:41)], as.numeric)

#perform calculations to obtained weighted means -if a duplicated variant is in a pam codon, the weight for that library contribution will be 0 and the alternative library will contribute 100% to the combined value
sga_sgb_combined2 <-sga_sgb_collapse %>%
  mutate(weight_a_D4_D7=case_when(str_detect(PAM_status_simplified_collated,"A")~ 0, TRUE~ 1/(sga_lfcSE_D4_D7)^2)) %>% mutate(weight_b_D4_D7= case_when(str_detect(PAM_status_simplified_collated,"B")~ 0, TRUE~1/(sgb_lfcSE_D4_D7)^2)) %>%
  mutate(weight_a_D4_D10=case_when(str_detect(PAM_status_simplified_collated,"A")~ 0, TRUE~ 1/(sga_lfcSE_D4_D10)^2)) %>% mutate(weight_b_D4_D10= case_when(str_detect(PAM_status_simplified_collated,"B")~ 0, TRUE~1/(sgb_lfcSE_D4_D10)^2)) %>%
  mutate(weight_a_D4_D14=case_when(str_detect(PAM_status_simplified_collated,"A")~ 0, TRUE~ 1/(sga_lfcSE_D4_D14)^2)) %>% mutate(weight_b_D4_D14= case_when(str_detect(PAM_status_simplified_collated,"B")~ 0, TRUE~1/(sgb_lfcSE_D4_D14)^2)) %>%
  mutate(weight_a_D4_D21=case_when(str_detect(PAM_status_simplified_collated,"A")~ 0, TRUE~ 1/(sga_lfcSE_D4_D21)^2)) %>% mutate(weight_b_D4_D21= case_when(str_detect(PAM_status_simplified_collated,"B")~ 0, TRUE~1/(sgb_lfcSE_D4_D21)^2)) %>%
  mutate(weight_a_continuous=case_when(str_detect(PAM_status_simplified_collated,"A")~ 0, TRUE~ 1/(sga_lfcSE_continuous)^2)) %>% mutate(weight_b_continuous= case_when(str_detect(PAM_status_simplified_collated,"B")~ 0, TRUE~1/(sgb_lfcSE_continuous)^2)) %>%
  
  mutate(sum_of_weight_D4_D7=rowSums(cbind(weight_a_D4_D7,weight_b_D4_D7), na.rm=TRUE))%>%
  mutate(sum_of_weight_D4_D10=rowSums(cbind(weight_a_D4_D10,weight_b_D4_D10), na.rm=TRUE))%>%
  mutate(sum_of_weight_D4_D14=rowSums(cbind(weight_a_D4_D14,weight_b_D4_D14), na.rm=TRUE))%>%
  mutate(sum_of_weight_D4_D21=rowSums(cbind(weight_a_D4_D21,weight_b_D4_D21), na.rm=TRUE))%>%
  mutate(sum_of_weight_continuous=rowSums(cbind(weight_a_continuous,weight_b_continuous), na.rm=TRUE))%>%
  
  mutate(SE_bind_D4_D7 = (sum_of_weight_D4_D7)^(-0.5)) %>%
  mutate(SE_bind_D4_D10 = (sum_of_weight_D4_D10)^(-0.5)) %>%
  mutate(SE_bind_D4_D14 = (sum_of_weight_D4_D14)^(-0.5)) %>%
  mutate(SE_bind_D4_D21 = (sum_of_weight_D4_D21)^(-0.5)) %>%
  mutate(SE_bind_continuous = (sum_of_weight_continuous)^(-0.5)) %>%
  
  mutate(weighted_sga_LFC_D4_D7= weight_a_D4_D7*sga_adj_lfc_D4_D7) %>% mutate(weighted_sgb_LFC_D4_D7= weight_b_D4_D7*sgb_adj_lfc_D4_D7) %>%
  mutate(weighted_sga_LFC_D4_D10= weight_a_D4_D10*sga_adj_lfc_D4_D10) %>% mutate(weighted_sgb_LFC_D4_D10= weight_b_D4_D10*sgb_adj_lfc_D4_D10) %>%
  mutate(weighted_sga_LFC_D4_D14= weight_a_D4_D14*sga_adj_lfc_D4_D14) %>% mutate(weighted_sgb_LFC_D4_D14= weight_b_D4_D14*sgb_adj_lfc_D4_D14) %>%
  mutate(weighted_sga_LFC_D4_D21= weight_a_D4_D21*sga_adj_lfc_D4_D21) %>% mutate(weighted_sgb_LFC_D4_D21= weight_b_D4_D21*sgb_adj_lfc_D4_D21) %>%
  mutate(weighted_sga_LFC_continuous= weight_a_continuous*sga_adj_lfc_continuous) %>% mutate(weighted_sgb_LFC_continuous= weight_b_continuous*sgb_adj_lfc_continuous) %>%
  
  mutate(sum_of_weighted_LFC_D4_D7=rowSums(cbind(weighted_sga_LFC_D4_D7,weighted_sgb_LFC_D4_D7), na.rm=TRUE)) %>% mutate(combined_LFC_D4_D7=sum_of_weighted_LFC_D4_D7/sum_of_weight_D4_D7) %>% filter(combined_LFC_D4_D7 != "NaN") %>%
  mutate(sum_of_weighted_LFC_D4_D10=rowSums(cbind(weighted_sga_LFC_D4_D10,weighted_sgb_LFC_D4_D10), na.rm=TRUE)) %>% mutate(combined_LFC_D4_D10=sum_of_weighted_LFC_D4_D10/sum_of_weight_D4_D10) %>% filter(combined_LFC_D4_D10 != "NaN") %>%
  mutate(sum_of_weighted_LFC_D4_D14=rowSums(cbind(weighted_sga_LFC_D4_D14,weighted_sgb_LFC_D4_D14), na.rm=TRUE)) %>% mutate(combined_LFC_D4_D14=sum_of_weighted_LFC_D4_D14/sum_of_weight_D4_D14) %>% filter(combined_LFC_D4_D14 != "NaN") %>%
  mutate(sum_of_weighted_LFC_D4_D21=rowSums(cbind(weighted_sga_LFC_D4_D21,weighted_sgb_LFC_D4_D21), na.rm=TRUE)) %>% mutate(combined_LFC_D4_D21=sum_of_weighted_LFC_D4_D21/sum_of_weight_D4_D21) %>% filter(combined_LFC_D4_D21 != "NaN") %>%
  mutate(sum_of_weighted_LFC_continuous=rowSums(cbind(weighted_sga_LFC_continuous,weighted_sgb_LFC_continuous), na.rm=TRUE)) %>% mutate(combined_LFC_continuous=sum_of_weighted_LFC_continuous/sum_of_weight_continuous) %>% filter(combined_LFC_continuous != "NaN") %>%
  
  #can use the calculated z scores to obtain a new p value - however at this stage there are mseq duplications, in other words variants with different ids but with the same mseq and as
  #twist libraries are made to the unique for mseq there will be value duplications in the z distribution, which might affect the p value calculation and will definitely affect BH correction
  #so need to wait to the obtain the final p-values and corrected p-values in the condensed dataset that 
  mutate(combined_Z_D4_D7=combined_LFC_D4_D7/SE_bind_D4_D7) %>% mutate(two_tailed_p_D4_D7= pnorm(abs(combined_Z_D4_D7),lower.tail = FALSE) *2) %>% #mutate(BH_FDR_D4_D7 = p.adjust(two_tailed_p_D4_D7, method = "BH")) %>%
  mutate(combined_Z_D4_D10=combined_LFC_D4_D10/SE_bind_D4_D10) %>% mutate(two_tailed_p_D4_D10= pnorm(abs(combined_Z_D4_D10),lower.tail = FALSE) *2) %>% #mutate(BH_FDR_D4_D10 = p.adjust(two_tailed_p_D4_D10, method = "BH")) %>%
  mutate(combined_Z_D4_D14=combined_LFC_D4_D14/SE_bind_D4_D14) %>% mutate(two_tailed_p_D4_D14= pnorm(abs(combined_Z_D4_D14),lower.tail = FALSE) *2) %>% #mutate(BH_FDR_D4_D14 = p.adjust(two_tailed_p_D4_D14, method = "BH")) %>%
  mutate(combined_Z_D4_D21=combined_LFC_D4_D21/SE_bind_D4_D21) %>% mutate(two_tailed_p_D4_D21= pnorm(abs(combined_Z_D4_D21),lower.tail = FALSE) *2) %>% #mutate(BH_FDR_D4_D21 = p.adjust(two_tailed_p_D4_D21, method = "BH"))
  mutate(combined_Z_continuous=combined_LFC_continuous/SE_bind_continuous) %>% mutate(two_tailed_p_continuous= pnorm(abs(combined_Z_continuous),lower.tail = FALSE) *2) #%>% mutate(BH_FDR_continuous = p.adjust(two_tailed_p_continuous, method = "BH"))

write.csv(sga_sgb_combined2,"./sga_sgb_combined.csv", row.names = FALSE)


################# COMBINATION OF GUIDES ############# END
################# COMBINATION OF GUIDES ############# END
################# COMBINATION OF GUIDES ############# END
################# COMBINATION OF GUIDES ############# END
################# COMBINATION OF GUIDES ############# END
################# COMBINATION OF GUIDES ############# END

################# PREP FOR COMBINATION OF TARGETONS ############# START
################# PREP FOR COMBINATION OF TARGETONS ############# START
################# PREP FOR COMBINATION OF TARGETONS ############# START
################# PREP FOR COMBINATION OF TARGETONS ############# START
################# PREP FOR COMBINATION OF TARGETONS ############# START
################# PREP FOR COMBINATION OF TARGETONS ############# START

#produced dataset with reduced number of columns so that combined guides and uncombined guides datasets can be stacked together
combined_sub <- sga_sgb_combined2[ ,c(1:9,57:61,73,75,77,79,81:91)]
#make new field that allows identification as to whether the values were obtained from combined LFC or from a single, uncombined, targeton
combined_sub$dataset_process<-"combined"

#change the name of the logfold change standard error that was the product of a weighted process in the combination calculations into the lfcSE format given by Deseq2, as the uncombined data labelled this way
combined_sub<-combined_sub %>% dplyr::rename(lfcSE_D4_D7 = "SE_bind_D4_D7")
combined_sub<-combined_sub %>% dplyr::rename(lfcSE_D4_D10 = "SE_bind_D4_D10")
combined_sub<-combined_sub %>% dplyr::rename(lfcSE_D4_D14 = "SE_bind_D4_D14")
combined_sub<-combined_sub %>% dplyr::rename(lfcSE_D4_D21 = "SE_bind_D4_D21")
combined_sub<-combined_sub %>% dplyr::rename(lfcSE_continuous = "SE_bind_continuous")
#change the name of the combined LFCs which will be the combination of two guides adjusted log2fold change, into the Deseq2 process script output
combined_sub<-combined_sub %>% dplyr::rename(adj_lfc_D4_D7 = "combined_LFC_D4_D7")
combined_sub<-combined_sub %>% dplyr::rename(adj_lfc_D4_D10 = "combined_LFC_D4_D10")
combined_sub<-combined_sub %>% dplyr::rename(adj_lfc_D4_D14 = "combined_LFC_D4_D14")
combined_sub<-combined_sub %>% dplyr::rename(adj_lfc_D4_D21 = "combined_LFC_D4_D21")
combined_sub<-combined_sub %>% dplyr::rename(adj_lfc_continuous = "combined_LFC_continuous")
#change the name of the combined z scores, to be the same as the uncombined columns
combined_sub<-combined_sub %>% dplyr::rename(adj_score_D4_D7 = "combined_Z_D4_D7")
combined_sub<-combined_sub %>% dplyr::rename(adj_score_D4_D10 = "combined_Z_D4_D10")
combined_sub<-combined_sub %>% dplyr::rename(adj_score_D4_D14 = "combined_Z_D4_D14")
combined_sub<-combined_sub %>% dplyr::rename(adj_score_D4_D21 = "combined_Z_D4_D21")
combined_sub<-combined_sub %>% dplyr::rename(adj_score_continuous = "combined_Z_continuous")

#select the necessary columns from the uncombined data
uncombined_sub <- total_single_bap1_pam_corrected
#add a column and populate with description of this data as being from an uncombined source, that is one guide was used to get this data
uncombined_sub$Variant_duplication<-NA
uncombined_sub$dataset_process<-"not_combined"
#change the names of the following columns, so that they can be bound to the combined dataframe
uncombined_sub<-uncombined_sub %>% dplyr::rename(Variant_Sources = "SG")
uncombined_sub<-uncombined_sub %>% dplyr::rename(PAM_status = "pam_mut_sgrna_id")
uncombined_sub<-uncombined_sub %>% dplyr::rename(PAM_status_simplified_collated = "PAM_status_simplified")
uncombined_sub<-uncombined_sub %>% dplyr::rename(dup_mseq_collated = "dup_mseq")
uncombined_sub<-uncombined_sub %>% dplyr::rename(two_tailed_p_D4_D7 = "uncombined_two_tailed_p_D4_D7")
uncombined_sub<-uncombined_sub %>% dplyr::rename(two_tailed_p_D4_D10 = "uncombined_two_tailed_p_D4_D10")
uncombined_sub<-uncombined_sub %>% dplyr::rename(two_tailed_p_D4_D14 = "uncombined_two_tailed_p_D4_D14")
uncombined_sub<-uncombined_sub %>% dplyr::rename(two_tailed_p_D4_D21 = "uncombined_two_tailed_p_D4_D21")
uncombined_sub<-uncombined_sub %>% dplyr::rename(two_tailed_p_continuous = "uncombined_two_tailed_p_continuous")

#turn mseq in uncombined sub into mutseq_a and mutseq_b in order to bind tables
uncombined_sub <- uncombined_sub %>% mutate(mseq_b=case_when(Variant_Sources == "B" ~ mseq, TRUE ~ "NA")) %>%
  mutate(mseq_a=case_when(Variant_Sources == "A" ~ mseq, TRUE ~ "NA"))
#remove the mseq field from the combined as it is now redundant
uncombined_sub$mseq<-NULL  

#now bind the combined dataframe (that is weighted data, based on library a and library b) and the uncombined dataframe (unweighted, library a or library b only data)
#only columns in present in both uncombined and combined will be non-NA completely
total_bap1_combined_uncombined <- do.call("rbind.fill", list(combined_sub,uncombined_sub))
#remove columns in combined that are not in uncombined
total_bap1_combined_uncombined<-total_bap1_combined_uncombined %>% select(c(1:30))

#make NAs in mseqb column text based NA for consistency
total_bap1_combined_uncombined<-total_bap1_combined_uncombined %>%
  mutate_at(vars(mseq_b), ~replace_na(., "NA")) 

#output the dataframe to check - this will be input dataframe to combined tiled regions
write.csv(total_bap1_combined_uncombined,"./total_bap1_combined_uncombined.csv", row.names = FALSE)


#can also run the section below to get a dataframe of duplciated ids, that will be within overlapping regions, from the sga sgb combination set = QC
#duplicated<-total_bap1_combined_uncombined %>%
#group_by(id) %>%
#filter(n() > 1) %>%
#ungroup

################# PREP FOR COMBINATION OF TARGETONS ############# END
################# PREP FOR COMBINATION OF TARGETONS ############# END
################# PREP FOR COMBINATION OF TARGETONS ############# END
################# PREP FOR COMBINATION OF TARGETONS ############# END
################# PREP FOR COMBINATION OF TARGETONS ############# END
################# PREP FOR COMBINATION OF TARGETONS ############# END



################# COMBINATION OF TARGETONS ############# START
################# COMBINATION OF TARGETONS ############# START
################# COMBINATION OF TARGETONS ############# START
################# COMBINATION OF TARGETONS ############# START
################# COMBINATION OF TARGETONS ############# START
################# COMBINATION OF TARGETONS ############# START

#group the targetons by id - this collapse the different mseq_a categories that assess the same variant from different (overlapping) targetons


targeton_input<-total_bap1_combined_uncombined %>%
  mutate(REGION=case_when(str_detect(EXON_GROUP,"1")~"1",TRUE~"-")) %>%
  mutate(REGION=case_when(str_detect(EXON_GROUP,"2")~"2",TRUE~REGION)) %>%
  mutate(REGION=case_when(str_detect(EXON_GROUP,"3")~"3",TRUE~REGION)) %>%
  mutate(REGION=case_when(str_detect(EXON_GROUP,"4")~"4",TRUE~REGION)) %>%
  mutate(REGION=case_when(str_detect(EXON_GROUP,"5")~"5",TRUE~REGION)) %>%
  mutate(REGION=case_when(str_detect(EXON_GROUP,"6")~"6",TRUE~REGION)) %>%
  mutate(REGION=case_when(str_detect(EXON_GROUP,"7")~"7",TRUE~REGION)) %>%
  mutate(REGION=case_when(str_detect(EXON_GROUP,"8")~"8",TRUE~REGION)) %>%
  mutate(REGION=case_when(str_detect(EXON_GROUP,"9")~"9",TRUE~REGION)) %>%
  mutate(REGION=case_when(str_detect(EXON_GROUP,"10")~"10",TRUE~REGION)) %>%
  mutate(REGION=case_when(str_detect(EXON_GROUP,"11")~"11",TRUE~REGION)) %>%
  mutate(REGION=case_when(str_detect(EXON_GROUP,"12")~"12",TRUE~REGION)) %>%
  mutate(REGION=case_when(str_detect(EXON_GROUP,"13")~"13",TRUE~REGION)) %>%
  mutate(REGION=case_when(str_detect(EXON_GROUP,"14")~"14",TRUE~REGION)) %>%
  mutate(REGION=case_when(str_detect(EXON_GROUP,"15")~"15",TRUE~REGION)) %>%
  mutate(REGION=case_when(str_detect(EXON_GROUP,"16")~"16",TRUE~REGION)) %>%
  mutate(REGION=case_when(str_detect(EXON_GROUP,"17")~"17",TRUE~REGION))



targeton_collapse <- targeton_input %>%
  group_by(id,REGION) %>%
  dplyr::summarize(variant_n_instances_between_targetons= n(), 
                   #exon_specific_names = paste(combined_name, collapse = '|'),
                   source_dataset = paste(dataset_process, collapse ="|"),
                   source_targetons = paste(EXON_GROUP, collapse = '|'),
                   source_libraries = paste(Variant_Sources, collapse = '|'),
                   variant_n_instances_between_libraries = paste(Variant_duplication, collapse = '|'),
                   PAM_status_simplified_collated = paste(PAM_status_simplified_collated, collapse = '|'), 
                   PAM_status = paste(PAM_status, collapse='|'),
                   dup_mseq_collated = paste(dup_mseq_collated, collapse='|'),
                   arrayed_mseq_a = paste(mseq_a, collapse = '|'), 
                   arrayed_mseq_b = paste(mseq_b, collapse = '|'), 
                   
                   lfcSE_D4_D7 = paste(lfcSE_D4_D7, collapse = ','), 
                   lfcSE_D4_D10 = paste(lfcSE_D4_D10, collapse = ','), 
                   lfcSE_D4_D14 = paste(lfcSE_D4_D14, collapse = ','), 
                   lfcSE_D4_D21 = paste(lfcSE_D4_D21, collapse = ','), 
                   lfcSE_continuous = paste(lfcSE_continuous, collapse = ','), 
                   
                   adj_lfc_D4_D7 = paste(adj_lfc_D4_D7, collapse = ','),
                   adj_lfc_D4_D10 = paste(adj_lfc_D4_D10, collapse = ','),
                   adj_lfc_D4_D14 = paste(adj_lfc_D4_D14, collapse = ','),
                   adj_lfc_D4_D21 = paste(adj_lfc_D4_D21, collapse = ','),
                   adj_lfc_continuous = paste(adj_lfc_continuous, collapse = ',')) %>%
  
  ungroup() %>%
  #ungroup the key data into primary and secondary targetons fields
  separate(adj_lfc_D4_D7, sep=",", c("primary_tg_adj_lfc_D4_D7","secondary_tg_adj_lfc_D4_D7")) %>%
  separate(adj_lfc_D4_D10, sep=",", c("primary_tg_adj_lfc_D4_D10","secondary_tg_adj_lfc_D4_D10")) %>%
  separate(adj_lfc_D4_D14, sep=",", c("primary_tg_adj_lfc_D4_D14","secondary_tg_adj_lfc_D4_D14")) %>%
  separate(adj_lfc_D4_D21, sep=",", c("primary_tg_adj_lfc_D4_D21","secondary_tg_adj_lfc_D4_D21")) %>%
  separate(adj_lfc_continuous, sep=",", c("primary_tg_adj_lfc_continuous","secondary_tg_adj_lfc_continuous")) %>%
  
  separate(lfcSE_D4_D7, sep=",", c("primary_tg_lfcSE_D4_D7","secondary_tg_lfcSE_D4_D7")) %>%
  separate(lfcSE_D4_D10, sep=",", c("primary_tg_lfcSE_D4_D10","secondary_tg_lfcSE_D4_D10"))%>%
  separate(lfcSE_D4_D14, sep=",", c("primary_tg_lfcSE_D4_D14","secondary_tg_lfcSE_D4_D14"))%>%
  separate(lfcSE_D4_D21, sep=",", c("primary_tg_lfcSE_D4_D21","secondary_tg_lfcSE_D4_D21")) %>%
  separate(lfcSE_continuous, sep=",", c("primary_tg_lfcSE_continuous","secondary_tg_lfcSE_continuous"))

#make the newly created fields numeric in order to do weighted calculations
targeton_collapse[, c(13:32)] <- sapply(targeton_collapse[, c(13:32)], as.numeric)
#make the number of targeton source field character based to get subsequent string detection to work
targeton_collapse$variant_n_instances_between_targetons<-as.character(targeton_collapse$variant_n_instances_between_targetons)

#PERFORM THE WEIGHTED CALCULATIONS
targeton_processed <-targeton_collapse %>%
  mutate(weight_tg1_D4_D7=case_when(str_detect(variant_n_instances_between_targetons,"1")~as.numeric(NA), TRUE~ 1/(primary_tg_lfcSE_D4_D7)^2)) %>%
  mutate(weight_tg1_D4_D10=case_when(str_detect(variant_n_instances_between_targetons,"1")~as.numeric(NA), TRUE~ 1/(primary_tg_lfcSE_D4_D10)^2)) %>%
  mutate(weight_tg1_D4_D14=case_when(str_detect(variant_n_instances_between_targetons,"1")~as.numeric(NA), TRUE~ 1/(primary_tg_lfcSE_D4_D14)^2)) %>%
  mutate(weight_tg1_D4_D21=case_when(str_detect(variant_n_instances_between_targetons,"1")~as.numeric(NA), TRUE~ 1/(primary_tg_lfcSE_D4_D21)^2)) %>%
  mutate(weight_tg1_continuous=case_when(str_detect(variant_n_instances_between_targetons,"1")~as.numeric(NA), TRUE~ 1/(primary_tg_lfcSE_continuous)^2)) %>%
  
  mutate(weight_tg2_D4_D7=case_when(str_detect(variant_n_instances_between_targetons,"1")~as.numeric(NA), TRUE~ 1/(secondary_tg_lfcSE_D4_D7)^2)) %>%
  mutate(weight_tg2_D4_D10=case_when(str_detect(variant_n_instances_between_targetons,"1")~as.numeric(NA), TRUE~ 1/(secondary_tg_lfcSE_D4_D10)^2)) %>%
  mutate(weight_tg2_D4_D14=case_when(str_detect(variant_n_instances_between_targetons,"1")~as.numeric(NA), TRUE~ 1/(secondary_tg_lfcSE_D4_D14)^2)) %>%
  mutate(weight_tg2_D4_D21=case_when(str_detect(variant_n_instances_between_targetons,"1")~as.numeric(NA), TRUE~ 1/(secondary_tg_lfcSE_D4_D21)^2)) %>%
  mutate(weight_tg2_continuous=case_when(str_detect(variant_n_instances_between_targetons,"1")~as.numeric(NA), TRUE~ 1/(secondary_tg_lfcSE_continuous)^2)) %>%
  
  mutate(sum_of_weight_D4_D7=rowSums(cbind(weight_tg1_D4_D7,weight_tg2_D4_D7), na.rm=FALSE))%>%
  mutate(sum_of_weight_D4_D10=rowSums(cbind(weight_tg1_D4_D10,weight_tg2_D4_D10), na.rm=FALSE))%>%
  mutate(sum_of_weight_D4_D14=rowSums(cbind(weight_tg1_D4_D14,weight_tg2_D4_D14), na.rm=FALSE))%>%
  mutate(sum_of_weight_D4_D21=rowSums(cbind(weight_tg1_D4_D21,weight_tg2_D4_D21), na.rm=FALSE))%>%
  mutate(sum_of_weight_continuous=rowSums(cbind(weight_tg1_continuous,weight_tg2_continuous), na.rm=FALSE))%>%
  
  mutate(SE_bind_D4_D7 = case_when(str_detect(variant_n_instances_between_targetons,"1")~primary_tg_lfcSE_D4_D7, TRUE~ (sum_of_weight_D4_D7)^(-0.5))) %>%
  mutate(SE_bind_D4_D10 = case_when(str_detect(variant_n_instances_between_targetons,"1")~primary_tg_lfcSE_D4_D10, TRUE~ (sum_of_weight_D4_D10)^(-0.5))) %>%
  mutate(SE_bind_D4_D14 = case_when(str_detect(variant_n_instances_between_targetons,"1")~primary_tg_lfcSE_D4_D14, TRUE~ (sum_of_weight_D4_D14)^(-0.5))) %>%
  mutate(SE_bind_D4_D21 = case_when(str_detect(variant_n_instances_between_targetons,"1")~primary_tg_lfcSE_D4_D21, TRUE~ (sum_of_weight_D4_D21)^(-0.5))) %>%
  mutate(SE_bind_continuous = case_when(str_detect(variant_n_instances_between_targetons,"1")~primary_tg_lfcSE_continuous, TRUE~ (sum_of_weight_continuous)^(-0.5))) %>%
  
  mutate(weighted_primary_tg_LFC_D4_D7= weight_tg1_D4_D7*primary_tg_adj_lfc_D4_D7) %>% mutate(weighted_secondary_tg_LFC_D4_D7= weight_tg2_D4_D7*secondary_tg_adj_lfc_D4_D7) %>%
  mutate(weighted_primary_tg_LFC_D4_D10= weight_tg1_D4_D10*primary_tg_adj_lfc_D4_D10) %>% mutate(weighted_secondary_tg_LFC_D4_D10= weight_tg2_D4_D10*secondary_tg_adj_lfc_D4_D10) %>%
  mutate(weighted_primary_tg_LFC_D4_D14= weight_tg1_D4_D14*primary_tg_adj_lfc_D4_D14) %>% mutate(weighted_secondary_tg_LFC_D4_D14= weight_tg2_D4_D14*secondary_tg_adj_lfc_D4_D14) %>%
  mutate(weighted_primary_tg_LFC_D4_D21= weight_tg1_D4_D21*primary_tg_adj_lfc_D4_D21) %>% mutate(weighted_secondary_tg_LFC_D4_D21= weight_tg2_D4_D21*secondary_tg_adj_lfc_D4_D21) %>%
  mutate(weighted_primary_tg_LFC_continuous= weight_tg1_continuous*primary_tg_adj_lfc_continuous) %>% mutate(weighted_secondary_tg_LFC_continuous= weight_tg2_continuous*secondary_tg_adj_lfc_continuous) %>%
  
  mutate(sum_of_weighted_LFC_D4_D7=rowSums(cbind(weighted_primary_tg_LFC_D4_D7,weighted_secondary_tg_LFC_D4_D7), na.rm=FALSE)) %>%
  mutate(processed_LFC_D4_D7=case_when(str_detect(variant_n_instances_between_targetons,"1")~primary_tg_adj_lfc_D4_D7, TRUE~ sum_of_weighted_LFC_D4_D7/sum_of_weight_D4_D7)) %>%
  mutate(sum_of_weighted_LFC_D4_D10=rowSums(cbind(weighted_primary_tg_LFC_D4_D10,weighted_secondary_tg_LFC_D4_D10), na.rm=FALSE)) %>%
  mutate(processed_LFC_D4_D10=case_when(str_detect(variant_n_instances_between_targetons,"1")~primary_tg_adj_lfc_D4_D10, TRUE~ sum_of_weighted_LFC_D4_D10/sum_of_weight_D4_D10)) %>%
  mutate(sum_of_weighted_LFC_D4_D14=rowSums(cbind(weighted_primary_tg_LFC_D4_D14,weighted_secondary_tg_LFC_D4_D14), na.rm=FALSE)) %>%
  mutate(processed_LFC_D4_D14=case_when(str_detect(variant_n_instances_between_targetons,"1")~primary_tg_adj_lfc_D4_D14, TRUE~ sum_of_weighted_LFC_D4_D14/sum_of_weight_D4_D14)) %>%
  mutate(sum_of_weighted_LFC_D4_D21=rowSums(cbind(weighted_primary_tg_LFC_D4_D21,weighted_secondary_tg_LFC_D4_D21), na.rm=FALSE)) %>%
  mutate(processed_LFC_D4_D21=case_when(str_detect(variant_n_instances_between_targetons,"1")~primary_tg_adj_lfc_D4_D21, TRUE~ sum_of_weighted_LFC_D4_D21/sum_of_weight_D4_D21)) %>%
  mutate(sum_of_weighted_LFC_continuous=rowSums(cbind(weighted_primary_tg_LFC_continuous,weighted_secondary_tg_LFC_continuous), na.rm=FALSE)) %>%
  mutate(processed_LFC_continuous=case_when(str_detect(variant_n_instances_between_targetons,"1")~primary_tg_adj_lfc_continuous, TRUE~ sum_of_weighted_LFC_continuous/sum_of_weight_continuous)) #%>%

#check that no duplicate ids remain
duplicated <- targeton_processed  %>% 
  group_by(id) %>%
  filter(n() > 1)

#now make a new mseq identifier for the appropriate mseq if either mseq is ok then can use a or b annotation

#we now need to derive the correct annotation source for each variant
#we can use the mseq as an identifier for the variant that can be made unique so that BH FDR calculation is correct for the number of unique observations
#we also need to annotate correctly, taking into account whether the variant scores were derived from library a or library b in the case of the PPE codons

#first = is the variant derived from library a or library b alone? if yes then the correct annotation source for each variant is simply the library used:
targeton_processed_req_anno<-targeton_processed %>% mutate(required_annotation_source=case_when(source_dataset=="not_combined"~source_libraries,TRUE~"-")) %>%
  relocate(required_annotation_source,.after=id) %>%
  #second = if the variant was observed in library a and library b did the variant fall in a PPE codon? if yes then the variant calculations were based on the alternative library alone, so get that annotation
  mutate(required_annotation_source=case_when(source_dataset=="combined" & str_detect(PAM_status_simplified_collated,"B")~"A",TRUE~required_annotation_source)) %>%
  mutate(required_annotation_source=case_when(source_dataset=="combined" & str_detect(PAM_status_simplified_collated,"A")~"B",TRUE~required_annotation_source)) %>%
  #third = is the variant observed in both libraries and also in an overlapping region of a tiled targeton and a PPE codon? if yes the annotation will again be from the alternative libary to the PPE containing library 
  mutate(required_annotation_source=case_when(source_dataset=="combined|combined" & str_detect(PAM_status_simplified_collated,"B")~"A",TRUE~required_annotation_source)) %>%
  mutate(required_annotation_source=case_when(source_dataset=="combined|combined" & str_detect(PAM_status_simplified_collated,"A")~"B",TRUE~required_annotation_source)) %>%
  #fourth = is the variant observed in one targeton where library a and b werd combined, and in one targeton where the only one library was used, and does the variant fall in a PPE? if so then take the non-PPE annotaiton source  
  mutate(required_annotation_source=case_when(source_dataset=="combined|not_combined" & str_detect(PAM_status_simplified_collated,"B")~"A",TRUE~required_annotation_source)) %>%
  mutate(required_annotation_source=case_when(source_dataset=="combined|not_combined" & str_detect(PAM_status_simplified_collated,"A")~"B",TRUE~required_annotation_source)) %>%
  #fifth = is the variant not so far allocated an annotation source? therefore it is not a PPE codon. Is the variant only observed in library a or library b, if so then that library is the correct source  
  mutate(required_annotation_source=case_when(required_annotation_source=="-" & source_libraries=="A"~"A",TRUE~required_annotation_source)) %>%
  mutate(required_annotation_source=case_when(required_annotation_source=="-" & source_libraries=="B"~"B",TRUE~required_annotation_source)) %>%
  #sixth = is the variant observed in only one library a or b? and is the variant in an overlapping region between targetons? if yes, then again the annotaiton source is that library (a or b)  
  mutate(required_annotation_source=case_when(required_annotation_source=="-" & source_libraries=="A|A"~"A",TRUE~required_annotation_source)) %>%
  mutate(required_annotation_source=case_when(required_annotation_source=="-" & source_libraries=="B|B"~"B",TRUE~required_annotation_source)) %>%
  #finally = if the variant has so far not ben allocated an annotation source, then either library a or library b mseq can be used as an annotation as they have identical VEP and VaLiAnT annotations  
  mutate(required_annotation_source=case_when(required_annotation_source=="-" ~ "either",TRUE~required_annotation_source))


#seperate out the arrayed mseq field (that were combined into one string at the weighted overlapping tile step, and separated by "|" with targeton 1 on the left, targeton 2 on the right)
targeton_processed_req_anno_2<-targeton_processed_req_anno %>% separate(arrayed_mseq_a, sep='\\|',c("mseq_a_tg1","mseq_a_tg2")) %>%
  separate(arrayed_mseq_b, sep='\\|',c("mseq_b_tg1","mseq_b_tg2")) 


##dealing with specific deletion annotations - START
##dealing with specific deletion annotations - START

#16 duplicates (8 variants duplicated in the final mseq unique dataframe - theses are all 1bp deletions that produce the same genomic edit at a polynucleotide
#1 duplicate is in a pam codon and hence takes the mseq from the alternative library
#1 duplicate is not in the pam codon but produces the same mseq, so can be taken from either and is taken from a
#hence the duplication
#each has the same classificaiton of processed_LFC_continuous and processed_BH_FDR so can be made safely unique without problems - however for accurate annotation the following was done

#the following are specific annotations for BAP1 at polynucleotide sequences adjacent to select PPE codons
#example a deletion of either G at [GG] would give the same genomic edit. If  one G is in a sgRNA PPE codon then mseq_b would be chosen
#however if the deletion of the next G is made, this gives the same genomic edit, but no PPE annotation as not in a PPE codon - but adjacent 
#this results in one deletion being assigned one mseq and the other a different mseq even though they result in the same change
#this occurs at 8 deletions in the library and must be addressed specifically 
#if this next code block is ignored then the deletions can be found at the following step : #check for any duplicates in annotation --- SEE ABOVE*

#LIBRARY MSEQ correction
targeton_processed_req_anno_mseq_deletions<-targeton_processed_req_anno_2 %>% 
  mutate(required_annotation_source=case_when(id=="ENST00000460680.6.ENSG00000163930.10_chr3:52402844_1del_rc" ~ "B",TRUE~required_annotation_source)) %>%
  mutate(required_annotation_source=case_when(id=="ENST00000460680.6.ENSG00000163930.10_chr3:52403611_1del_rc" ~ "B",TRUE~required_annotation_source)) %>%
  mutate(required_annotation_source=case_when(id=="ENST00000460680.6.ENSG00000163930.10_chr3:52405218_52405220_inframe_rc" ~ "B",TRUE~required_annotation_source)) %>%
  mutate(required_annotation_source=case_when(id=="ENST00000460680.6.ENSG00000163930.10_chr3:52405220_1del_rc" ~ "B",TRUE~required_annotation_source)) %>%
  mutate(required_annotation_source=case_when(id=="ENST00000460680.6.ENSG00000163930.10_chr3:52408584_1del_rc" ~ "B",TRUE~required_annotation_source)) %>%
  mutate(required_annotation_source=case_when(id=="ENST00000460680.6.ENSG00000163930.10_chr3:52409854_1del_rc" ~ "B",TRUE~required_annotation_source))

#now create a field with the appropriate mseq for the variant - most variants are not in overlapping region so take targeton 1 (either a or b) for the overlapping take also take targeton 1 as the unqiue identifer
targeton_processed_req_anno_mseq<-targeton_processed_req_anno_mseq_deletions %>% mutate(mseq_for_filter=case_when(required_annotation_source=="A" ~ mseq_a_tg1,TRUE~"-")) %>%
  mutate(mseq_for_filter=case_when(required_annotation_source=="B" ~ mseq_b_tg1,TRUE~mseq_for_filter)) %>%
  mutate(mseq_for_filter=case_when(required_annotation_source=="either" ~ mseq_a_tg1,TRUE~mseq_for_filter)) %>%
  relocate(mseq_for_filter,.after=id)

#TARGETON library correction 
#these are two deletion variants that need annotation correction for a different reason than above
#they require annotation from the mseq_a in targeton 17_2 and 13_1, respectively
#currently an adjacent 1bp deletion from 17_1 and 13_2, which gives the same sequence is being annotated as well - but pulling a different mseq for the annotation
#chosing a single targeton source for the duplicate annotations solves the issue

targeton_processed_req_anno_mseq<-targeton_processed_req_anno_mseq %>%
  mutate(mseq_for_filter=case_when(id=="ENST00000460680.6.ENSG00000163930.10_chr3:52402303_1del_rc" ~ mseq_a_tg2,TRUE~mseq_for_filter)) %>%
  mutate(mseq_for_filter=case_when(id=="ENST00000460680.6.ENSG00000163930.10_chr3:52403779_1del_rc" ~ mseq_a_tg2,TRUE~mseq_for_filter))

##dealing with specific deletion annotations - END
##dealing with specific deletion annotations - END
##dealing with specific deletion annotations - END

#make unique for the mseq
targeton_processed_req_anno_mseq_unique<-targeton_processed_req_anno_mseq[!duplicated(targeton_processed_req_anno_mseq$mseq_for_filter), ]

#check that no duplicate mseqs remain
duplicated_mseq <- targeton_processed_req_anno_mseq_unique  %>% 
  group_by(mseq_for_filter) %>%
  filter(n() > 1)

#produce the final statistical model now that the variant rows are truly unique 
final_stats_frame<-targeton_processed_req_anno_mseq_unique %>%
  #PERFORM THE CALCULATION OF Z SCORES > P VALUES > BH_FDR - this is the final statistics that will be used 
  mutate(processed_Z_D4_D7=processed_LFC_D4_D7/SE_bind_D4_D7) %>% mutate(two_tailed_p_D4_D7= pnorm(abs(processed_Z_D4_D7),lower.tail = FALSE) *2) %>% mutate(processed_BH_FDR_D4_D7 = p.adjust(two_tailed_p_D4_D7, method = "BH")) %>%
  mutate(processed_Z_D4_D10=processed_LFC_D4_D10/SE_bind_D4_D10) %>% mutate(two_tailed_p_D4_D10= pnorm(abs(processed_Z_D4_D10),lower.tail = FALSE) *2) %>% mutate(processed_BH_FDR_D4_D10 = p.adjust(two_tailed_p_D4_D10, method = "BH")) %>%
  mutate(processed_Z_D4_D14=processed_LFC_D4_D14/SE_bind_D4_D14) %>% mutate(two_tailed_p_D4_D14= pnorm(abs(processed_Z_D4_D14),lower.tail = FALSE) *2) %>% mutate(processed_BH_FDR_D4_D14 = p.adjust(two_tailed_p_D4_D14, method = "BH")) %>%
  mutate(processed_Z_D4_D21=processed_LFC_D4_D21/SE_bind_D4_D21) %>% mutate(two_tailed_p_D4_D21= pnorm(abs(processed_Z_D4_D21),lower.tail = FALSE) *2) %>% mutate(processed_BH_FDR_D4_D21 = p.adjust(two_tailed_p_D4_D21, method = "BH")) %>%
  mutate(processed_Z_continuous=processed_LFC_continuous/SE_bind_continuous) %>% mutate(two_tailed_p_continuous= pnorm(abs(processed_Z_continuous),lower.tail = FALSE) *2) %>% mutate(processed_BH_FDR_continuous = p.adjust(two_tailed_p_continuous, method = "BH"))

#metrics calculation
final_stats_frame$functional_classification<-NA
final_stats_frame$functional_classification[final_stats_frame$processed_BH_FDR_continuous<0.01 & final_stats_frame$processed_LFC_continuous<0] <- 'depleted'
final_stats_frame$functional_classification[final_stats_frame$processed_BH_FDR_continuous<0.01 & final_stats_frame$processed_LFC_continuous>0] <- 'enriched'
final_stats_frame$functional_classification[final_stats_frame$processed_BH_FDR_continuous>=0.01] <- 'unchanged'

#check the categories
table(final_stats_frame$functional_classification)
#check there are no NA rows in new classification column
is_na_class<-final_stats_frame %>% filter(functional_classification %in% NA)


#single library stats and scores to check for library a and library b concordance where possible 
concordance_frame_a<-sga_sgb_combined2 %>% select(8,30,31,40,41) %>% dplyr::rename(mseq_for_filter= "mseq_a")
concordance_frame_b<-sga_sgb_combined2 %>% select(9,30,31,40,41) %>% dplyr::rename(mseq_for_filter= "mseq_b")
concordance_frame_both<-do.call("rbind.fill", list(concordance_frame_a,concordance_frame_b))
concordance_frame_both<-concordance_frame_both[!duplicated(concordance_frame_both$mseq_for_filter), ]
concordance_frame<-merge(final_stats_frame,concordance_frame_both,by="mseq_for_filter",all.x=TRUE)
concordance_frame<-concordance_frame %>%
  mutate(processed_Z_continuous_A=sga_adj_lfc_continuous/sga_lfcSE_continuous) %>% mutate(two_tailed_p_continuous_A= pnorm(abs(processed_Z_continuous_A),lower.tail = FALSE) *2) %>% mutate(processed_BH_FDR_continuous_A = p.adjust(two_tailed_p_continuous_A, method = "BH")) %>%
  mutate(processed_Z_continuous_B=sgb_adj_lfc_continuous/sgb_lfcSE_continuous) %>% mutate(two_tailed_p_continuous_B= pnorm(abs(processed_Z_continuous_B),lower.tail = FALSE) *2) %>% mutate(processed_BH_FDR_continuous_B = p.adjust(two_tailed_p_continuous_B, method = "BH"))
#metics calculation based on single libraries
library_concordance_calculated<-concordance_frame

library_concordance_calculated<-library_concordance_calculated %>% mutate(functional_classification_a = case_when(is.na(sga_adj_lfc_continuous) ~ "single_screen", TRUE ~ "fill"))
library_concordance_calculated<-library_concordance_calculated %>% mutate(functional_classification_b = case_when(is.na(sgb_adj_lfc_continuous) ~ "single_screen", TRUE ~ "fill"))


library_concordance_calculated$functional_classification_a[library_concordance_calculated$processed_BH_FDR_continuous_A<0.01 & library_concordance_calculated$sga_adj_lfc_continuous<0] <- 'depleted'
library_concordance_calculated$functional_classification_a[library_concordance_calculated$processed_BH_FDR_continuous_A<0.01 & library_concordance_calculated$sga_adj_lfc_continuous>0] <- 'enriched'
library_concordance_calculated$functional_classification_a[library_concordance_calculated$processed_BH_FDR_continuous_A>=0.01] <- 'unchanged'

library_concordance_calculated$functional_classification_b[library_concordance_calculated$processed_BH_FDR_continuous_B<0.01 & library_concordance_calculated$sgb_adj_lfc_continuous<0] <- 'depleted'
library_concordance_calculated$functional_classification_b[library_concordance_calculated$processed_BH_FDR_continuous_B<0.01 & library_concordance_calculated$sgb_adj_lfc_continuous>0] <- 'enriched'
library_concordance_calculated$functional_classification_b[library_concordance_calculated$processed_BH_FDR_continuous_B>=0.01] <- 'unchanged'


library_concordance_calculated$concordance_type[#library_concordance_calculated$sga_adj_lfc_continuous<0 & 
  #library_concordance_calculated$sgb_adj_lfc_continuous<0 &
  library_concordance_calculated$functional_classification_a== "depleted" &
    library_concordance_calculated$functional_classification_b == "depleted"] <- "concordant_depleted"


library_concordance_calculated$concordance_type[#library_concordance_calculated$sga_adj_lfc_continuous<0 & 
  #library_concordance_calculated$sgb_adj_lfc_continuous<0 &
  library_concordance_calculated$functional_classification_a== "unchanged" &
    library_concordance_calculated$functional_classification_b == "unchanged"] <- "concordant_unchanged"

library_concordance_calculated$concordance_type[#library_concordance_calculated$sga_adj_lfc_continuous<0 & 
  #library_concordance_calculated$sgb_adj_lfc_continuous<0 &
  library_concordance_calculated$functional_classification_a== "enriched" &
    library_concordance_calculated$functional_classification_b == "enriched"] <- "concordant_enriched"

library_concordance_calculated$concordance_type[#library_concordance_calculated$sga_adj_lfc_continuous<0 & 
  #library_concordance_calculated$sgb_adj_lfc_continuous<0 &
  library_concordance_calculated$functional_classification_a!=library_concordance_calculated$functional_classification_b] <- "discordant_classification"


library_concordance_calculated$concordance_type[#library_concordance_calculated$sga_adj_lfc_continuous<0 & 
  #library_concordance_calculated$sgb_adj_lfc_continuous<0 &
  library_concordance_calculated$functional_classification_a== "single_screen" |
    library_concordance_calculated$functional_classification_b == "single_screen"] <- "single_observation"


#WRITE OUT this dataframe - there are many fields that are not necessary to keep but good for checking
write.csv(library_concordance_calculated,"./sga_sgb_combined_and_targeton_processed.csv", row.names = FALSE)

#CLEANUP THE DATAFRAME AND ADD FUNCTIONAL METRICS
#CLEANUP THE DATAFRAME AND ADD FUNCTIONAL METRICS
#CLEANUP THE DATAFRAME AND ADD FUNCTIONAL METRICS
#expand back to full annotation
clean_one<-library_concordance_calculated %>% dplyr::rename(mseq= "mseq_for_filter")
clean_one<-merge(clean_one,total_META_VEP,by='mseq',all.x=FALSE)

#check for any duplicates in annotation
#check for any duplicates in annotation
#check for any duplicates in annotation --- SEE ABOVE*
duplicated <- clean_one  %>% 
  group_by(id.y) %>%
  filter(n() > 1)
write.csv(duplicated,"./duplicated.csv", row.names = FALSE)
#check for any duplicates in annotation
#check for any duplicates in annotation
#check for any duplicates in annotation --- SEE ABOVE*

#simplify the pam status field and create a flag field for variants that only come from one dataset and are pam variant codon
clean_one<- clean_one %>% 
  mutate(pam_codon=case_when(str_detect(PAM_status_simplified_collated,"A") ~ "Y", TRUE ~ "N")) %>%
  mutate(pam_codon=case_when(str_detect(PAM_status_simplified_collated,"B") ~ "Y", TRUE ~ pam_codon)) %>%
  mutate(pam_flag=case_when(source_dataset=="not_combined" & str_detect(PAM_status_simplified_collated,"A") ~ "Y", TRUE ~ "N")) %>%
  mutate(pam_flag=case_when(source_dataset=="not_combined" & str_detect(PAM_status_simplified_collated,"B") ~ "Y", TRUE ~ pam_flag))

#reorder the dataframe by mutator factor
x_levels<-c("snv","1del","custom","2del0","2del1","inframe","stop","ala","snvre")

#re-order by mutator according to set levels
clean_one<-clean_one %>% mutate(mutator =  factor(mutator, levels = x_levels)) %>%
  arrange(mutator)

write.csv(clean_one,"./clean_one.csv", row.names = FALSE)

#select and reorder key columns from  the dataframe
clean_final_annotated_dataset<-clean_one[, c(
  "mseq",
  "id.y",
  "HGVSc",
  "HGVSp",
  "functional_classification",
  "EXON",
  "INTRON",
  "Consequence",
  "REGION",
  "mut_position",
  "ref",
  "new",
  "source_dataset",
  "source_targetons",
  "pam_codon",
  "pam_flag",
  "ref_chr",
  "ref_start",
  "ref_end",
  "ref_aa",
  "alt_aa",
  "cDNA_position",
  "CDS_position",
  "Protein_position",
  "Amino_acids",
  "Codons",
  "mut_type",
  "mutator",
  "processed_LFC_D4_D7",
  "processed_LFC_D4_D10",
  "processed_LFC_D4_D14",
  "processed_LFC_D4_D21",
  "processed_LFC_continuous",
  "SE_bind_D4_D7",
  "SE_bind_D4_D10",
  "SE_bind_D4_D14",
  "SE_bind_D4_D21",
  "SE_bind_continuous",
  "processed_Z_D4_D7",
  "processed_Z_D4_D10",
  "processed_Z_D4_D14",
  "processed_Z_D4_D21",
  "processed_Z_continuous",
  "processed_BH_FDR_D4_D7",
  "processed_BH_FDR_D4_D10",
  "processed_BH_FDR_D4_D14",
  "processed_BH_FDR_D4_D21",
  "processed_BH_FDR_continuous",
  "sga_adj_lfc_continuous",
  "sgb_adj_lfc_continuous",
  "sga_lfcSE_continuous",
  "sgb_lfcSE_continuous",
  "processed_Z_continuous_A",
  "processed_Z_continuous_B",
  "processed_BH_FDR_continuous_A",
  "processed_BH_FDR_continuous_B",
  "concordance",
  "SIFT",
  "PolyPhen",
  "vcf_alias",
  "vcf_var_id",
  "vcf_var_in_const",
  "mave_nt",
  "DOMAINS",
  "Existing_variation",
  "HGVS_OFFSET",
  "AF",
  "gnomAD_AF",
  "gnomAD_AFR_AF",
  "gnomAD_AMR_AF",
  "gnomAD_ASJ_AF",
  "gnomAD_EAS_AF",
  "gnomAD_FIN_AF",
  "gnomAD_NFE_AF",
  "gnomAD_OTH_AF",
  "gnomAD_SAS_AF",
  "gnomAD",
  "gnomAD_FLAG",
  "gnomAD_AF.1")]

#rename columns

clean_final_annotated_dataset_reordered <- clean_final_annotated_dataset %>% dplyr::rename(id = "id.y") %>%
  dplyr::rename(source_target_regions= "source_targetons") %>%
  dplyr::rename(LFC_continuous_A= "sga_adj_lfc_continuous") %>%
  dplyr::rename(LFC_continuous_B= "sgb_adj_lfc_continuous") %>%
  dplyr::rename(SE_continuous_A= "sga_lfcSE_continuous") %>%
  dplyr::rename(SE_continuous_B= "sgb_lfcSE_continuous")

#Consequence simplification
clean_bap1_v1<- clean_final_annotated_dataset_reordered  %>% mutate(slim_consequence=case_when(str_detect(Consequence,"stop_gained") ~ "stop_gained")) %>%
  mutate(slim_consequence=case_when(str_detect(Consequence,"UTR_variant") ~ "UTR", TRUE ~ slim_consequence)) %>%
  mutate(slim_consequence=case_when(str_detect(Consequence,"frameshift_variant") ~ "frameshift", TRUE ~ slim_consequence)) %>% 
  mutate(slim_consequence=case_when(str_detect(mutator,"inframe") ~ "codon_deletion", TRUE ~ slim_consequence)) %>% 
  mutate(slim_consequence=case_when(str_detect(Consequence,"intron_variant") ~ "intron", TRUE ~ slim_consequence)) %>% 
  mutate(slim_consequence=case_when(str_detect(Consequence,"missense_variant") ~ "missense", TRUE ~ slim_consequence)) %>% 
  mutate(slim_consequence=case_when(str_detect(Consequence,"synonymous_variant") ~ "synonymous", TRUE ~ slim_consequence)) %>% 
  mutate(slim_consequence=case_when(str_detect(Consequence,"splice_acceptor_variant") ~ "splice_acceptor", TRUE ~ slim_consequence)) %>% 
  mutate(slim_consequence=case_when(str_detect(Consequence,"splice_donor_variant") ~ "splice_donor", TRUE ~ slim_consequence)) %>%
  
  mutate(slim_consequence=case_when(str_detect(Consequence,"start_lost") & mutator== "snvre" |
                                      str_detect(Consequence,"start_lost") & mutator== "snv" |
                                      str_detect(Consequence,"start_lost") & mutator== "custom" |
                                      str_detect(Consequence,"start_lost") & mutator== "ala"~ "start_lost", TRUE ~ slim_consequence)) %>% 
  
  mutate(slim_consequence=case_when(str_detect(Consequence,"stop_lost") & mutator== "snvre" |
                                      str_detect(Consequence,"stop_lost") & mutator== "snv" |
                                      str_detect(Consequence,"stop_lost") & mutator== "custom" |
                                      str_detect(Consequence,"stop_lost") & mutator== "ala"~ "stop_lost", TRUE ~ slim_consequence)) %>% 
  
  
  mutate(slim_consequence=case_when(str_detect(Consequence,"stop_retained_variant") & mutator=="snvre" |
                                      str_detect(Consequence,"stop_retained_variant") & mutator=="snv" |
                                      str_detect(Consequence,"stop_retained_variant") & mutator=="custom" ~ "synonymous", TRUE ~ slim_consequence)) %>% 
  
  mutate(slim_consequence=case_when(str_detect(Consequence,"inframe_deletion") & mutator=="custom" ~ "clinical_inframe_deletion", TRUE ~ slim_consequence)) %>% 
  mutate(slim_consequence=case_when(str_detect(Consequence,"inframe_insertion") & mutator=="custom" ~ "clinical_inframe_insertion", TRUE ~ slim_consequence)) 
#relocate new column
clean_bap1_v1 <- clean_bap1_v1 %>% relocate(slim_consequence, .after = Consequence)

#add LFC_continuous as a new column called functional socre
clean_bap1_v1$functional_score <- clean_bap1_v1$processed_LFC_continuous
#relocate new column
clean_bap1_v1 <- clean_bap1_v1 %>% relocate(functional_score, .after = functional_classification)


#reorder the dataframe by mutator factor
x_levels<-c("snv","1del","custom","2del0","2del1","inframe","stop","ala","snvre")

#re-order by mutator according to set levels
clean_bap1_v1<-clean_bap1_v1 %>% mutate(mutator =  factor(mutator, levels = x_levels)) %>%
  arrange(mutator)

write.csv(clean_bap1_v1,"./clean_bap1_v1.csv", row.names = FALSE)

#bind EVE scores - 'EVE' is a dataframe from VEP which includes EVE score and classification for BAP1 (ENST00000460680)
evo_EVE<-EVE %>% filter (Feature  %in% "ENST00000460680") %>% select('X.Uploaded_variation',"EVE_CLASS","EVE_SCORE") %>% dplyr::rename(HGVSc = "X.Uploaded_variation")
evo_EVE <- evo_EVE[!duplicated(evo_EVE$HGVSc), ]
eve_EVE_1<-merge(clean_bap1_v1,evo_EVE,by="HGVSc",all.x=TRUE)
eve_EVE_1$EVE_SCORE<-as.numeric(eve_EVE_1$EVE_SCORE)
#rename
clean_bap1_v2<-eve_EVE_1

#bind clinvar release 1 star only - 'one_star_040923' is a dataframe downloaded on 4th September 2023 which contains >=1* review status variants
one_star_040923$CLIN_SIG<-gsub("\\s*\\([^\\)]+\\)","",as.character(one_star_040923$Clinical.significance..Last.reviewed.))
one_star_040923$HGVSc<-gsub("\\s*\\([^\\)]+\\)","",as.character(one_star_040923$Name))
one_star_040923$HGVSc<-str_replace(one_star_040923$HGVSc, "NM_004656.4", "ENST00000460680.6")
sept_23_clinvar<-one_star_040923 %>% select(c(12,14,17,18))
clean_bap1_v2<-merge(clean_bap1_v2,sept_23_clinvar, by="HGVSc",all.x=TRUE)

#test<-merge(clean_bap1_v2,sept_23_clinvar, by="HGVSc",all.x=FALSE)
#test <- test[!duplicated(test$mseq), ]
#test <- test %>% filter(!pam_flag %in% "Y")

#dataframe for the plotting on database information (clinvar and gnomad)
databases<-clean_bap1_v2 %>% mutate(is_in_clinvar=case_when(is.na(VariationID) ~ "N", TRUE ~ "Y")) %>%
  mutate(is_in_gnomad=case_when(gnomAD == "-" ~ "N", TRUE ~ "Y")) %>%
  mutate(present_accession=case_when(is_in_clinvar == "N" ~ "Unobserved", TRUE ~ "ClinVar only")) %>%
  mutate(present_accession=case_when(is_in_clinvar == "N" & is_in_gnomad =="Y" ~ "gnomAD only", TRUE ~present_accession)) %>%
  mutate(present_accession=case_when(is_in_clinvar == "Y" & is_in_gnomad =="Y" ~ "ClinVar and gnomAD", TRUE ~present_accession)) %>%
  mutate(ClinVar_CLNSIG_slim=case_when(CLIN_SIG=="Uncertain significance" ~ "Uncertain significance")) %>%
  mutate(ClinVar_CLNSIG_slim=case_when(is.na(CLIN_SIG) ~ "Unobserved", TRUE ~ClinVar_CLNSIG_slim)) %>%
  mutate(ClinVar_CLNSIG_slim=case_when(CLIN_SIG=="Pathogenic/Likely pathogenic" ~ "Pathogenic/Likely pathogenic", TRUE ~ClinVar_CLNSIG_slim)) %>%
  mutate(ClinVar_CLNSIG_slim=case_when(CLIN_SIG=="Pathogenic" ~ "Pathogenic/Likely pathogenic", TRUE ~ClinVar_CLNSIG_slim)) %>%
  mutate(ClinVar_CLNSIG_slim=case_when(CLIN_SIG=="Likely pathogenic" ~ "Pathogenic/Likely pathogenic", TRUE ~ClinVar_CLNSIG_slim)) %>%
  mutate(ClinVar_CLNSIG_slim=case_when(CLIN_SIG=="Likely benign" ~ "Benign/Likely benign", TRUE ~ClinVar_CLNSIG_slim)) %>%
  mutate(ClinVar_CLNSIG_slim=case_when(CLIN_SIG=="Benign/Likely benign" ~ "Benign/Likely benign", TRUE ~ClinVar_CLNSIG_slim)) %>%
  mutate(ClinVar_CLNSIG_slim=case_when(CLIN_SIG=="Benign" ~ "Benign/Likely benign", TRUE ~ClinVar_CLNSIG_slim)) %>%
  mutate(ClinVar_CLNSIG_slim=case_when(CLIN_SIG=="Conflicting interpretations of pathogenicity" ~ "Conflicting interpretation", TRUE ~ClinVar_CLNSIG_slim)) %>%
  mutate(present_accession_filter=paste0(present_accession,"_",mseq))
#reorder the dataframe by present_accession factor
database_levels<-c("ClinVar and gnomAD","ClinVar only","gnomAD only","Unobserved")
#re-order by present_accession to set levels
databases<-databases %>% mutate(present_accession =  factor(present_accession, levels = database_levels)) %>%
  arrange(present_accession)
#save full dataframe
databases_expanded <- databases
write.csv(databases_expanded,"./databases_expanded.csv", row.names = FALSE)
#make unique for variants in each database
databases <- databases %>% filter(!pam_flag %in% "Y")
databases <- databases[!duplicated(databases$mseq), ]    
databases <- databases[!duplicated(databases$HGVSc), ]
#databases <- databases[!duplicated(databases$present_accession_filter), ]    
#filter for clinvar only for some plots
databases_clinvar <- databases %>% filter(present_accession %in% "ClinVar only" | present_accession %in% "ClinVar and gnomAD")
databases_clinvar <- databases_clinvar[!duplicated(databases_clinvar$mseq), ]    
#separate out gnomad3 and gnomad2 allele frequencies
databases_expanded<-databases_expanded %>% separate(gnomAD_AF, sep=',',c("gnomAD_AF_v3","gnomAD_AF_v2")) 
#gnomad_only_df<-databases_expanded %>% filter(present_accession %in% "gnomAD only")
#gnomad_only_df <- gnomad_only_df[!duplicated(gnomad_only_df$HGVSc), ]    
#gnomad_only_df$gnomAD_AF_v3<-as.numeric(gnomad_only_df$gnomAD_AF_v3)

#create dataframe for the plotting on aminoacid information
clean_bap1_v3<-databases_expanded
#reorder the dataframe by mutator factor
x_levels<-c("snv","1del","custom","2del0","2del1","inframe","stop","ala","snvre")

#re-order by mutator according to set levels
clean_bap1_v3<-clean_bap1_v3 %>% mutate(mutator =  factor(mutator, levels = x_levels)) %>%
  arrange(mutator,-mut_position) %>% filter(!HGVSc %in% "-") #%>% filter(!pam_flag %in% "Y")

write.csv(clean_bap1_v3,"./clean_bap1_v3.csv", row.names = FALSE)

#make the dataframe unique for varaints edited into the genome
clean_bap1_v3_unique <- clean_bap1_v3[!duplicated(clean_bap1_v3$mseq), ]  
clean_bap1_v3_unique <- clean_bap1_v3_unique[!duplicated(clean_bap1_v3_unique$HGVSc), ]  

#clean_bap1_v3 - was fixed for correct nomenclature at chrom_pos_ref_alt and SpliceAI values computed, see SpliceAI code, by Prashant Gupta and returned as 'fixed_261023'.
#clean_bap1_v3 - was fixed for correct nomenclature at chrom_pos_ref_alt and SpliceAI values computed, see SpliceAI code, by Prashant Gupta and returned as 'fixed_261023'.
#clean_bap1_v3 - was fixed for correct nomenclature at chrom_pos_ref_alt and SpliceAI values computed, see SpliceAI code, by Prashant Gupta and returned as 'fixed_261023'.

#final data set clean up 
expanded_dataset<-fixed_261023
#reorder the dataframe by mutator factor
x_levels<-c("snv","1del","custom","2del0","2del1","inframe","stop","ala","snvre")

#re-order by mutator according to set levels
expanded_dataset<-expanded_dataset %>% mutate(mutator =  factor(mutator, levels = x_levels)) %>%
  arrange(mutator,-mut_position) %>% filter(!HGVSc %in% "-")

#merge concordance 

concordance_new<-library_concordance_calculated %>% select(mseq_for_filter,functional_classification_a,functional_classification_b, concordance_type) %>% 
  dplyr::rename(mseq= "mseq_for_filter")

clean_bap1_v4 <- merge(expanded_dataset,concordance_new, by="mseq", all.x=FALSE)

clean_bap1_v4 <- clean_bap1_v4 %>% mutate(mutator =  factor(mutator, levels = x_levels)) %>%
  arrange(mutator,-mut_position) %>% filter(!HGVSc %in% "-")


#addition of concordance column - simplified version
clean_bap1_v4$concordance_a_b<-NA
clean_bap1_v4<-clean_bap1_v4 %>% mutate(concordance_a_b=case_when(str_detect(concordance_type,"concordant_enriched") ~ "concordant_classification", TRUE ~ concordance_a_b)) %>%
  mutate(concordance_a_b=case_when(str_detect(concordance_type,"concordant_depleted") ~ "concordant_classification", TRUE ~ concordance_a_b)) %>%
  mutate(concordance_a_b=case_when(str_detect(concordance_type,"concordant_unchanged") ~ "concordant_classification", TRUE ~ concordance_a_b)) %>%
  mutate(concordance_a_b=case_when(str_detect(concordance_type,"discordant_classification") ~ "discordant_classification", TRUE ~ concordance_a_b)) %>%
  mutate(concordance_a_b=case_when(str_detect(concordance_type,"single_observation") ~ "single_observation", TRUE ~ concordance_a_b))

write.csv(clean_bap1_v4,"clean_bap1_v4.csv", row.names = FALSE)

#rearrange the dataframe

#selected rows and changed order in excel
clean_bap1_v4_filtered<-clean_bap1_v4[, c(
  "HGVSc",
  "HGVSp",
  "ref_chr",
  "fixed_pos",
  "fixed_ref",
  "fixed_alt",
  "functional_classification",
  "functional_score",
  "EXON",
  "INTRON",
  "Consequence",
  "slim_consequence",
  "ref_aa",
  "alt_aa",
  "cDNA_position",
  "CDS_position",
  "Protein_position",
  "Amino_acids",
  "Codons",
  "DOMAINS",
  "EVE_CLASS",
  "EVE_SCORE",
  "SIFT",
  "PolyPhen",
  "mutator",
  "mseq",
  "pam_codon",
  "pam_flag",
  "REGION",
  "source_dataset",
  "source_target_regions",
  "concordance_a_b",
  "processed_LFC_D4_D7",
  "processed_LFC_D4_D10",
  "processed_LFC_D4_D14",
  "processed_LFC_D4_D21",
  "processed_LFC_continuous",
  "SE_bind_D4_D7",
  "SE_bind_D4_D10",
  "SE_bind_D4_D14",
  "SE_bind_D4_D21",
  "SE_bind_continuous",
  "processed_Z_D4_D7",
  "processed_Z_D4_D10",
  "processed_Z_D4_D14",
  "processed_Z_D4_D21",
  "processed_Z_continuous",
  "processed_BH_FDR_D4_D7",
  "processed_BH_FDR_D4_D10",
  "processed_BH_FDR_D4_D14",
  "processed_BH_FDR_D4_D21",
  "processed_BH_FDR_continuous",
  "is_in_gnomad",
  "gnomAD",
  "gnomAD_AF_v3",
  "gnomAD_AF_v2",
  "gnomAD_AFR_AF",
  "gnomAD_AMR_AF",
  "gnomAD_ASJ_AF",
  "gnomAD_EAS_AF",
  "gnomAD_FIN_AF",
  "gnomAD_NFE_AF",
  "gnomAD_OTH_AF",
  "gnomAD_SAS_AF",
  "is_in_clinvar",
  "VariationID",
  "CLIN_SIG",
  "ClinVar_CLNSIG_slim",
  "dbSNP.ID",
  "DS_AG",
  "DS_AL",
  "DS_DG",
  "DS_DL",
  "DP_AG",
  "DP_AL",
  "DP_DG",
  "DP_DL",
  "REF_DS_AG",
  "REF_DS_AL",
  "REF_DS_DG",
  "REF_DS_DL",
  "ALT_DS_AG",
  "ALT_DS_AL",
  "ALT_DS_DG",
  "ALT_DS_DL",
  "mave_nt")]


#rename columns
clean_bap1_v4_filtered$chrom_pos_ref_alt<-paste0("3_",clean_bap1_v4_filtered$fixed_pos,"_",clean_bap1_v4_filtered$fixed_ref,"_",clean_bap1_v4_filtered$fixed_alt)
clean_bap1_v4_filtered <- clean_bap1_v4_filtered %>% relocate(chrom_pos_ref_alt, .before = ref_chr)

clean_bap1_v4_filtered_renamed <-clean_bap1_v4_filtered %>% dplyr::rename(pos= "fixed_pos") %>% 
  dplyr::rename(ref= "fixed_ref") %>% 
  dplyr::rename(alt= "fixed_alt") %>% 
  dplyr::rename(exon= "EXON") %>% 
  dplyr::rename(intron= "INTRON") %>% 
  dplyr::rename(vep_consequence= "Consequence") %>%
  dplyr::rename(vep_consequence_slim= "slim_consequence") %>%
  dplyr::rename(protein_position= "Protein_position") %>%
  dplyr::rename(amino_acids= "Amino_acids") %>%
  dplyr::rename(codons= "Codons") %>% 
  dplyr::rename(domains= "DOMAINS") %>% 
  dplyr::rename(EVE_class= "EVE_CLASS") %>% 
  dplyr::rename(EVE_score= "EVE_SCORE") %>% 
  dplyr::rename(region= "REGION") %>% 
  dplyr::rename(concordance= "concordance_a_b") %>%
  
  dplyr::rename(variation_id= "VariationID") %>% 
  dplyr::rename(clinvar_clinical_significance= "CLIN_SIG") %>% 
  dplyr::rename(clinvar_clinical_significance_slim= "ClinVar_CLNSIG_slim")


clean_bap1_v4_unique <- clean_bap1_v4_filtered_renamed %>% filter(!pam_flag %in% "Y")
clean_bap1_v4_unique <- clean_bap1_v4_unique[!duplicated(clean_bap1_v4_unique$mseq), ]  
clean_bap1_v4_unique <- clean_bap1_v4_unique[!duplicated(clean_bap1_v4_unique$HGVSc), ]  


###output the final sge_bap1_dataset and sge_bap1_expanded_dataset
write.csv(clean_bap1_v4_filtered_renamed,"sge_bap1_expanded_dataset.csv", row.names = FALSE)

write.csv(clean_bap1_v4_unique,"sge_bap1_dataset.csv", row.names = FALSE)

# for writing a data.frame or list of data.frames to an xlsx file
write.xlsx(clean_bap1_v4_unique, 'sge_bap1_dataset.xlsx')

write.xlsx(clean_bap1_v4_filtered_renamed, 'sge_bap1_expanded_dataset.xlsx')

#END
#