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
  dir=paste0("/Users/aw28/Desktop/BAP1_analysis/counts_sge/","W",exon,"_",sg,"/results/pycroquet/")
  dir_plasmid=paste0("/Users/aw28/Desktop/BAP1_analysis/counts_plasmid/","W",exon,"_",sg,"/results/pycroquet/")
  out=paste0(dir, "E", exon, "_SG", sg, "_count_frame.csv")
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
dir_meta=paste0("/Users/aw28/Desktop/BAP1_analysis/meta/")
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
dir_vep=paste0("/Users/aw28/Desktop/BAP1_analysis/vep_output/")
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
annot_a<-read.csv("/Users/aw28/Desktop/BAP1_analysis/experimental_design/annot_a.csv")
annot_a<-as.matrix.data.frame(annot_a)
annot_a_continuous<-read.csv("/Users/aw28/Desktop/BAP1_analysis/experimental_design/annot_a_continuous.csv")
annot_a_continuous$condition<-as.numeric(annot_a_continuous$condition)


#read in experimental designs SGB
#annot_b<-read.csv("/Users/aw28/Desktop/BAP1_analysis/experimental_design/annot_b.csv")
#annot_b_continuous<-read.csv("/Users/aw28/Desktop/BAP1_analysis/experimental_design/annot_b_continuous.csv")
#annot_b<-as.matrix.data.frame(annot_b)
#annot_b_continuous$condition<-as.numeric(annot_b_continuous$condition)


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
  table_wald <- degComps(dds, combs = "condition", contrast = list("condition_D7_vs_D4","condition_D10_vs_D4" , "condition_D14_vs_D4", "condition_D21_vs_D4"), alpha = 0.05, skip = FALSE, type = "apeglm", pairs = FALSE, fdr = "default")
  
  #Summary Table
  #production of summary table - CHANGE RAW TO SHRUNKEN IF USING SYNONYNOUS NORMALIZATION FACTORS
  summary <- purrr::reduce(c(deg(table_wald[[1]], "shrunken") %>% appendDataFrameColumns(suffix="_D4_D7") %>% as.data.frame() %>% rownames_to_column(var="mseq") %>% list(),
                             deg(table_wald[[2]], "shrunken") %>% appendDataFrameColumns(suffix="_D4_D10") %>% as.data.frame() %>% rownames_to_column(var="mseq") %>% list(),
                             deg(table_wald[[3]], "shrunken") %>% appendDataFrameColumns(suffix="_D4_D14") %>% as.data.frame() %>% rownames_to_column(var="mseq") %>% list(),
                             deg(table_wald[[4]], "shrunken") %>% appendDataFrameColumns(suffix="_D4_D21") %>% as.data.frame() %>% rownames_to_column(var="mseq") %>% list()), left_join, by="mseq") %>% select(-c("baseMean_D4_D10","baseMean_D4_D14","baseMean_D4_D21")) %>% dplyr::rename(baseMean = "baseMean_D4_D7") %>% data.frame()
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
#Run Deseq - SHRUNKEN LFC - CHANGE TO table_wald[[1]] if want raw
dds <- DESeq(dds)
table_wald <- degComps(dds, combs = "condition", alpha = 0.05, skip = FALSE, type = "apeglm", pairs = FALSE, fdr = "default")
rate <- as.data.frame(table_wald[[2]])%>%rownames_to_column(var="mseq")%>%appendDataFrameColumns(suffix="_continuous")
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
annot_a<-read.csv("/Users/aw28/Desktop/BAP1_analysis/experimental_design/annot_a_exclude_12.csv")
annot_a_continuous<-read.csv("/Users/aw28/Desktop/BAP1_analysis/experimental_design/annot_a_continuous_exclude_12.csv")
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


################# COMBINATION OF GUIDES ############# START
################# COMBINATION OF GUIDES ############# START
################# COMBINATION OF GUIDES ############# START
#in this step we take only the libraries where there is a possibility of library duplication (ie. the same variant assayed with guide library a and library b)
#we group by the exon group in which the variant is found - so there will be duplicated ids as a whole
#fields of interest for each duplicated variant are put into a new field and separated by comma then key fields for calculation are ungrouped by comma into new fields
################# COMBINATION OF GUIDES ############# START
################# COMBINATION OF GUIDES ############# START
################# COMBINATION OF GUIDES ############# START
sga_sgb_collapse <- total_bap1  %>% 
  group_by(id,EXON_GROUP,combined_name) %>%
  dplyr::summarize(Variant_duplication= n(), 
                   Variant_Sources = paste(SG, collapse = ','),
                   PAM_status = paste(pam_mut_sgrna_id, collapse = ','),
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
sga_sgb_collapse[, c(9:40)] <- sapply(sga_sgb_collapse[, c(9:40)], as.numeric)

#perform calculations to obtained weighted means -if a duplicated variant is in a pam codon, the weight for that library contribution will be 0 and the alternative library will contribute 100% to the combined value
sga_sgb_combined2 <-sga_sgb_collapse %>%
  mutate(weight_a_D4_D7=case_when(str_detect(PAM_status,"_a")~ 0, TRUE~ 1/(sga_lfcSE_D4_D7)^2)) %>% mutate(weight_b_D4_D7= case_when(str_detect(PAM_status,"_b")~ 0, TRUE~1/(sgb_lfcSE_D4_D7)^2)) %>%
  mutate(weight_a_D4_D10=case_when(str_detect(PAM_status,"_a")~ 0, TRUE~ 1/(sga_lfcSE_D4_D10)^2)) %>% mutate(weight_b_D4_D10= case_when(str_detect(PAM_status,"_b")~ 0, TRUE~1/(sgb_lfcSE_D4_D10)^2)) %>%
  mutate(weight_a_D4_D14=case_when(str_detect(PAM_status,"_a")~ 0, TRUE~ 1/(sga_lfcSE_D4_D14)^2)) %>% mutate(weight_b_D4_D14= case_when(str_detect(PAM_status,"_b")~ 0, TRUE~1/(sgb_lfcSE_D4_D14)^2)) %>%
  mutate(weight_a_D4_D21=case_when(str_detect(PAM_status,"_a")~ 0, TRUE~ 1/(sga_lfcSE_D4_D21)^2)) %>% mutate(weight_b_D4_D21= case_when(str_detect(PAM_status,"_b")~ 0, TRUE~1/(sgb_lfcSE_D4_D21)^2)) %>%
  mutate(weight_a_continuous=case_when(str_detect(PAM_status,"_a")~ 0, TRUE~ 1/(sga_lfcSE_continuous)^2)) %>% mutate(weight_b_continuous= case_when(str_detect(PAM_status,"_b")~ 0, TRUE~1/(sgb_lfcSE_continuous)^2)) %>%
  
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


##########QC combination of guides#################### START
#simplify the PAM_status field so that each specific guide is not colored, but each guide per targeton (ie. not 'sgRNA_1_A' but 'A') 
sga_sgb_combined_guide<-sga_sgb_combined2 %>% mutate(PAM_status_simplified=case_when(str_detect(PAM_status,"_a")~ "A", str_detect(PAM_status,"_b")~ "B", str_detect(PAM_status,"")~ "-", str_detect(PAM_status,",")~ "-"))


#opens a new PDF in the working directory and pastes each plot comparison - only shows those targetons that have a posibility of duplication, some guides will only be present in library a or library b - filered by variant duplicatio field
pdf("correlated_adj_LFC_sga_sgb_combined2.pdf")
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D7 != "NA") %>% filter(sgb_adj_lfc_D4_D7!="NA") %>% mutate(sga_adj_lfc_D4_D7 = as.numeric(sga_adj_lfc_D4_D7)) %>% mutate(sgb_adj_lfc_D4_D7 = as.numeric(sgb_adj_lfc_D4_D7)) %>% ggscatter(x="sga_adj_lfc_D4_D7",y="sgb_adj_lfc_D4_D7", color = "PAM_status_simplified", alpha= 0.35, size=2, palette = c("black", "red", "green"),  add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D7 != "NA") %>% filter(sgb_adj_lfc_D4_D7!="NA") %>% mutate(sga_adj_lfc_D4_D7 = as.numeric(sga_adj_lfc_D4_D7)) %>% mutate(sgb_adj_lfc_D4_D7 = as.numeric(sgb_adj_lfc_D4_D7)) %>% ggscatter(x="sga_adj_lfc_D4_D7",y="combined_LFC_D4_D7", color = "PAM_status_simplified", alpha= 0.35, size=2, palette = c("black", "red", "green"),  add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D7 != "NA") %>% filter(sgb_adj_lfc_D4_D7!="NA") %>% mutate(sga_adj_lfc_D4_D7 = as.numeric(sga_adj_lfc_D4_D7)) %>% mutate(sgb_adj_lfc_D4_D7 = as.numeric(sgb_adj_lfc_D4_D7)) %>% ggscatter(x="combined_LFC_D4_D7",y="sgb_adj_lfc_D4_D7", color = "PAM_status_simplified", alpha= 0.35, size=2, palette = c("black", "red", "green"),  add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D10 != "NA") %>% filter(sgb_adj_lfc_D4_D10!="NA") %>% mutate(sga_adj_lfc_D4_D10 = as.numeric(sga_adj_lfc_D4_D10)) %>% mutate(sgb_adj_lfc_D4_D10 = as.numeric(sgb_adj_lfc_D4_D10)) %>% ggscatter(x="sga_adj_lfc_D4_D10",y="sgb_adj_lfc_D4_D10", color = "PAM_status_simplified", alpha= 0.35, size=2, palette = c("black", "red", "green"),  add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D10 != "NA") %>% filter(sgb_adj_lfc_D4_D10!="NA") %>% mutate(sga_adj_lfc_D4_D10 = as.numeric(sga_adj_lfc_D4_D10)) %>% mutate(sgb_adj_lfc_D4_D10 = as.numeric(sgb_adj_lfc_D4_D10)) %>% ggscatter(x="sga_adj_lfc_D4_D10",y="combined_LFC_D4_D10", color = "PAM_status_simplified", alpha= 0.35, size=2, palette = c("black", "red", "green"),  add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D10 != "NA") %>% filter(sgb_adj_lfc_D4_D10!="NA") %>% mutate(sga_adj_lfc_D4_D10 = as.numeric(sga_adj_lfc_D4_D10)) %>% mutate(sgb_adj_lfc_D4_D10 = as.numeric(sgb_adj_lfc_D4_D10)) %>% ggscatter(x="combined_LFC_D4_D10",y="sgb_adj_lfc_D4_D10", color = "PAM_status_simplified", alpha= 0.35, size=2, palette = c("black", "red", "green"),  add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D14 != "NA") %>% filter(sgb_adj_lfc_D4_D14!="NA") %>% mutate(sga_adj_lfc_D4_D14 = as.numeric(sga_adj_lfc_D4_D14)) %>% mutate(sgb_adj_lfc_D4_D14 = as.numeric(sgb_adj_lfc_D4_D14)) %>% ggscatter(x="sga_adj_lfc_D4_D14",y="sgb_adj_lfc_D4_D14", color = "PAM_status_simplified", alpha= 0.35, size=2, palette = c("black", "red", "green"),  add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D14 != "NA") %>% filter(sgb_adj_lfc_D4_D14!="NA") %>% mutate(sga_adj_lfc_D4_D14 = as.numeric(sga_adj_lfc_D4_D14)) %>% mutate(sgb_adj_lfc_D4_D14 = as.numeric(sgb_adj_lfc_D4_D14)) %>% ggscatter(x="sga_adj_lfc_D4_D14",y="combined_LFC_D4_D14", color = "PAM_status_simplified", alpha= 0.35, size=2, palette = c("black", "red", "green"),  add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D14 != "NA") %>% filter(sgb_adj_lfc_D4_D14!="NA") %>% mutate(sga_adj_lfc_D4_D14 = as.numeric(sga_adj_lfc_D4_D14)) %>% mutate(sgb_adj_lfc_D4_D14 = as.numeric(sgb_adj_lfc_D4_D14)) %>% ggscatter(x="combined_LFC_D4_D14",y="sgb_adj_lfc_D4_D14", color = "PAM_status_simplified", alpha= 0.35, size=2, palette = c("black", "red", "green"),  add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D21 != "NA") %>% filter(sgb_adj_lfc_D4_D21!="NA") %>% mutate(sga_adj_lfc_D4_D21 = as.numeric(sga_adj_lfc_D4_D21)) %>% mutate(sgb_adj_lfc_D4_D21 = as.numeric(sgb_adj_lfc_D4_D21)) %>% ggscatter(x="sga_adj_lfc_D4_D21",y="sgb_adj_lfc_D4_D21", color = "PAM_status_simplified", alpha= 0.35, size=2, palette = c("black", "red", "green"),  add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D21 != "NA") %>% filter(sgb_adj_lfc_D4_D21!="NA") %>% mutate(sga_adj_lfc_D4_D21 = as.numeric(sga_adj_lfc_D4_D21)) %>% mutate(sgb_adj_lfc_D4_D21 = as.numeric(sgb_adj_lfc_D4_D21)) %>% ggscatter(x="sga_adj_lfc_D4_D21",y="combined_LFC_D4_D21", color = "PAM_status_simplified", alpha= 0.35, size=2, palette = c("black", "red", "green"),  add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D21 != "NA") %>% filter(sgb_adj_lfc_D4_D21!="NA") %>% mutate(sga_adj_lfc_D4_D21 = as.numeric(sga_adj_lfc_D4_D21)) %>% mutate(sgb_adj_lfc_D4_D21 = as.numeric(sgb_adj_lfc_D4_D21)) %>% ggscatter(x="combined_LFC_D4_D21",y="sgb_adj_lfc_D4_D21", color = "PAM_status_simplified", alpha= 0.35, size=2, palette = c("black", "red", "green"),  add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_continuous != "NA") %>% filter(combined_LFC_continuous!="NA") %>% mutate(sga_adj_lfc_continuous = as.numeric(sga_adj_lfc_continuous)) %>% mutate(combined_LFC_continuous = as.numeric(combined_LFC_continuous)) %>% ggscatter(x="combined_LFC_continuous",y="sga_adj_lfc_continuous", color = "PAM_status_simplified", alpha= 0.35, size=2, palette = c("black", "red", "green"),  add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sgb_adj_lfc_continuous != "NA") %>% filter(combined_LFC_continuous!="NA") %>% mutate(sgb_adj_lfc_continuous = as.numeric(sgb_adj_lfc_continuous)) %>% mutate(combined_LFC_continuous = as.numeric(combined_LFC_continuous)) %>% ggscatter(x="sgb_adj_lfc_continuous",y="combined_LFC_continuous", color = "PAM_status_simplified", alpha= 0.35, size=2, palette = c("black", "red", "green"),  add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
dev.off()

#AS ABOVE WITH Z SCORE RATHER THAN LFC pdf("correlated_adj_LFC_sga_sgb_combined2.pdf")
pdf("correlated_adj_Z_sga_sgb_combined2.pdf")

sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D7 != "NA") %>% filter(sgb_adj_lfc_D4_D7!="NA") %>% mutate(sga_adj_lfc_D4_D7 = as.numeric(sga_adj_lfc_D4_D7)) %>% mutate(sgb_adj_lfc_D4_D7 = as.numeric(sgb_adj_lfc_D4_D7)) %>% ggscatter(x="sga_adj_lfc_D4_D7",y="sgb_adj_lfc_D4_D7", color = "PAM_status_simplified", alpha= 0.35, size=2, palette = c("black", "red", "green"),  add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D7 != "NA") %>% filter(sgb_adj_lfc_D4_D7!="NA") %>% mutate(sga_adj_lfc_D4_D7 = as.numeric(sga_adj_lfc_D4_D7)) %>% mutate(sgb_adj_lfc_D4_D7 = as.numeric(sgb_adj_lfc_D4_D7)) %>% ggscatter(x="sga_adj_lfc_D4_D7",y="combined_LFC_D4_D7", color = "PAM_status_simplified", alpha= 0.35, size=2, palette = c("black", "red", "green"),  add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D7 != "NA") %>% filter(sgb_adj_lfc_D4_D7!="NA") %>% mutate(sga_adj_lfc_D4_D7 = as.numeric(sga_adj_lfc_D4_D7)) %>% mutate(sgb_adj_lfc_D4_D7 = as.numeric(sgb_adj_lfc_D4_D7)) %>% ggscatter(x="combined_LFC_D4_D7",y="sgb_adj_lfc_D4_D7", color = "PAM_status_simplified", alpha= 0.35, size=2, palette = c("black", "red", "green"),  add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D10 != "NA") %>% filter(sgb_adj_lfc_D4_D10!="NA") %>% mutate(sga_adj_lfc_D4_D10 = as.numeric(sga_adj_lfc_D4_D10)) %>% mutate(sgb_adj_lfc_D4_D10 = as.numeric(sgb_adj_lfc_D4_D10)) %>% ggscatter(x="sga_adj_lfc_D4_D10",y="sgb_adj_lfc_D4_D10", color = "PAM_status_simplified", alpha= 0.35, size=2, palette = c("black", "red", "green"),  add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D10 != "NA") %>% filter(sgb_adj_lfc_D4_D10!="NA") %>% mutate(sga_adj_lfc_D4_D10 = as.numeric(sga_adj_lfc_D4_D10)) %>% mutate(sgb_adj_lfc_D4_D10 = as.numeric(sgb_adj_lfc_D4_D10)) %>% ggscatter(x="sga_adj_lfc_D4_D10",y="combined_LFC_D4_D10", color = "PAM_status_simplified", alpha= 0.35, size=2, palette = c("black", "red", "green"),  add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D10 != "NA") %>% filter(sgb_adj_lfc_D4_D10!="NA") %>% mutate(sga_adj_lfc_D4_D10 = as.numeric(sga_adj_lfc_D4_D10)) %>% mutate(sgb_adj_lfc_D4_D10 = as.numeric(sgb_adj_lfc_D4_D10)) %>% ggscatter(x="combined_LFC_D4_D10",y="sgb_adj_lfc_D4_D10", color = "PAM_status_simplified", alpha= 0.35, size=2, palette = c("black", "red", "green"),  add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D14 != "NA") %>% filter(sgb_adj_lfc_D4_D14!="NA") %>% mutate(sga_adj_lfc_D4_D14 = as.numeric(sga_adj_lfc_D4_D14)) %>% mutate(sgb_adj_lfc_D4_D14 = as.numeric(sgb_adj_lfc_D4_D14)) %>% ggscatter(x="sga_adj_lfc_D4_D14",y="sgb_adj_lfc_D4_D14", color = "PAM_status_simplified", alpha= 0.35, size=2, palette = c("black", "red", "green"),  add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D14 != "NA") %>% filter(sgb_adj_lfc_D4_D14!="NA") %>% mutate(sga_adj_lfc_D4_D14 = as.numeric(sga_adj_lfc_D4_D14)) %>% mutate(sgb_adj_lfc_D4_D14 = as.numeric(sgb_adj_lfc_D4_D14)) %>% ggscatter(x="sga_adj_lfc_D4_D14",y="combined_LFC_D4_D14", color = "PAM_status_simplified", alpha= 0.35, size=2, palette = c("black", "red", "green"),  add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D14 != "NA") %>% filter(sgb_adj_lfc_D4_D14!="NA") %>% mutate(sga_adj_lfc_D4_D14 = as.numeric(sga_adj_lfc_D4_D14)) %>% mutate(sgb_adj_lfc_D4_D14 = as.numeric(sgb_adj_lfc_D4_D14)) %>% ggscatter(x="combined_LFC_D4_D14",y="sgb_adj_lfc_D4_D14", color = "PAM_status_simplified", alpha= 0.35, size=2, palette = c("black", "red", "green"),  add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D21 != "NA") %>% filter(sgb_adj_lfc_D4_D21!="NA") %>% mutate(sga_adj_lfc_D4_D21 = as.numeric(sga_adj_lfc_D4_D21)) %>% mutate(sgb_adj_lfc_D4_D21 = as.numeric(sgb_adj_lfc_D4_D21)) %>% ggscatter(x="sga_adj_lfc_D4_D21",y="sgb_adj_lfc_D4_D21", color = "PAM_status_simplified", alpha= 0.35, size=2, palette = c("black", "red", "green"),  add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D21 != "NA") %>% filter(sgb_adj_lfc_D4_D21!="NA") %>% mutate(sga_adj_lfc_D4_D21 = as.numeric(sga_adj_lfc_D4_D21)) %>% mutate(sgb_adj_lfc_D4_D21 = as.numeric(sgb_adj_lfc_D4_D21)) %>% ggscatter(x="sga_adj_lfc_D4_D21",y="combined_LFC_D4_D21", color = "PAM_status_simplified", alpha= 0.35, size=2, palette = c("black", "red", "green"),  add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D21 != "NA") %>% filter(sgb_adj_lfc_D4_D21!="NA") %>% mutate(sga_adj_lfc_D4_D21 = as.numeric(sga_adj_lfc_D4_D21)) %>% mutate(sgb_adj_lfc_D4_D21 = as.numeric(sgb_adj_lfc_D4_D21)) %>% ggscatter(x="combined_LFC_D4_D21",y="sgb_adj_lfc_D4_D21", color = "PAM_status_simplified", alpha= 0.35, size=2, palette = c("black", "red", "green"),  add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_continuous != "NA") %>% filter(combined_LFC_continuous!="NA") %>% mutate(sga_adj_lfc_continuous = as.numeric(sga_adj_lfc_continuous)) %>% mutate(combined_LFC_continuous = as.numeric(combined_LFC_continuous)) %>% ggscatter(x="combined_LFC_continuous",y="sga_adj_lfc_continuous", color = "PAM_status_simplified", alpha= 0.35, size=2, palette = c("black", "red", "green"),  add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sgb_adj_lfc_continuous != "NA") %>% filter(combined_LFC_continuous!="NA") %>% mutate(sgb_adj_lfc_continuous = as.numeric(sgb_adj_lfc_continuous)) %>% mutate(combined_LFC_continuous = as.numeric(combined_LFC_continuous)) %>% ggscatter(x="sgb_adj_lfc_continuous",y="combined_LFC_continuous", color = "PAM_status_simplified", alpha= 0.35, size=2, palette = c("black", "red", "green"),  add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
dev.off()



#same plot as above but coloured by exon according to VEP annotation
pdf("correlated_adj_LFC_sga_sgb_combined2_EXON.pdf")
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D7 != "NA") %>% filter(sgb_adj_lfc_D4_D7!="NA") %>% mutate(sga_adj_lfc_D4_D7 = as.numeric(sga_adj_lfc_D4_D7)) %>% mutate(sgb_adj_lfc_D4_D7 = as.numeric(sgb_adj_lfc_D4_D7)) %>% ggscatter(x="sga_adj_lfc_D4_D7",y="sgb_adj_lfc_D4_D7", color = "EXON_GROUP", shape = "PAM_status_simplified", alpha= 0.35, size=2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n")) %>% ggpar(font.legend = c(6, "plain", "black"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D7 != "NA") %>% filter(sgb_adj_lfc_D4_D7!="NA") %>% mutate(sga_adj_lfc_D4_D7 = as.numeric(sga_adj_lfc_D4_D7)) %>% mutate(sgb_adj_lfc_D4_D7 = as.numeric(sgb_adj_lfc_D4_D7)) %>% ggscatter(x="sga_adj_lfc_D4_D7",y="combined_LFC_D4_D7", color = "EXON_GROUP", shape = "PAM_status_simplified", alpha= 0.35, size=2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n")) %>% ggpar(font.legend = c(6, "plain", "black"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D7 != "NA") %>% filter(sgb_adj_lfc_D4_D7!="NA") %>% mutate(sga_adj_lfc_D4_D7 = as.numeric(sga_adj_lfc_D4_D7)) %>% mutate(sgb_adj_lfc_D4_D7 = as.numeric(sgb_adj_lfc_D4_D7)) %>% ggscatter(x="combined_LFC_D4_D7",y="sgb_adj_lfc_D4_D7", color = "EXON_GROUP", shape = "PAM_status_simplified", alpha= 0.35, size=2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))  %>% ggpar(font.legend = c(6, "plain", "black"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D10 != "NA") %>% filter(sgb_adj_lfc_D4_D10!="NA") %>% mutate(sga_adj_lfc_D4_D10 = as.numeric(sga_adj_lfc_D4_D10)) %>% mutate(sgb_adj_lfc_D4_D10 = as.numeric(sgb_adj_lfc_D4_D10)) %>% ggscatter(x="sga_adj_lfc_D4_D10",y="sgb_adj_lfc_D4_D10", color = "EXON_GROUP", shape = "PAM_status_simplified", alpha= 0.35, size=2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))  %>% ggpar(font.legend = c(6, "plain", "black"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D10 != "NA") %>% filter(sgb_adj_lfc_D4_D10!="NA") %>% mutate(sga_adj_lfc_D4_D10 = as.numeric(sga_adj_lfc_D4_D10)) %>% mutate(sgb_adj_lfc_D4_D10 = as.numeric(sgb_adj_lfc_D4_D10)) %>% ggscatter(x="sga_adj_lfc_D4_D10",y="combined_LFC_D4_D10", color = "EXON_GROUP", shape = "PAM_status_simplified", alpha= 0.35, size=2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))  %>% ggpar(font.legend = c(6, "plain", "black"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D10 != "NA") %>% filter(sgb_adj_lfc_D4_D10!="NA") %>% mutate(sga_adj_lfc_D4_D10 = as.numeric(sga_adj_lfc_D4_D10)) %>% mutate(sgb_adj_lfc_D4_D10 = as.numeric(sgb_adj_lfc_D4_D10)) %>% ggscatter(x="combined_LFC_D4_D10",y="sgb_adj_lfc_D4_D10", color = "EXON_GROUP", shape = "PAM_status_simplified", alpha= 0.35, size=2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))  %>% ggpar(font.legend = c(6, "plain", "black"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D14 != "NA") %>% filter(sgb_adj_lfc_D4_D14!="NA") %>% mutate(sga_adj_lfc_D4_D14 = as.numeric(sga_adj_lfc_D4_D14)) %>% mutate(sgb_adj_lfc_D4_D14 = as.numeric(sgb_adj_lfc_D4_D14)) %>% ggscatter(x="sga_adj_lfc_D4_D14",y="sgb_adj_lfc_D4_D14", color = "EXON_GROUP", shape = "PAM_status_simplified", alpha= 0.35, size=2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))  %>% ggpar(font.legend = c(6, "plain", "black"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D14 != "NA") %>% filter(sgb_adj_lfc_D4_D14!="NA") %>% mutate(sga_adj_lfc_D4_D14 = as.numeric(sga_adj_lfc_D4_D14)) %>% mutate(sgb_adj_lfc_D4_D14 = as.numeric(sgb_adj_lfc_D4_D14)) %>% ggscatter(x="sga_adj_lfc_D4_D14",y="combined_LFC_D4_D14", color = "EXON_GROUP", shape = "PAM_status_simplified", alpha= 0.35, size=2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))  %>% ggpar(font.legend = c(6, "plain", "black"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D14 != "NA") %>% filter(sgb_adj_lfc_D4_D14!="NA") %>% mutate(sga_adj_lfc_D4_D14 = as.numeric(sga_adj_lfc_D4_D14)) %>% mutate(sgb_adj_lfc_D4_D14 = as.numeric(sgb_adj_lfc_D4_D14)) %>% ggscatter(x="combined_LFC_D4_D14",y="sgb_adj_lfc_D4_D14", color = "EXON_GROUP", shape = "PAM_status_simplified", alpha= 0.35, size=2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))  %>% ggpar(font.legend = c(6, "plain", "black"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D21 != "NA") %>% filter(sgb_adj_lfc_D4_D21!="NA") %>% mutate(sga_adj_lfc_D4_D21 = as.numeric(sga_adj_lfc_D4_D21)) %>% mutate(sgb_adj_lfc_D4_D21 = as.numeric(sgb_adj_lfc_D4_D21)) %>% ggscatter(x="sga_adj_lfc_D4_D21",y="sgb_adj_lfc_D4_D21", color = "EXON_GROUP", shape = "PAM_status_simplified", alpha= 0.35, size=2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))  %>% ggpar(font.legend = c(6, "plain", "black"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D21 != "NA") %>% filter(sgb_adj_lfc_D4_D21!="NA") %>% mutate(sga_adj_lfc_D4_D21 = as.numeric(sga_adj_lfc_D4_D21)) %>% mutate(sgb_adj_lfc_D4_D21 = as.numeric(sgb_adj_lfc_D4_D21)) %>% ggscatter(x="sga_adj_lfc_D4_D21",y="combined_LFC_D4_D21", color = "EXON_GROUP", shape = "PAM_status_simplified", alpha= 0.35, size=2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))  %>% ggpar(font.legend = c(6, "plain", "black"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_D4_D21 != "NA") %>% filter(sgb_adj_lfc_D4_D21!="NA") %>% mutate(sga_adj_lfc_D4_D21 = as.numeric(sga_adj_lfc_D4_D21)) %>% mutate(sgb_adj_lfc_D4_D21 = as.numeric(sgb_adj_lfc_D4_D21)) %>% ggscatter(x="combined_LFC_D4_D21",y="sgb_adj_lfc_D4_D21", color = "EXON_GROUP", shape = "PAM_status_simplified", alpha= 0.35, size=2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))  %>% ggpar(font.legend = c(6, "plain", "black"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sga_adj_lfc_continuous != "NA") %>% filter(combined_LFC_continuous!="NA") %>% mutate(sga_adj_lfc_continuous = as.numeric(sga_adj_lfc_continuous)) %>% mutate(combined_LFC_continuous = as.numeric(combined_LFC_continuous)) %>% ggscatter(x="combined_LFC_continuous",y="sga_adj_lfc_continuous", color = "EXON_GROUP", shape = "PAM_status_simplified", alpha= 0.35, size=2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))  %>% ggpar(font.legend = c(6, "plain", "black"))
sga_sgb_combined_guide %>% filter(Variant_duplication=="2") %>% filter(sgb_adj_lfc_continuous != "NA") %>% filter(combined_LFC_continuous!="NA") %>% mutate(sgb_adj_lfc_continuous = as.numeric(sgb_adj_lfc_continuous)) %>% mutate(combined_LFC_continuous = as.numeric(combined_LFC_continuous)) %>% ggscatter(x="sgb_adj_lfc_continuous",y="combined_LFC_continuous", color = "EXON_GROUP", shape = "PAM_status_simplified", alpha= 0.35, size=2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))  %>% ggpar(font.legend = c(6, "plain", "black"))
dev.off()


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
combined_sub <- sga_sgb_combined2[ ,c(1:8,56:60,72,74,76,78,80:90)]
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
uncombined_sub <- total_single_bap1[ ,c(1:4,6,7,10,14,18,22,42,50:59,60,62,64,66,68)]
#add a column and populate with description of this data as being from an uncombined source, that is one guide was used to get this data
uncombined_sub$Variant_duplication<-NA
uncombined_sub$dataset_process<-"not_combined"
#change the names of the following columns, so that they can be bound to the combined dataframe
uncombined_sub<-uncombined_sub %>% dplyr::rename(Variant_Sources = "SG")
uncombined_sub<-uncombined_sub %>% dplyr::rename(PAM_status = "pam_mut_sgrna_id")
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
total_bap1_combined_uncombined <- do.call("rbind.fill", list(combined_sub,uncombined_sub))

#output the dataframe to check - this will be input dataframe to combined tiled regions
write.csv(total_bap1_combined_uncombined,"./total_bap1_combined_uncombined.csv", row.names = FALSE)

#can 'uncomment' the line below to check how many unique ids are going to be processed - ids will be duplicated between targetons as we have grouped by id AND targeton so far 
                #total_bap1_combined_uncombined_count<-total_bap1_combined_uncombined[!duplicated(total_bap1_combined_uncombined$id), ]

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

#NOTES START 
# the previous sectrions produce a dataframe that contains combined and uncombined variant counts
# many mseqs will be duplicated even though the ids are unique 
# however even the ids are not unique as they are still grouped by id AND targeton from which they came, so overlapping targetons will contain duplicates
# we now need to combine these overlapping regions in a way similar to the weighting of the two sources of sgRNAs
# as have already accounted for the PAM sites in the combination of guides we do not need to worry about PAM sites now
# but as variants could be derived from library a or library b alone we still need to annotate with the appropriate metadata 
#NOTES END

#First step: we filter the inherited dataset (the output from the above processes) to obtain variants that are not NA in mseq_a 
#(some will be because they are from the unconbined source and from library b)
is_not_NA_mseq_a<- total_bap1_combined_uncombined %>% filter(!mseq_a %in% "NA") %>% filter(!is.na(mseq_a))
#group by these non-NA mseq_a sequences to make them unique
#we also create two new fields dup_mseqa, which is the number of instances that mseq_a is seen
#also dup_PAM_status_mseqa - this is similar to arrayed_mseq_a seen for the combination of guides, except here is captures and vectorizes each pam_status for each mseqa
#this is important because some mseq_a duplicates will not be listed as having a pam status and other identical sequences will be, particular to 1bp deletions of repetitive sequence
is_not_NA_mseq_a_unique <- is_not_NA_mseq_a  %>% 
  group_by(mseq_a) %>%
  dplyr::summarize(dup_mseqa= n(), 
                   dup_PAM_status_mseqa = paste(PAM_status, collapse = ','))
#join back to inherited dataframe to get the data
is_not_NA_mseq_a_unique <- inner_join(is_not_NA_mseq_a_unique, total_bap1_combined_uncombined, by="mseq_a", all.x=FALSE)
#now it's safe to make unique for mseq_a
is_not_NA_mseq_a_unique <- is_not_NA_mseq_a_unique[!duplicated(is_not_NA_mseq_a_unique$mseq_a), ]
#now get the NA for mseqa
is_NA_mseq_a<- total_bap1_combined_uncombined %>% filter(mseq_a %in% "NA")
#now it's safe to make unique for mseq_b
is_NA_mseq_a <- is_NA_mseq_a[!duplicated(is_NA_mseq_a$mseq_b), ]
#fill the new fields appropriately 
is_NA_mseq_a$dup_mseqa<-NA
is_NA_mseq_a$dup_PAM_status_mseqa<-NA
#combine non-NA unique for mseq_a and NA mseq_b
total_mseq_a_unique<-rbind(is_not_NA_mseq_a_unique,is_NA_mseq_a)
#we cannot make unique for mseq_b and mseq_a together as some rows will be lost
#RENAME the same dataframe 
total_mseq_unique<-total_mseq_a_unique



#ids are not in themselves unique - each id can be linked to multiple mseqs
#each mseq can be linked to multiple ids
#when u make mseq_a unique it will link back to multiple ids
#when u make mseq_a unique and annotion b is required, mseq_b for that id could be multiple sequences which each link back to multiple ids 
#this is an issue in annotation - how do you solve it?
#required annotation id may beed to be refined
#make sure the arrayed pam status includes arrayed for mseq b
#do the above grouping for mseqa with mseqb also and see the interesction
#the targeton processing should have made everything unique for id - you make it unique for mseq a, then a bit for mseqb then you group by id
#when u do the targeton processing - think about this....




#can check mseq_a is not duplciated and look at mseqb

duplicated_b <- total_mseq_unique  %>% 
  group_by(mseq_b) %>%
  filter(n() > 1) %>%
  dplyr::summarize(dup_mseqb= n(),
                   Variant_duplication=Variant_duplication,
                   Variant_Sources=Variant_Sources,
                   mseq_a=mseq_a,
                   dup_PAM_status_mseqb = paste(PAM_status, collapse = ','))




                    duplicated_b<-total_mseq_unique %>%
                      group_by(mseq_b) %>%
                      filter(n() > 1) %>%
                     ungroup
                    duplicated_b <- duplicated_b[!duplicated(duplicated_b$mseq_b), ]
                    
################# COMBINATION OF TARGETONS ############# START
################# COMBINATION OF TARGETONS ############# START
################# COMBINATION OF TARGETONS ############# START
################# COMBINATION OF TARGETONS ############# START
################# COMBINATION OF TARGETONS ############# START
################# COMBINATION OF TARGETONS ############# START

#group the targetons by id - this collapse the different mseq_a categories that assess the same variant from different (overlapping) targetons
                    
                    
                    
                    
  targeton_collapse <- total_mseq_unique %>%
  group_by(id) %>%
  dplyr::summarize(variant_n_instances_between_targetons= n(), 
                   exon_specific_names = paste(combined_name, collapse = '|'),
                   source_dataset = paste(dataset_process, collapse ="|"),
                   source_targetons = paste(EXON_GROUP, collapse = '|'),
                   source_libraries = paste(Variant_Sources, collapse = '|'),
                   variant_n_instances_between_libraries = paste(Variant_duplication, collapse = '|'),
                   arrayed_PAM_status = paste(dup_PAM_status_mseqa, collapse = '|'), 
                   single_PAM_status = paste(PAM_status, collapse='|'),
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
targeton_collapse[, c(12:31)] <- sapply(targeton_collapse[, c(12:31)], as.numeric)
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
  mutate(processed_LFC_continuous=case_when(str_detect(variant_n_instances_between_targetons,"1")~primary_tg_adj_lfc_continuous, TRUE~ sum_of_weighted_LFC_continuous/sum_of_weight_continuous)) %>%
  
  #PERFORM THE CALCULATION OF Z SCORES > P VALUES > BH_FDR - this is the final statistics that will be used 
  mutate(processed_Z_D4_D7=processed_LFC_D4_D7/SE_bind_D4_D7) %>% mutate(two_tailed_p_D4_D7= pnorm(abs(processed_Z_D4_D7),lower.tail = FALSE) *2) %>% mutate(processed_BH_FDR_D4_D7 = p.adjust(two_tailed_p_D4_D7, method = "BH")) %>%
  mutate(processed_Z_D4_D10=processed_LFC_D4_D10/SE_bind_D4_D10) %>% mutate(two_tailed_p_D4_D10= pnorm(abs(processed_Z_D4_D10),lower.tail = FALSE) *2) %>% mutate(processed_BH_FDR_D4_D10 = p.adjust(two_tailed_p_D4_D10, method = "BH")) %>%
  mutate(processed_Z_D4_D14=processed_LFC_D4_D14/SE_bind_D4_D14) %>% mutate(two_tailed_p_D4_D14= pnorm(abs(processed_Z_D4_D14),lower.tail = FALSE) *2) %>% mutate(processed_BH_FDR_D4_D14 = p.adjust(two_tailed_p_D4_D14, method = "BH")) %>%
  mutate(processed_Z_D4_D21=processed_LFC_D4_D21/SE_bind_D4_D21) %>% mutate(two_tailed_p_D4_D21= pnorm(abs(processed_Z_D4_D21),lower.tail = FALSE) *2) %>% mutate(processed_BH_FDR_D4_D21 = p.adjust(two_tailed_p_D4_D21, method = "BH")) %>%
  mutate(processed_Z_continuous=processed_LFC_continuous/SE_bind_continuous) %>% mutate(two_tailed_p_continuous= pnorm(abs(processed_Z_continuous),lower.tail = FALSE) *2) %>% mutate(processed_BH_FDR_continuous = p.adjust(two_tailed_p_continuous, method = "BH"))

#WRITE OUT this dataframe - there are many fields that are not necessary to keep but good for checking
write.csv(targeton_processed,"./sga_sgb_combined_and_targeton_processed.csv", row.names = FALSE)

##########QC combination of targetons #################### START
#check correlation after combination of targeton overlapping reions
pdf("targeton_overlap_correlation.pdf")
targeton_processed %>% filter(variant_n_instances_between_targetons=="2") %>% filter(primary_tg_adj_lfc_D4_D7 != "NA") %>% filter(secondary_tg_adj_lfc_D4_D7!="NA") %>% mutate(primary_tg_adj_lfc_D4_D7 = as.numeric(primary_tg_adj_lfc_D4_D7)) %>% mutate(secondary_tg_adj_lfc_D4_D7 = as.numeric(secondary_tg_adj_lfc_D4_D7)) %>% ggscatter(x="primary_tg_adj_lfc_D4_D7",y="secondary_tg_adj_lfc_D4_D7",  alpha= 0.35, size=2, palette = c("black", "red", "green"), add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
targeton_processed %>% filter(variant_n_instances_between_targetons=="2") %>% filter(primary_tg_adj_lfc_D4_D7 != "NA") %>% filter(secondary_tg_adj_lfc_D4_D7!="NA") %>% mutate(primary_tg_adj_lfc_D4_D7 = as.numeric(primary_tg_adj_lfc_D4_D7)) %>% mutate(secondary_tg_adj_lfc_D4_D7 = as.numeric(secondary_tg_adj_lfc_D4_D7)) %>% ggscatter(x="primary_tg_adj_lfc_D4_D7",y="processed_LFC_D4_D7",  alpha= 0.35, size=2, palette = c("black", "red", "green"), add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
targeton_processed %>% filter(variant_n_instances_between_targetons=="2") %>% filter(primary_tg_adj_lfc_D4_D7 != "NA") %>% filter(secondary_tg_adj_lfc_D4_D7!="NA") %>% mutate(primary_tg_adj_lfc_D4_D7 = as.numeric(primary_tg_adj_lfc_D4_D7)) %>% mutate(secondary_tg_adj_lfc_D4_D7 = as.numeric(secondary_tg_adj_lfc_D4_D7)) %>% ggscatter(x="processed_LFC_D4_D7",y="secondary_tg_adj_lfc_D4_D7",  alpha= 0.35, size=2, palette = c("black", "red", "green"), add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
targeton_processed %>% filter(variant_n_instances_between_targetons=="2") %>% filter(primary_tg_adj_lfc_D4_D10 != "NA") %>% filter(secondary_tg_adj_lfc_D4_D10!="NA") %>% mutate(primary_tg_adj_lfc_D4_D10 = as.numeric(primary_tg_adj_lfc_D4_D10)) %>% mutate(secondary_tg_adj_lfc_D4_D10 = as.numeric(secondary_tg_adj_lfc_D4_D10)) %>% ggscatter(x="primary_tg_adj_lfc_D4_D10",y="secondary_tg_adj_lfc_D4_D10",  alpha= 0.35, size=2, palette = c("black", "red", "green"), add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
targeton_processed %>% filter(variant_n_instances_between_targetons=="2") %>% filter(primary_tg_adj_lfc_D4_D10 != "NA") %>% filter(secondary_tg_adj_lfc_D4_D10!="NA") %>% mutate(primary_tg_adj_lfc_D4_D10 = as.numeric(primary_tg_adj_lfc_D4_D10)) %>% mutate(secondary_tg_adj_lfc_D4_D10 = as.numeric(secondary_tg_adj_lfc_D4_D10)) %>% ggscatter(x="primary_tg_adj_lfc_D4_D10",y="processed_LFC_D4_D10",  alpha= 0.35, size=2, palette = c("black", "red", "green"), add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
targeton_processed %>% filter(variant_n_instances_between_targetons=="2") %>% filter(primary_tg_adj_lfc_D4_D10 != "NA") %>% filter(secondary_tg_adj_lfc_D4_D10!="NA") %>% mutate(primary_tg_adj_lfc_D4_D10 = as.numeric(primary_tg_adj_lfc_D4_D10)) %>% mutate(secondary_tg_adj_lfc_D4_D10 = as.numeric(secondary_tg_adj_lfc_D4_D10)) %>% ggscatter(x="processed_LFC_D4_D10",y="secondary_tg_adj_lfc_D4_D10",  alpha= 0.35, size=2, palette = c("black", "red", "green"), add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
targeton_processed %>% filter(variant_n_instances_between_targetons=="2") %>% filter(primary_tg_adj_lfc_D4_D14 != "NA") %>% filter(secondary_tg_adj_lfc_D4_D14!="NA") %>% mutate(primary_tg_adj_lfc_D4_D14 = as.numeric(primary_tg_adj_lfc_D4_D14)) %>% mutate(secondary_tg_adj_lfc_D4_D14 = as.numeric(secondary_tg_adj_lfc_D4_D14)) %>% ggscatter(x="primary_tg_adj_lfc_D4_D14",y="secondary_tg_adj_lfc_D4_D14",  alpha= 0.35, size=2, palette = c("black", "red", "green"), add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
targeton_processed %>% filter(variant_n_instances_between_targetons=="2") %>% filter(primary_tg_adj_lfc_D4_D14 != "NA") %>% filter(secondary_tg_adj_lfc_D4_D14!="NA") %>% mutate(primary_tg_adj_lfc_D4_D14 = as.numeric(primary_tg_adj_lfc_D4_D14)) %>% mutate(secondary_tg_adj_lfc_D4_D14 = as.numeric(secondary_tg_adj_lfc_D4_D14)) %>% ggscatter(x="primary_tg_adj_lfc_D4_D14",y="processed_LFC_D4_D14",  alpha= 0.35, size=2, palette = c("black", "red", "green"), add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
targeton_processed %>% filter(variant_n_instances_between_targetons=="2") %>% filter(primary_tg_adj_lfc_D4_D14 != "NA") %>% filter(secondary_tg_adj_lfc_D4_D14!="NA") %>% mutate(primary_tg_adj_lfc_D4_D14 = as.numeric(primary_tg_adj_lfc_D4_D14)) %>% mutate(secondary_tg_adj_lfc_D4_D14 = as.numeric(secondary_tg_adj_lfc_D4_D14)) %>% ggscatter(x="processed_LFC_D4_D14",y="secondary_tg_adj_lfc_D4_D14",  alpha= 0.35, size=2, palette = c("black", "red", "green"), add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
targeton_processed %>% filter(variant_n_instances_between_targetons=="2") %>% filter(primary_tg_adj_lfc_D4_D21 != "NA") %>% filter(secondary_tg_adj_lfc_D4_D21!="NA") %>% mutate(primary_tg_adj_lfc_D4_D21 = as.numeric(primary_tg_adj_lfc_D4_D21)) %>% mutate(secondary_tg_adj_lfc_D4_D21 = as.numeric(secondary_tg_adj_lfc_D4_D21)) %>% ggscatter(x="primary_tg_adj_lfc_D4_D21",y="secondary_tg_adj_lfc_D4_D21",  alpha= 0.35, size=2, palette = c("black", "red", "green"), add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
targeton_processed %>% filter(variant_n_instances_between_targetons=="2") %>% filter(primary_tg_adj_lfc_D4_D21 != "NA") %>% filter(secondary_tg_adj_lfc_D4_D21!="NA") %>% mutate(primary_tg_adj_lfc_D4_D21 = as.numeric(primary_tg_adj_lfc_D4_D21)) %>% mutate(secondary_tg_adj_lfc_D4_D21 = as.numeric(secondary_tg_adj_lfc_D4_D21)) %>% ggscatter(x="primary_tg_adj_lfc_D4_D21",y="processed_LFC_D4_D21",  alpha= 0.35, size=2, palette = c("black", "red", "green"), add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
targeton_processed %>% filter(variant_n_instances_between_targetons=="2") %>% filter(primary_tg_adj_lfc_D4_D21 != "NA") %>% filter(secondary_tg_adj_lfc_D4_D21!="NA") %>% mutate(primary_tg_adj_lfc_D4_D21 = as.numeric(primary_tg_adj_lfc_D4_D21)) %>% mutate(secondary_tg_adj_lfc_D4_D21 = as.numeric(secondary_tg_adj_lfc_D4_D21)) %>% ggscatter(x="processed_LFC_D4_D21",y="secondary_tg_adj_lfc_D4_D21",  alpha= 0.35, size=2, palette = c("black", "red", "green"), add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
targeton_processed %>% filter(variant_n_instances_between_targetons=="2") %>% filter(primary_tg_adj_lfc_continuous != "NA") %>% filter(secondary_tg_adj_lfc_continuous!="NA") %>% mutate(primary_tg_adj_lfc_continuous = as.numeric(primary_tg_adj_lfc_continuous)) %>% mutate(secondary_tg_adj_lfc_continuous = as.numeric(secondary_tg_adj_lfc_continuous)) %>% ggscatter(x="processed_LFC_continuous",y="primary_tg_adj_lfc_continuous",  alpha= 0.35, size=2, palette = c("black", "red", "green"), add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
targeton_processed %>% filter(variant_n_instances_between_targetons=="2") %>% filter(primary_tg_adj_lfc_continuous != "NA") %>% filter(secondary_tg_adj_lfc_continuous!="NA") %>% mutate(primary_tg_adj_lfc_continuous = as.numeric(primary_tg_adj_lfc_continuous)) %>% mutate(secondary_tg_adj_lfc_continuous = as.numeric(secondary_tg_adj_lfc_continuous)) %>% ggscatter(x="secondary_tg_adj_lfc_continuous",y="processed_LFC_continuous", alpha= 0.35, size=2, palette = c("black", "red", "green"), add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
dev.off()

#check correlation after combination of targeton overlapping reions - coloured by exon 
pdf("targeton_overlap_correlation_EXON.pdf")
targeton_processed %>% filter(variant_n_instances_between_targetons=="2") %>% filter(primary_tg_adj_lfc_D4_D7 != "NA") %>% filter(secondary_tg_adj_lfc_D4_D7!="NA") %>% mutate(primary_tg_adj_lfc_D4_D7 = as.numeric(primary_tg_adj_lfc_D4_D7)) %>% mutate(secondary_tg_adj_lfc_D4_D7 = as.numeric(secondary_tg_adj_lfc_D4_D7)) %>% ggscatter(x="primary_tg_adj_lfc_D4_D7",y="secondary_tg_adj_lfc_D4_D7",  color = "collected_EXON_GROUP", alpha= 0.35, size=2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))  %>% ggpar(font.legend = c(6, "plain", "black"))
targeton_processed %>% filter(variant_n_instances_between_targetons=="2") %>% filter(primary_tg_adj_lfc_D4_D7 != "NA") %>% filter(secondary_tg_adj_lfc_D4_D7!="NA") %>% mutate(primary_tg_adj_lfc_D4_D7 = as.numeric(primary_tg_adj_lfc_D4_D7)) %>% mutate(secondary_tg_adj_lfc_D4_D7 = as.numeric(secondary_tg_adj_lfc_D4_D7)) %>% ggscatter(x="primary_tg_adj_lfc_D4_D7",y="processed_LFC_D4_D7",  color = "collected_EXON_GROUP", alpha= 0.35, size=2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))  %>% ggpar(font.legend = c(6, "plain", "black"))
targeton_processed %>% filter(variant_n_instances_between_targetons=="2") %>% filter(primary_tg_adj_lfc_D4_D7 != "NA") %>% filter(secondary_tg_adj_lfc_D4_D7!="NA") %>% mutate(primary_tg_adj_lfc_D4_D7 = as.numeric(primary_tg_adj_lfc_D4_D7)) %>% mutate(secondary_tg_adj_lfc_D4_D7 = as.numeric(secondary_tg_adj_lfc_D4_D7)) %>% ggscatter(x="processed_LFC_D4_D7",y="secondary_tg_adj_lfc_D4_D7",  color = "collected_EXON_GROUP", alpha= 0.35, size=2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))  %>% ggpar(font.legend = c(6, "plain", "black"))
targeton_processed %>% filter(variant_n_instances_between_targetons=="2") %>% filter(primary_tg_adj_lfc_D4_D10 != "NA") %>% filter(secondary_tg_adj_lfc_D4_D10!="NA") %>% mutate(primary_tg_adj_lfc_D4_D10 = as.numeric(primary_tg_adj_lfc_D4_D10)) %>% mutate(secondary_tg_adj_lfc_D4_D10 = as.numeric(secondary_tg_adj_lfc_D4_D10)) %>% ggscatter(x="primary_tg_adj_lfc_D4_D10",y="secondary_tg_adj_lfc_D4_D10",  color = "collected_EXON_GROUP", alpha= 0.35, size=2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))  %>% ggpar(font.legend = c(6, "plain", "black"))
targeton_processed %>% filter(variant_n_instances_between_targetons=="2") %>% filter(primary_tg_adj_lfc_D4_D10 != "NA") %>% filter(secondary_tg_adj_lfc_D4_D10!="NA") %>% mutate(primary_tg_adj_lfc_D4_D10 = as.numeric(primary_tg_adj_lfc_D4_D10)) %>% mutate(secondary_tg_adj_lfc_D4_D10 = as.numeric(secondary_tg_adj_lfc_D4_D10)) %>% ggscatter(x="primary_tg_adj_lfc_D4_D10",y="processed_LFC_D4_D10",  color = "collected_EXON_GROUP", alpha= 0.35, size=2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))  %>% ggpar(font.legend = c(6, "plain", "black"))
targeton_processed %>% filter(variant_n_instances_between_targetons=="2") %>% filter(primary_tg_adj_lfc_D4_D10 != "NA") %>% filter(secondary_tg_adj_lfc_D4_D10!="NA") %>% mutate(primary_tg_adj_lfc_D4_D10 = as.numeric(primary_tg_adj_lfc_D4_D10)) %>% mutate(secondary_tg_adj_lfc_D4_D10 = as.numeric(secondary_tg_adj_lfc_D4_D10)) %>% ggscatter(x="processed_LFC_D4_D10",y="secondary_tg_adj_lfc_D4_D10",  color = "collected_EXON_GROUP", alpha= 0.35, size=2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))  %>% ggpar(font.legend = c(6, "plain", "black"))
targeton_processed %>% filter(variant_n_instances_between_targetons=="2") %>% filter(primary_tg_adj_lfc_D4_D14 != "NA") %>% filter(secondary_tg_adj_lfc_D4_D14!="NA") %>% mutate(primary_tg_adj_lfc_D4_D14 = as.numeric(primary_tg_adj_lfc_D4_D14)) %>% mutate(secondary_tg_adj_lfc_D4_D14 = as.numeric(secondary_tg_adj_lfc_D4_D14)) %>% ggscatter(x="primary_tg_adj_lfc_D4_D14",y="secondary_tg_adj_lfc_D4_D14",  color = "collected_EXON_GROUP", alpha= 0.35, size=2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))  %>% ggpar(font.legend = c(6, "plain", "black"))
targeton_processed %>% filter(variant_n_instances_between_targetons=="2") %>% filter(primary_tg_adj_lfc_D4_D14 != "NA") %>% filter(secondary_tg_adj_lfc_D4_D14!="NA") %>% mutate(primary_tg_adj_lfc_D4_D14 = as.numeric(primary_tg_adj_lfc_D4_D14)) %>% mutate(secondary_tg_adj_lfc_D4_D14 = as.numeric(secondary_tg_adj_lfc_D4_D14)) %>% ggscatter(x="primary_tg_adj_lfc_D4_D14",y="processed_LFC_D4_D14",  color = "collected_EXON_GROUP", alpha= 0.35, size=2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))  %>% ggpar(font.legend = c(6, "plain", "black"))
targeton_processed %>% filter(variant_n_instances_between_targetons=="2") %>% filter(primary_tg_adj_lfc_D4_D14 != "NA") %>% filter(secondary_tg_adj_lfc_D4_D14!="NA") %>% mutate(primary_tg_adj_lfc_D4_D14 = as.numeric(primary_tg_adj_lfc_D4_D14)) %>% mutate(secondary_tg_adj_lfc_D4_D14 = as.numeric(secondary_tg_adj_lfc_D4_D14)) %>% ggscatter(x="processed_LFC_D4_D14",y="secondary_tg_adj_lfc_D4_D14",  color = "collected_EXON_GROUP", alpha= 0.35, size=2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))  %>% ggpar(font.legend = c(6, "plain", "black"))
targeton_processed %>% filter(variant_n_instances_between_targetons=="2") %>% filter(primary_tg_adj_lfc_D4_D21 != "NA") %>% filter(secondary_tg_adj_lfc_D4_D21!="NA") %>% mutate(primary_tg_adj_lfc_D4_D21 = as.numeric(primary_tg_adj_lfc_D4_D21)) %>% mutate(secondary_tg_adj_lfc_D4_D21 = as.numeric(secondary_tg_adj_lfc_D4_D21)) %>% ggscatter(x="primary_tg_adj_lfc_D4_D21",y="secondary_tg_adj_lfc_D4_D21",  color = "collected_EXON_GROUP", alpha= 0.35, size=2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))  %>% ggpar(font.legend = c(6, "plain", "black"))
targeton_processed %>% filter(variant_n_instances_between_targetons=="2") %>% filter(primary_tg_adj_lfc_D4_D21 != "NA") %>% filter(secondary_tg_adj_lfc_D4_D21!="NA") %>% mutate(primary_tg_adj_lfc_D4_D21 = as.numeric(primary_tg_adj_lfc_D4_D21)) %>% mutate(secondary_tg_adj_lfc_D4_D21 = as.numeric(secondary_tg_adj_lfc_D4_D21)) %>% ggscatter(x="primary_tg_adj_lfc_D4_D21",y="processed_LFC_D4_D21",  color = "collected_EXON_GROUP", alpha= 0.35, size=2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))  %>% ggpar(font.legend = c(6, "plain", "black"))
targeton_processed %>% filter(variant_n_instances_between_targetons=="2") %>% filter(primary_tg_adj_lfc_D4_D21 != "NA") %>% filter(secondary_tg_adj_lfc_D4_D21!="NA") %>% mutate(primary_tg_adj_lfc_D4_D21 = as.numeric(primary_tg_adj_lfc_D4_D21)) %>% mutate(secondary_tg_adj_lfc_D4_D21 = as.numeric(secondary_tg_adj_lfc_D4_D21)) %>% ggscatter(x="processed_LFC_D4_D21",y="secondary_tg_adj_lfc_D4_D21",  color = "collected_EXON_GROUP", alpha= 0.35, size=2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))  %>% ggpar(font.legend = c(6, "plain", "black"))
targeton_processed %>% filter(variant_n_instances_between_targetons=="2") %>% filter(primary_tg_adj_lfc_continuous != "NA") %>% filter(secondary_tg_adj_lfc_continuous!="NA") %>% mutate(primary_tg_adj_lfc_continuous = as.numeric(primary_tg_adj_lfc_continuous)) %>% mutate(secondary_tg_adj_lfc_continuous = as.numeric(secondary_tg_adj_lfc_continuous)) %>% ggscatter(x="processed_LFC_continuous",y="primary_tg_adj_lfc_continuous",  color = "collected_EXON_GROUP", alpha= 0.35, size=2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))  %>% ggpar(font.legend = c(6, "plain", "black"))
targeton_processed %>% filter(variant_n_instances_between_targetons=="2") %>% filter(primary_tg_adj_lfc_continuous != "NA") %>% filter(secondary_tg_adj_lfc_continuous!="NA") %>% mutate(primary_tg_adj_lfc_continuous = as.numeric(primary_tg_adj_lfc_continuous)) %>% mutate(secondary_tg_adj_lfc_continuous = as.numeric(secondary_tg_adj_lfc_continuous)) %>% ggscatter(x="secondary_tg_adj_lfc_continuous",y="processed_LFC_continuous",  color = "collected_EXON_GROUP", alpha= 0.35, size=2, add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))  %>% ggpar(font.legend = c(6, "plain", "black"))
dev.off()

##########QC combination of targetons #################### END

#final data cleaning and annotation
#final_metrics <- targeton_processed[ ,c(1:7,43:47,59,61,63,65,67:82)]

################# FINAL ANNOTATION ############# START
################# FINAL ANNOTATION ############# START
################# FINAL ANNOTATION ############# START
################# FINAL ANNOTATION ############# START
################# FINAL ANNOTATION ############# START
################# FINAL ANNOTATION ############# START

#annotation formatting - TAKE VALIANT OUTPUT and VEP output derrived from PAM-specific VCFs - SGA SET
E1_SGA_META_VEP<-merge(chr3_52409736_52409980_minus_sgRNA_w1_a_meta.csv,VEP_OUTPUT_chr3_52409736_52409980_minus_sgRNA_w1_a_pam.tsv, by="id", all.x=TRUE)
E2_SGA_META_VEP<-merge(chr3_52409602_52409846_minus_sgRNA_w2_a_meta.csv,VEP_OUTPUT_chr3_52409602_52409846_minus_sgRNA_w2_a_pam.tsv, by="id", all.x=TRUE)
E3_SGA_META_VEP<-merge(chr3_52409466_52409711_minus_sgRNA_w3_a_meta.csv,VEP_OUTPUT_chr3_52409466_52409711_minus_sgRNA_w3_a_pam.tsv, by="id", all.x=TRUE)
E4_SGA_META_VEP<-merge(chr3_52408417_52408661_minus_sgRNA_w4_a_meta.csv,VEP_OUTPUT_chr3_52408417_52408661_minus_sgRNA_w4_a_pam.tsv, by="id", all.x=TRUE)
E5_SGA_META_VEP<-merge(chr3_52407887_52408131_minus_sgRNA_w5_a_meta.csv,VEP_OUTPUT_chr3_52407887_52408131_minus_sgRNA_w5_a_pam.tsv, by="id", all.x=TRUE)
E6_SGA_META_VEP<-merge(chr3_52407284_52407528_minus_sgRNA_w6_a_meta.csv,VEP_OUTPUT_chr3_52407284_52407528_minus_sgRNA_w6_a_pam.tsv, by="id", all.x=TRUE)
E7_SGA_META_VEP<-merge(chr3_52407123_52407367_minus_sgRNA_w7_a_meta.csv,VEP_OUTPUT_chr3_52407123_52407367_minus_sgRNA_w7_a_pam.tsv, by="id", all.x=TRUE)
E9_SGA_META_VEP<-merge(chr3_52406199_52406443_minus_sgRNA_w9_a_meta.csv,VEP_OUTPUT_chr3_52406199_52406443_minus_sgRNA_w9_a_pam.tsv, by="id", all.x=TRUE)
E10_SGA_META_VEP<-merge(chr3_52405725_52405969_minus_sgRNA_w10_a_meta.csv,VEP_OUTPUT_chr3_52405725_52405969_minus_sgRNA_w10_a_pam.tsv, by="id", all.x=TRUE)
E11_1_SGA_META_VEP<-merge(chr3_52405129_52405373_minus_sgRNA_w11_x_a_meta.csv,VEP_OUTPUT_chr3_52405129_52405373_minus_sgRNA_w11_x_a_pam.tsv, by="id", all.x=TRUE)
E11_2_SGA_META_VEP<-merge(chr3_52405033_52405277_minus_sgRNA_w11_x_a_meta.csv,VEP_OUTPUT_chr3_52405033_52405277_minus_sgRNA_w11_x_a_pam.tsv, by="id", all.x=TRUE)
E12_1_SGA_META_VEP<-merge(chr3_52404452_52404696_minus_sgRNA_w12_x_a_meta.csv,VEP_OUTPUT_chr3_52404452_52404696_minus_sgRNA_w12_x_a_pam.tsv, by="id", all.x=TRUE)
E12_2_SGA_META_VEP<-merge(chr3_52404348_52404592_minus_sgRNA_w12_x_a_meta.csv,VEP_OUTPUT_chr3_52404348_52404592_minus_sgRNA_w12_x_a_pam.tsv, by="id", all.x=TRUE)
E13_1_SGA_META_VEP<-merge(chr3_52403738_52403982_minus_sgRNA_w13_1_a_meta.csv,VEP_OUTPUT_chr3_52403738_52403982_minus_sgRNA_w13_1_a_pam.tsv, by="id", all.x=TRUE)
E13_2_SGA_META_VEP<-merge(chr3_52403556_52403800_minus_sgRNA_w13_2_a_meta.csv,VEP_OUTPUT_chr3_52403556_52403800_minus_sgRNA_w13_2_a_pam.tsv, by="id", all.x=TRUE)
E13_3_SGA_META_VEP<-merge(chr3_52403378_52403622_minus_sgRNA_w13_3_a_meta.csv,VEP_OUTPUT_chr3_52403378_52403622_minus_sgRNA_w13_3_a_pam.tsv, by="id", all.x=TRUE)
E14_SGA_META_VEP<-merge(chr3_52403097_52403341_minus_sgRNA_w14_a_meta.csv,VEP_OUTPUT_chr3_52403097_52403341_minus_sgRNA_w14_a_pam.tsv, by="id", all.x=TRUE)
E15_SGA_META_VEP<-merge(chr3_52402703_52402947_minus_sgRNA_w15_a_meta.csv,VEP_OUTPUT_chr3_52402703_52402947_minus_sgRNA_w15_a_pam.tsv, by="id", all.x=TRUE)
E16_SGA_META_VEP<-merge(chr3_52402511_52402755_minus_sgRNA_w16_a_meta.csv,VEP_OUTPUT_chr3_52402511_52402755_minus_sgRNA_w16_a_pam.tsv, by="id", all.x=TRUE)
E17_1_SGA_META_VEP<-merge(chr3_52402282_52402526_minus_sgRNA_w17_x_a_meta.csv,VEP_OUTPUT_chr3_52402282_52402526_minus_sgRNA_w17_x_a_pam.tsv, by="id", all.x=TRUE)
E17_2_SGA_META_VEP<-merge(chr3_52402173_52402418_minus_sgRNA_w17_x_a_meta.csv,VEP_OUTPUT_chr3_52402173_52402418_minus_sgRNA_w17_x_a_pam.tsv, by="id", all.x=TRUE)

#annotation formatting - TAKE VALIANT OUTPUT and VEP output derrived from PAM-specific VCFs - SGB SET (deselected libraries not present in final output)
E1_SGB_META_VEP<-merge(chr3_52409736_52409980_minus_sgRNA_w1_b_meta.csv,VEP_OUTPUT_chr3_52409736_52409980_minus_sgRNA_w1_b_pam.tsv, by="id", all.x=TRUE)
#E2_SGB_META_VEP<-merge(chr3_52409602_52409846_minus_sgRNA_w2_b_meta.csv,VEP_OUTPUT_chr3_52409602_52409846_minus_sgRNA_w2_b_pam.tsv, by="id", all.x=TRUE)
E4_SGB_META_VEP<-merge(chr3_52408417_52408661_minus_sgRNA_w4_b_meta.csv,VEP_OUTPUT_chr3_52408417_52408661_minus_sgRNA_w4_b_pam.tsv, by="id", all.x=TRUE)
E5_SGB_META_VEP<-merge(chr3_52407887_52408131_minus_sgRNA_w5_b_meta.csv,VEP_OUTPUT_chr3_52407887_52408131_minus_sgRNA_w5_b_pam.tsv, by="id", all.x=TRUE)
E6_SGB_META_VEP<-merge(chr3_52407284_52407528_minus_sgRNA_w6_b_meta.csv,VEP_OUTPUT_chr3_52407284_52407528_minus_sgRNA_w6_b_pam.tsv, by="id", all.x=TRUE)
E7_SGB_META_VEP<-merge(chr3_52407123_52407367_minus_sgRNA_w7_b_meta.csv,VEP_OUTPUT_chr3_52407123_52407367_minus_sgRNA_w7_b_pam.tsv, by="id", all.x=TRUE)
E8_SGB_META_VEP<-merge(chr3_52406736_52406980_minus_sgRNA_w8_b_meta.csv,VEP_OUTPUT_chr3_52406736_52406980_minus_sgRNA_w8_b_pam.tsv, by="id", all.x=TRUE)
E9_SGB_META_VEP<-merge(chr3_52406199_52406443_minus_sgRNA_w9_b_meta.csv,VEP_OUTPUT_chr3_52406199_52406443_minus_sgRNA_w9_b_pam.tsv, by="id", all.x=TRUE)
E10_SGB_META_VEP<-merge(chr3_52405725_52405969_minus_sgRNA_w10_b_meta.csv,VEP_OUTPUT_chr3_52405725_52405969_minus_sgRNA_w10_b_pam.tsv, by="id", all.x=TRUE)
E11_1_SGB_META_VEP<-merge(chr3_52405129_52405373_minus_sgRNA_w11_x_b_meta.csv,VEP_OUTPUT_chr3_52405129_52405373_minus_sgRNA_w11_x_b_pam.tsv, by="id", all.x=TRUE)
E11_2_SGB_META_VEP<-merge(chr3_52405033_52405277_minus_sgRNA_w11_x_b_meta.csv,VEP_OUTPUT_chr3_52405033_52405277_minus_sgRNA_w11_x_b_pam.tsv, by="id", all.x=TRUE)
E12_1_SGB_META_VEP<-merge(chr3_52404452_52404696_minus_sgRNA_w12_x_b_meta.csv,VEP_OUTPUT_chr3_52404452_52404696_minus_sgRNA_w12_x_b_pam.tsv, by="id", all.x=TRUE)
E12_2_SGB_META_VEP<-merge(chr3_52404348_52404592_minus_sgRNA_w12_x_b_meta.csv,VEP_OUTPUT_chr3_52404348_52404592_minus_sgRNA_w12_x_b_pam.tsv, by="id", all.x=TRUE)
#E13_1_SGB_META_VEP<-merge(chr3_52403738_52403982_minus_sgRNA_w13_1_b_meta.csv,VEP_OUTPUT_chr3_52403738_52403982_minus_sgRNA_w13_1_b_pam.tsv, by="id", all.x=TRUE)
E13_2_SGB_META_VEP<-merge(chr3_52403556_52403800_minus_sgRNA_w13_2_b_meta.csv,VEP_OUTPUT_chr3_52403556_52403800_minus_sgRNA_w13_2_b_pam.tsv, by="id", all.x=TRUE)
E13_3_SGB_META_VEP<-merge(chr3_52403378_52403622_minus_sgRNA_w13_3_b_meta.csv,VEP_OUTPUT_chr3_52403378_52403622_minus_sgRNA_w13_3_b_pam.tsv, by="id", all.x=TRUE)
E14_SGB_META_VEP<-merge(chr3_52403097_52403341_minus_sgRNA_w14_b_meta.csv,VEP_OUTPUT_chr3_52403097_52403341_minus_sgRNA_w14_b_pam.tsv, by="id", all.x=TRUE)
E15_SGB_META_VEP<-merge(chr3_52402703_52402947_minus_sgRNA_w15_b_meta.csv,VEP_OUTPUT_chr3_52402703_52402947_minus_sgRNA_w15_b_pam.tsv, by="id", all.x=TRUE)
E16_SGB_META_VEP<-merge(chr3_52402511_52402755_minus_sgRNA_w16_b_meta.csv,VEP_OUTPUT_chr3_52402511_52402755_minus_sgRNA_w16_b_pam.tsv, by="id", all.x=TRUE)
E17_1_SGB_META_VEP<-merge(chr3_52402282_52402526_minus_sgRNA_w17_x_b_meta.csv,VEP_OUTPUT_chr3_52402282_52402526_minus_sgRNA_w17_x_b_pam.tsv, by="id", all.x=TRUE)
E17_2_SGB_META_VEP<-merge(chr3_52402173_52402418_minus_sgRNA_w17_x_b_meta.csv,VEP_OUTPUT_chr3_52402173_52402418_minus_sgRNA_w17_x_b_pam.tsv, by="id", all.x=TRUE)

#FUNCTION to add an additional field which is ID and SG combined - to label a variant annotation as either from SGA or SGB source - which is needed subsequently 
#label annotations with unique id
label <- function(exon, sg){
  input<-paste0("E",exon,"_SG",sg,"_META_VEP")
  input<-get(input)
  input$SG<-sg
  input$exon_group<-exon
  input$required_annotation_id<-paste0(input$id,"_",input$SG)
  countdfname<-paste0("E",exon,"_SG",sg,"_META_VEP")
  assign(countdfname, input, envir = .GlobalEnv)
}

#Function to apply to listed exons and designated sgRNA_A
lapply (c(1,2,3,4,5,6,7,9,10,"11_1","11_2","12_1","12_2","13_1","13_2","13_3",14,15,16,"17_1","17_2"), function(xx) {label(exon=xx,sg="A")
})
#Function to apply to listed exons and designated sgRNA_B
lapply (c(1,4,5,6,7,8,9,10,"11_1","11_2","12_1","12_2","13_2","13_3",14,15,16,"17_1","17_2"), function(xx) {label(exon=xx,sg="B")
})


#take the inherited dataset (post combination of guides and combination of overlapping targetons)
#if the variant is observed once (ie is not from combined from datasets, or is only observed in lirbary a or library b) then the required source of annoration
#will be whichever library is the source of the data (either A or B)
annotation_processed_single<-targeton_processed %>% filter(variant_n_instances_between_libraries %in% "1"|
                                                             variant_n_instances_between_libraries %in% "1|1" |
                                                             variant_n_instances_between_libraries %in% "NA") %>%
                                              mutate(required_annotation_source=case_when(str_detect(source_libraries,"A")~ "A", TRUE~ "B"))

#if the variant is derived from two libraries (A and B) then take metadata from source B unless the arrayed pam status contains a B guide, then take A
#THIS WORKS FOR BAP1 BUT MAY NOT WORK FOR OTHER GENES / LIRBARY COMBINATIONS
annotation_processed_duplicate<-targeton_processed %>% filter(variant_n_instances_between_libraries%in% "2"|
                                                                variant_n_instances_between_libraries %in% "2|2"|
                                                                variant_n_instances_between_libraries %in% "2|NA") %>%
                                                mutate(required_annotation_source=case_when(str_detect(arrayed_PAM_status,"_b")~ "A", TRUE~ "B"))
#BIND together these two dataframes
req_anno<-do.call("rbind.fill", list(annotation_processed_single, annotation_processed_duplicate))
#CREATE a new required annotation field in the dataframe
req_anno$required_annotation_id<-paste0(req_anno$id,"_",req_anno$required_annotation_source)
req_anno$exon_specific_names<-NULL
req_anno$arrayed_mseq_a<-NULL
req_anno$arrayed_mseq_b<-NULL

#COMBINE ALL META AND VEP FRAMES THAT ARE NEEDED TO ANNOTATE THE DATA YOU HAVE
total_META_VEP <- do.call("rbind.fill", list(E1_SGA_META_VEP,
                                             E2_SGA_META_VEP,
                                             E3_SGA_META_VEP,
                                             E4_SGA_META_VEP,
                                             E5_SGA_META_VEP,
                                             E6_SGA_META_VEP,
                                             E7_SGA_META_VEP,
                                             E9_SGA_META_VEP,
                                             E10_SGA_META_VEP,
                                             E11_1_SGA_META_VEP,
                                             E11_2_SGA_META_VEP,
                                             E12_1_SGA_META_VEP,
                                             E12_2_SGA_META_VEP,
                                             E13_1_SGA_META_VEP,
                                             E13_2_SGA_META_VEP,
                                             E13_3_SGA_META_VEP,
                                             E14_SGA_META_VEP,
                                             E15_SGA_META_VEP,
                                             E16_SGA_META_VEP,
                                             E17_1_SGA_META_VEP,
                                             E17_2_SGA_META_VEP,
                                             E1_SGB_META_VEP,
                                             #E2_SGB_META_VEP,
                                             E4_SGB_META_VEP,
                                             E5_SGB_META_VEP,
                                             E6_SGB_META_VEP,
                                             E7_SGB_META_VEP,
                                             E8_SGB_META_VEP,
                                             E9_SGB_META_VEP,
                                             E10_SGB_META_VEP,
                                             E11_1_SGB_META_VEP,
                                             E11_2_SGB_META_VEP,
                                             E12_1_SGB_META_VEP,
                                             E12_2_SGB_META_VEP,
                                             #E13_1_SGB_META_VEP,
                                             E13_2_SGB_META_VEP,
                                             E13_3_SGB_META_VEP,
                                             E14_SGB_META_VEP,
                                             E15_SGB_META_VEP,
                                             E16_SGB_META_VEP,
                                             E17_1_SGB_META_VEP,
                                             E17_2_SGB_META_VEP))

#remove any variant annotations that are variants in the constant regions 
total_META_VEP<- total_META_VEP %>% filter(vcf_var_in_const==0 | is.na(vcf_var_in_const))
total_META_VEP$pam_mut_sgrna_id<-NULL  
total_META_VEP$exon_group<-NULL  
#create minimal fields needed to link mseq to all (mseq redundant ids)
minimal_annotation<-total_META_VEP[ ,c("required_annotation_id","mseq")]
#REQ_ANNO has one id per unique mseq - need to get all ids per unique mseq
#1.fish for the required mseq using the required_annotation_id
annotation_key_association<-merge(minimal_annotation, req_anno, by="required_annotation_id", all.x=FALSE)
annotation_key_association$required_annotation_id<-NULL
annotation_key_association$id<-NULL
#2.get all the ids associated with the appropriate mseqs by merging the data by mseq
annotation_full<-merge(total_META_VEP, annotation_key_association,  by="mseq", all.x=FALSE)
#check which ids are missing from the final dataset versus the metadata
check<-subset(total_unique_id_count, !(id %in% filtered_annotation_full$id))
#all of the variants that are missing are likely to be pam site edits that further modify the mseq to be something unique, that would not occur in wild-type mutation 
#most will also be duplciated due to merging with expanded metadata annotation
duplicated <- annotation_full  %>% 
  group_by(id) %>%
  filter(n() > 1)
####DUPLICATE ACCESSION PROCESSING - START 
####DUPLICATE ACCESSION PROCESSING
####DUPLICATE ACCESSION PROCESSING
####DUPLICATE ACCESSION PROCESSING
####DUPLICATE ACCESSION PROCESSING
####DUPLICATE ACCESSION PROCESSING
####DUPLICATE ACCESSION PROCESSING

#There are some mseqs associated with an ID that we don't want because they are an meseq that is a duplicate of a mseq that is annotated as having a pam edit, but isn't itself annotated as a pam containing mseq
#this happens with 1del deletions, and in some cases 2del if the pam site is in a 2 del location (intron as is the case for BAP1 exon1 library b)
#some ids are represented by multiple mseqs 
#these mseqs will have remained distinct in the annotation process
#for library a and library b the required annotation id will have pulled the correct mseq(s)
#some required annotation ids will have pulled mseqs that are duplications

##remove accessions that are duplicated for id where the annotation states that is not observed twice in the dataset - this is an artifact of merging mseqs and then pulling by mseqs again

##can delete these - start 
#annotation_full$variant_n_instances_between_libraries<-as.character(annotation_full$variant_n_instances_between_libraries)
#annotation_full$variant_n_instances_between_targetons<-as.character(annotation_full$variant_n_instances_between_targetons)

#filtered_annotation_full <- annotation_full[!(duplicated(annotation_full$id) & (annotation_full$variant_n_instances_between_libraries == 1)), ]
#filtered_annotation_full <- filtered_annotation_full[!(duplicated(filtered_annotation_full$id) & is.na(filtered_annotation_full$variant_n_instances_between_libraries)), ]
##can delete these - end

#create a column that counts the number of duplicates (1=not duplciated)
filtered_annotation_full  <-annotation_full %>%
  group_by(id) %>% mutate(duplicated_id= n())
#remove accessions where the duplication is not due to targeton overlap + where states has been seen only once or is NA, both of these are artifacts of merging and not true
duplicated_set<-filtered_annotation_full %>% filter(variant_n_instances_between_targetons==1 & duplicated_id >1) %>% filter(!variant_n_instances_between_libraries==1) %>% filter(!variant_n_instances_between_libraries=="NA")
#check again for which acccessions are still duplicated
duplicated_2  <-duplicated_set %>%
  group_by(id) %>% mutate(duplicated_id_2= n())
#annotated this new dataframe with a column to state if the mseq contains pam site annotation, duplicates which don't contain this annotation are an artifact of merging
duplicated_set_2<- duplicated_2  %>% mutate(is_guide=case_when(str_detect(arrayed_PAM_status,"sgRNA") ~ "Y", TRUE ~ "N"))
#annotate the datafram with a removal indication column. Remove the row if it is duplicate for id and contains no guide annotation and is not exon17 (special case) or 
#is exon17 duplication without a pam status that indicates the row is from an array of values - this is a special case becuase the GG at boarder of overlap for 17_1 and 17_2, one G is in overlapping region and one is not, which leads to the same mseq being annotated with multiple ids
duplicated_set_2<- duplicated_set_2  %>% mutate(remove=case_when(is_guide=="N" & duplicated_id_2==2 & !str_detect(source_targetons,"17") | str_detect(source_targetons,"17") & arrayed_PAM_status=="," ~ "delete", TRUE ~ "keep"))
#remove rows based on the removal indication column
duplicated_set_3<-duplicated_set_2 %>% filter(!remove %in% 'delete')
#remove the columns that were generated to process the filtering  
duplicated_set_3$duplicated_id<-NULL
duplicated_set_3$duplicated_id_2<-NULL
duplicated_set_3$remove<-NULL
duplicated_set_3$is_guide<-NULL
#MAKE UNIQUE and count the number of unique ids in the final list, should be the same as 'annotation_full' > 'variant_n_instances_between_targetons==1' > make unique
single_targeton_duplication_unique<-duplicated_set_3[!duplicated(duplicated_set_3$id), ]
#now attend to duplications in annotation due to overlapping targetons being pulled by mseq from the anotation table
filtered_annotation_full  <-annotation_full %>%
  group_by(id) %>% mutate(duplicated_id= n())
#get duplicated targetons that are duplicated due to overlapping regions
duplicated_set_targeton_set <- filtered_annotation_full %>% filter(variant_n_instances_between_targetons >1 & duplicated_id >1)
#remove x4 duplciates that do not have pam status annotation where they should
duplicated_set_targeton_set_2 <- duplicated_set_targeton_set %>% filter(duplicated_id >1) %>% mutate(remove=case_when(duplicated_id==4 & arrayed_PAM_status=="NA|NA" | duplicated_id==4 & arrayed_PAM_status==",|," ~ "delete", TRUE ~ "keep"))
#delete those rows that do not contain the full annotation of pam status - an artifact of merging by mseq
duplicated_set_targeton_set_3<-duplicated_set_targeton_set_2 %>% filter(!remove %in% 'delete')
#remove the columns that were generated to process the filtering  
duplicated_set_targeton_set_3$duplicated_id<-NULL
duplicated_set_targeton_set_3$remove<-NULL
#MAKE UNIQUE and count the number of unique ids in the final list, should be the same as 'annotation_full' > 'variant_n_instances_between_targetons==2' > make unique
double_targeton_duplication_unique<-duplicated_set_targeton_set_3[!duplicated(duplicated_set_targeton_set_3$id), ]

###creation of final dataset - post duplicate removal

#Rbind to get a list of the correct accessions, that were previously duplicated but now are unique and contain the appropriate values
filtered_duplicates <- do.call("rbind.fill", list(single_targeton_duplication_unique, double_targeton_duplication_unique))
#remove those accessions that are duplicated from the main dataframe
annotation_full_minus_duplicates<-annotation_full[!annotation_full$id %in% filtered_duplicates$id,]
#add back in the filtered previously duplicated accesssions
final_annotated_dataset<-do.call("rbind.fill", list(annotation_full_minus_duplicates, filtered_duplicates))

#check that no duplicate ids remain
duplicated <- final_annotated_dataset  %>% 
  group_by(id) %>%
  filter(n() > 1)
#check which variants are not in the final dataset which were in the orignal post DEseq2 list
check<-subset(total_unique_id_count, !(id %in% final_annotated_dataset$id))
#are these accessions unique or are there duplicates - shouldn't be an any duplciates
check_unique<-check[!duplicated(check$id), ]

#this is your final annotated dataset: 'final_annotated_dataset', write it to file
write.csv(final_annotated_dataset,"./final_annotated_dataset.csv", row.names = FALSE)

####DUPLICATE ACCESSION PROCESSING - END
####DUPLICATE ACCESSION PROCESSING
####DUPLICATE ACCESSION PROCESSING
####DUPLICATE ACCESSION PROCESSING
####DUPLICATE ACCESSION PROCESSING
####DUPLICATE ACCESSION PROCESSING
####DUPLICATE ACCESSION PROCESSING

#CLEANUP THE DATAFRAME ADD USEFUL COLUMNS - SECTION 1
#CLEANUP THE DATAFRAME ADD USEFUL COLUMNS
#CLEANUP THE DATAFRAME ADD USEFUL COLUMNS
#CLEANUP THE DATAFRAME ADD USEFUL COLUMNS
#CLEANUP THE DATAFRAME ADD USEFUL COLUMNS

#simplify the pam status field and create a flag field for variants that only come from one dataset and are pam variant codon
clean_final_annotated_dataset<- final_annotated_dataset %>% 
  mutate(pam_codon=case_when(str_detect(arrayed_PAM_status,"sgRNA") ~ "Y", TRUE ~ "N")) %>%
  mutate(pam_codon=case_when(str_detect(single_PAM_status,"sgRNA") ~ "Y", TRUE ~ pam_codon)) %>%
  mutate(pam_flag=case_when(source_dataset=="not_combined" & str_detect(arrayed_PAM_status,"sgRNA") ~ "Y", TRUE ~ "N")) %>%
mutate(pam_flag=case_when(source_dataset=="not_combined" & str_detect(single_PAM_status,"sgRNA") ~ "Y", TRUE ~ pam_flag))

#some HGVSc fields are missing in pam_flag codons (b/c non_ref allele as starting nucleotide)- these should be added but will likely be disregared in analyses becasue of pam_flag status
clean_final_annotated_dataset<- clean_final_annotated_dataset %>% 
mutate(HGVSc=case_when(!HGVSc=="-" ~ HGVSc, TRUE ~ paste0(Feature,":c.",CDS_position,ref,">",new)))

#write out the full dataset
write.csv(clean_final_annotated_dataset,"./clean_final_annotated_dataset.csv", row.names = FALSE)

  
#METRICS SECTION - START
#METRICS SECTION 
#METRICS SECTION 
#METRICS SECTION 
#METRICS SECTION
#METRICS SECTION
#METRICS SECTION 
#METRICS SECTION 
#METRICS SECTION 
#METRICS SECTION 

#make a column and fill with NA
x_final_frame_class<-clean_final_annotated_dataset
x_final_frame_class$class <-NA
#create the categories
x_final_frame_class$class[x_final_frame_class$processed_BH_FDR_D4_D21<0.01 & x_final_frame_class$processed_Z_D4_D21<0] <- 'depleted'
x_final_frame_class$class[x_final_frame_class$processed_BH_FDR_D4_D21<0.01 & x_final_frame_class$processed_Z_D4_D21>0] <- 'enriched'
x_final_frame_class$class[x_final_frame_class$processed_BH_FDR_D4_D21>0.01] <- 'unchanged'
#check the categories
table(x_final_frame_class$class)
#check there are no NA rows in new class column
is_na_class<-x_final_frame_class %>% filter(class %in% NA)

x_final_frame_class <- x_final_frame_class %>% relocate(class, .after = id)

#write the df to file
write.csv(x_final_frame_class,"./clean_bap1_v2.csv", row.names = FALSE)


#create columns to give an idea of library_a vs library_b concordant where possible - useful perhaps for plots and for assessing rare population variants where strong evidence needed beyond simmple metric
#take all variants that have reads for
library_concordance_frame<-do.call("rbind.fill", list(total_combined_bap1_a,total_combined_bap1_b,total_single_bap1))
#select the median scaled lfc and p-value 
library_concordance_frame<-library_concordance_frame %>% select(id,adj_lfc_D4_D21,uncombined_two_tailed_p_D4_D21)
#collapse based on the id which will be shared between libraries and targetons
library_concordance_collapse <- library_concordance_frame %>%
  group_by(id) %>%
  dplyr::summarize(adj_lfc_D4_D21 = paste(adj_lfc_D4_D21, collapse = ','), 
                   uncombined_two_tailed_p_D4_D21 = paste(uncombined_two_tailed_p_D4_D21, collapse = ',')) %>%
  ungroup() %>%
  #ungroup the key data into primary and secondary, tertiary and quartenary targetons fields - most will be 2 observations, some will be 3 or 4 if targeton overlapping reion
  separate(adj_lfc_D4_D21, sep=",", c("adj_lfc_D4_D21_obvs_1","adj_lfc_D4_D21_obvs_2","adj_lfc_D4_D21_obvs_3","adj_lfc_D4_D21_obvs_4")) %>%
  separate(uncombined_two_tailed_p_D4_D21, sep=",", c("uncombined_two_tailed_p_D4_D21_obvs_1","uncombined_two_tailed_p_D4_D21_obvs_2","uncombined_two_tailed_p_D4_D21_obvs_3","uncombined_two_tailed_p_D4_D21_obvs_4"))
#create new frame for processing
library_concordance_calculated<-library_concordance_collapse
#make numeric for columns
library_concordance_calculated[, c(2:9)] <- sapply(library_concordance_calculated[, c(2:9)], as.numeric)

#create new concordance column - fill and perform calulations based on unprocessed lfc and p-value
library_concordance_calculated$concordance<-NA
library_concordance_calculated$concordance[library_concordance_calculated$adj_lfc_D4_D21_obvs_1<0 & 
                                             library_concordance_calculated$adj_lfc_D4_D21_obvs_2<0 &
                                             library_concordance_calculated$uncombined_two_tailed_p_D4_D21_obvs_1<0.01 & 
                                             library_concordance_calculated$uncombined_two_tailed_p_D4_D21_obvs_2<0.01] <- 'concordant_depleted'

library_concordance_calculated$concordance[library_concordance_calculated$adj_lfc_D4_D21_obvs_1>0 & 
                                             library_concordance_calculated$adj_lfc_D4_D21_obvs_2>0 &
                                             library_concordance_calculated$uncombined_two_tailed_p_D4_D21_obvs_1<0.01 & 
                                             library_concordance_calculated$uncombined_two_tailed_p_D4_D21_obvs_2<0.01] <- 'concordant_enriched'

library_concordance_calculated$concordance[library_concordance_calculated$uncombined_two_tailed_p_D4_D21_obvs_1>0.01 & 
                                             library_concordance_calculated$uncombined_two_tailed_p_D4_D21_obvs_2>0.01] <- 'concordant_unchanged'

library_concordance_calculated$concordance[is.na(library_concordance_calculated$concordance)] <- 'non_concordant_or_single_obvs'

#check numbers in table
table(library_concordance_calculated$concordance)
#merge with the final dataset to add columns
concordance_x_final_frame_class<-merge(x_final_frame_class,library_concordance_calculated, by="id", all.x=TRUE)
#check numbers in table
table(concordance_x_final_frame_class$concordance)

#add a simplified PAM STATUS field useful for concordnace plots
test<-concordance_x_final_frame_class 

  test<-concordance_x_final_frame_class %>% mutate(PAM_status_simplified=case_when(pam_codon=="Y" & str_detect(arrayed_PAM_status,"_a") ~ "A", TRUE ~ "-")) %>%
    mutate(PAM_status_simplified=case_when(pam_codon=="Y" & str_detect(arrayed_PAM_status,"_b") ~ "B", TRUE ~ PAM_status_simplified))
  

#SANITY CHECKING FRAME BASED MUTATORS 
#730 codons inc start and stop
#10 are split codons between exons # codons: 13,23,41,146,194,220,311,417,577,686 are all split codons
#stop has 1 stops missing this is becuase "ENST00000460680.6.ENSG00000163930.10_chr3:52407250_52407252_GAA>TCA_stop_rc" was filtered out for low counts
#ala has 21 missing but there are probably 21 wildtype ala GCC = THIS IS TRUE

#ROUGH PLOT TO CHECK CONCORDANCE AND CORRELATION OF LFC
#opens a new PDF in the working directory and pastes each plot comparison - only shows those targetons that have a posibility of duplication, some guides will only be present in library a or library b - filered by variant duplicatio field
pdf("sga_sgb_lfc_guide_comparison.pdf")
test %>% filter(variant_n_instances_between_libraries=="2"|variant_n_instances_between_libraries=="2|2"|variant_n_instances_between_libraries=="2|NA") %>% ggscatter(x="adj_lfc_D4_D21_obvs_1",y="adj_lfc_D4_D21_obvs_2", color = "PAM_status_simplified", alpha= 0.35, size=2,  add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
test %>% filter(variant_n_instances_between_libraries=="2"|variant_n_instances_between_libraries=="2|2"|variant_n_instances_between_libraries=="2|NA") %>% ggscatter(x="adj_lfc_D4_D21_obvs_1",y="processed_LFC_D4_D21", color = "PAM_status_simplified", alpha= 0.35, size=2,  add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
test %>% filter(variant_n_instances_between_libraries=="2"|variant_n_instances_between_libraries=="2|2"|variant_n_instances_between_libraries=="2|NA") %>% ggscatter(x="processed_LFC_D4_D21",y="adj_lfc_D4_D21_obvs_2", color = "PAM_status_simplified", alpha= 0.35, size=2,   add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))

test %>% filter(variant_n_instances_between_libraries=="2"|variant_n_instances_between_libraries=="2|2"|variant_n_instances_between_libraries=="2|NA") %>% ggscatter(x="adj_lfc_D4_D21_obvs_1",y="adj_lfc_D4_D21_obvs_2", color = "concordance", shape = "PAM_status_simplified", alpha= 0.35, size=2,  add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
test %>% filter(variant_n_instances_between_libraries=="2"|variant_n_instances_between_libraries=="2|2"|variant_n_instances_between_libraries=="2|NA") %>% ggscatter(x="adj_lfc_D4_D21_obvs_1",y="processed_LFC_D4_D21", color = "concordance", shape = "PAM_status_simplified", alpha= 0.35, size=2,  add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
test %>% filter(variant_n_instances_between_libraries=="2"|variant_n_instances_between_libraries=="2|2"|variant_n_instances_between_libraries=="2|NA") %>% ggscatter(x="processed_LFC_D4_D21",y="adj_lfc_D4_D21_obvs_2", color = "concordance", shape = "PAM_status_simplified", alpha= 0.35, size=2,   add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n"))
dev.off()



#METRICS SECTION - END
#METRICS SECTION 
#METRICS SECTION 
#METRICS SECTION 
#METRICS SECTION
#METRICS SECTION
#METRICS SECTION 
#METRICS SECTION 
#METRICS SECTION 
#METRICS SECTION 


#CLEANUP THE DATAFRAME ADD USEFUL COLUMNS - SECTION 2
#CLEANUP THE DATAFRAME ADD USEFUL COLUMNS
#CLEANUP THE DATAFRAME ADD USEFUL COLUMNS
#CLEANUP THE DATAFRAME ADD USEFUL COLUMNS
#CLEANUP THE DATAFRAME ADD USEFUL COLUMNS

#select key columns

#add a pos_ref_new field 
test$chrom_pos_ref_alt <- paste0(test$ref_chr,"_",test$mut_position,"_",test$ref,"_",test$new)
test_colnames<-as.data.frame(colnames(test))

write.csv(test_colnames,"./test_colnames.csv", row.names = FALSE)

full_fields<-test

#reorder the dataframe
clean_final_annotated_dataset_reordered <-full_fields[, c("id",
                                                          "HGVSc",
                                                          "HGVSp",
                                                          "chrom_pos_ref_alt",
                                                          "class",
                                                          "mutator",
                                                          "processed_LFC_D4_D7",
                                                          "processed_LFC_D4_D10",
                                                          "processed_LFC_D4_D14",
                                                          "processed_LFC_D4_D21",
                                                          "processed_LFC_continuous",
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
                                                          "two_tailed_p_D4_D7",
                                                          "two_tailed_p_D4_D10",
                                                          "two_tailed_p_D4_D14",
                                                          "two_tailed_p_D4_D21",
                                                          "two_tailed_p_continuous",
                                                          "SE_bind_D4_D7",
                                                          "SE_bind_D4_D10",
                                                          "SE_bind_D4_D14",
                                                          "SE_bind_D4_D21",
                                                          "SE_bind_continuous",
                                                          "concordance",
                                                          "pam_codon",
                                                          "pam_flag",
                                                          "PAM_status_simplified",
                                                          "variant_n_instances_between_libraries",
                                                          "variant_n_instances_between_targetons",
                                                          "source_dataset",
                                                          "source_targetons",
                                                          "source_libraries",
                                                          "Consequence",
                                                          "cDNA_position",
                                                          "CDS_position",
                                                          "Protein_position",
                                                          "Amino_acids",
                                                          "Codons",
                                                          "mut_position",
                                                          "ref",
                                                          "new",
                                                          "ref_aa",
                                                          "alt_aa",
                                                          "mut_type",
                                                          "adj_lfc_D4_D21_obvs_1",
                                                          "adj_lfc_D4_D21_obvs_2",
                                                          "adj_lfc_D4_D21_obvs_3",
                                                          "adj_lfc_D4_D21_obvs_4",
                                                          "uncombined_two_tailed_p_D4_D21_obvs_1",
                                                          "uncombined_two_tailed_p_D4_D21_obvs_2",
                                                          "uncombined_two_tailed_p_D4_D21_obvs_3",
                                                          "uncombined_two_tailed_p_D4_D21_obvs_4",
                                                          "mseq",
                                                          "species",
                                                          "assembly",
                                                          "gene_id",
                                                          "transcript_id",
                                                          "src_type",
                                                          "ref_chr",
                                                          "ref_strand",
                                                          "ref_start",
                                                          "ref_end",
                                                          "revc",
                                                          "ref_seq",
                                                          "pam_seq",
                                                          "vcf_alias",
                                                          "vcf_var_id",
                                                          "oligo_length",
                                                          "mseq_no_adapt",
                                                          "pam_mut_annot",
                                                          "mave_nt",
                                                          "vcf_var_in_const",
                                                          "Location",
                                                          "Allele",
                                                          "Gene",
                                                          "Feature",
                                                          "Feature_type",
                                                          "Existing_variation",
                                                          "IMPACT",
                                                          "DISTANCE",
                                                          "STRAND",
                                                          "FLAGS",
                                                          "PICK",
                                                          "SYMBOL",
                                                          "SYMBOL_SOURCE",
                                                          "HGNC_ID",
                                                          "BIOTYPE",
                                                          "CANONICAL",
                                                          "MANE_SELECT",
                                                          "MANE_PLUS_CLINICAL",
                                                          "TSL",
                                                          "APPRIS",
                                                          "CCDS",
                                                          "ENSP",
                                                          "SWISSPROT",
                                                          "TREMBL",
                                                          "UNIPARC",
                                                          "UNIPROT_ISOFORM",
                                                          "SOURCE",
                                                          "SIFT",
                                                          "PolyPhen",
                                                          "EXON",
                                                          "INTRON",
                                                          "DOMAINS",
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
                                                          "CLIN_SIG",
                                                          "SOMATIC",
                                                          "PHENO",
                                                          "PUBMED",
                                                          "MOTIF_NAME",
                                                          "MOTIF_POS",
                                                          "HIGH_INF_POS",
                                                          "MOTIF_SCORE_CHANGE",
                                                          "TRANSCRIPTION_FACTORS",
                                                          "gnomAD",
                                                          "gnomAD_FLAG",
                                                          "gnomAD_AF.1",
                                                          "ClinVar",
                                                          "ClinVar_CLNSIG",
                                                          "ClinVar_CLNREVSTAT",
                                                          "dbSNP",
                                                          "COSMIC",
                                                          "COSMIC_CNT",
                                                          "SG",
                                                          "required_annotation_id",
                                                          "arrayed_PAM_status",
                                                          "single_PAM_status",
                                                          "primary_tg_lfcSE_D4_D7",
                                                          "secondary_tg_lfcSE_D4_D7",
                                                          "primary_tg_lfcSE_D4_D10",
                                                          "secondary_tg_lfcSE_D4_D10",
                                                          "primary_tg_lfcSE_D4_D14",
                                                          "secondary_tg_lfcSE_D4_D14",
                                                          "primary_tg_lfcSE_D4_D21",
                                                          "secondary_tg_lfcSE_D4_D21",
                                                          "primary_tg_lfcSE_continuous",
                                                          "secondary_tg_lfcSE_continuous",
                                                          "primary_tg_adj_lfc_D4_D7",
                                                          "secondary_tg_adj_lfc_D4_D7",
                                                          "primary_tg_adj_lfc_D4_D10",
                                                          "secondary_tg_adj_lfc_D4_D10",
                                                          "primary_tg_adj_lfc_D4_D14",
                                                          "secondary_tg_adj_lfc_D4_D14",
                                                          "primary_tg_adj_lfc_D4_D21",
                                                          "secondary_tg_adj_lfc_D4_D21",
                                                          "primary_tg_adj_lfc_continuous",
                                                          "secondary_tg_adj_lfc_continuous",
                                                          "weight_tg1_D4_D7",
                                                          "weight_tg1_D4_D10",
                                                          "weight_tg1_D4_D14",
                                                          "weight_tg1_D4_D21",
                                                          "weight_tg1_continuous",
                                                          "weight_tg2_D4_D7",
                                                          "weight_tg2_D4_D10",
                                                          "weight_tg2_D4_D14",
                                                          "weight_tg2_D4_D21",
                                                          "weight_tg2_continuous",
                                                          "sum_of_weight_D4_D7",
                                                          "sum_of_weight_D4_D10",
                                                          "sum_of_weight_D4_D14",
                                                          "sum_of_weight_D4_D21",
                                                          "sum_of_weight_continuous",
                                                          "weighted_primary_tg_LFC_D4_D7",
                                                          "weighted_secondary_tg_LFC_D4_D7",
                                                          "weighted_primary_tg_LFC_D4_D10",
                                                          "weighted_secondary_tg_LFC_D4_D10",
                                                          "weighted_primary_tg_LFC_D4_D14",
                                                          "weighted_secondary_tg_LFC_D4_D14",
                                                          "weighted_primary_tg_LFC_D4_D21",
                                                          "weighted_secondary_tg_LFC_D4_D21",
                                                          "weighted_primary_tg_LFC_continuous",
                                                          "weighted_secondary_tg_LFC_continuous",
                                                          "sum_of_weighted_LFC_D4_D7",
                                                          "sum_of_weighted_LFC_D4_D10",
                                                          "sum_of_weighted_LFC_D4_D14",
                                                          "sum_of_weighted_LFC_D4_D21",
                                                          "sum_of_weighted_LFC_continuous",
                                                          "required_annotation_source")]
#write the df to file
write.csv(clean_final_annotated_dataset_reordered,"./clean_bap1.csv", row.names = FALSE)

################# FINAL ANNOTATION ############# END
################# FINAL ANNOTATION ############# END
################# FINAL ANNOTATION ############# END
################# FINAL ANNOTATION ############# END
################# FINAL ANNOTATION ############# END 
################# FINAL ANNOTATION ############# END


#Further refinement of columns and data cleaning for plots #### START
#Further refinement of columns and data cleaning for plots #### START
#Further refinement of columns and data cleaning for plots #### START
#Further refinement of columns and data cleaning for plots #### START
#Further refinement of columns and data cleaning for plots #### START
   
#read in the dataframe
bap1_sge_v2 <- read.csv("/Users/aw28/Desktop/BAP1_analysis/mseq_run_final/clean_bap1.csv")

#METRICS AND THRESHOLDS - START 
bap1_sge_v3<-bap1_sge_v2 %>% mutate(D4_D21_cross_threshold = case_when(processed_Z_D4_D21 < -4.450931 ~ "T", TRUE ~ "F")) %>%
  mutate(D4_D21_specific_classified = case_when(processed_Z_D4_D21 < -4.450931 & processed_BH_FDR_D4_D21 <0.01 ~ "depleted", TRUE ~ "not_depleted")) %>%
  mutate(D4_D21_general_classified = class)
bap1_sge_v3$class<-NULL

#Consequence simplification
bap1_sge_v4<- bap1_sge_v3  %>% mutate(slim_consequence=case_when(str_detect(Consequence,"stop_gained") ~ "stop_gained")) %>%
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

#check categories
count<-as.data.frame(table(bap1_sge_v4$slim_consequence))
#write df out
write.csv(bap1_sge_v4,"./bap1_sge_v4.csv", row.names = FALSE)


#targton simplification
bap1_sge_v5 <- bap1_sge_v4 %>% mutate(region=case_when(source_targetons=="1" ~ "1")) %>%
  mutate(region=case_when(source_targetons=="2" ~ "2", TRUE ~ region)) %>%
  mutate(region=case_when(source_targetons=="3" ~ "3", TRUE ~ region)) %>% 
  mutate(region=case_when(source_targetons=="4" ~ "4", TRUE ~ region)) %>% 
  mutate(region=case_when(source_targetons=="5" ~ "5", TRUE ~ region)) %>%
  mutate(region=case_when(source_targetons=="6" ~ "6", TRUE ~ region)) %>%
  mutate(region=case_when(source_targetons=="7" ~ "7", TRUE ~ region)) %>% 
  mutate(region=case_when(source_targetons=="8" ~ "8", TRUE ~ region)) %>% 
  mutate(region=case_when(source_targetons=="9" ~ "9", TRUE ~ region)) %>% 
  mutate(region=case_when(source_targetons=="10" ~ "10", TRUE ~ region)) %>% 
  mutate(region=case_when(str_detect(source_targetons,"11_") ~ "11", TRUE ~ region)) %>% 
  mutate(region=case_when(str_detect(source_targetons,"12_") ~ "12", TRUE ~ region)) %>% 
  mutate(region=case_when(str_detect(source_targetons,"13_") ~ "13", TRUE ~ region)) %>% 
  mutate(region=case_when(source_targetons=="14" ~ "14", TRUE ~ region)) %>% 
  mutate(region=case_when(source_targetons=="15" ~ "15", TRUE ~ region)) %>% 
  mutate(region=case_when(source_targetons=="16" ~ "16", TRUE ~ region)) %>% 
  mutate(region=case_when(str_detect(source_targetons,"17_") ~ "17", TRUE ~ region))

region_test<-as.data.frame(table(bap1_sge_v5$region))


write.csv(bap1_sge_v5,"./bap1_sge_v5.csv", row.names = FALSE)

#####REORDER THE DATAFRAME
bap1_sge_v5_reordered <- bap1_sge_v5 %>% relocate(D4_D21_general_classified, .after = mutator) %>%
  relocate(D4_D21_specific_classified, .after = D4_D21_general_classified) %>% 
  relocate(slim_consequence, .after = Consequence) %>% 
  relocate(region, .before = source_dataset) 
bap1_sge_v5<-bap1_sge_v5_reordered

#factor order to make regions appear in plots numerically increasing
bap1_sge_v5$region<-as.numeric(bap1_sge_v5$region)
bap1_sge_v5<-bap1_sge_v5[order(bap1_sge_v5$region),]
bap1_sge_v5$region<-as.factor(bap1_sge_v5$region)


#some plots require unique HGVSc
bap1_sge_v5_unique<-bap1_sge_v5[!duplicated(bap1_sge_v5$HGVSc),]

#density plot - START 
density_frame<- bap1_sge_v5_unique %>% filter(slim_consequence %in% "codon_deletion" |
                                                slim_consequence %in% "frameshift" |
                                                slim_consequence %in% "intron" |
                                                slim_consequence %in% "missense" |
                                                slim_consequence %in% "stop_gained" |
                                                slim_consequence %in% "splice_acceptor" |
                                                slim_consequence %in% "splice_donor" |
                                                slim_consequence %in% "synonymous" |
                                                slim_consequence %in% "UTR")

density_frame$slim_consequence<-as.factor(density_frame$slim_consequence)
density_frame$slim_consequence<- factor(density_frame$slim_consequence, 
                                        levels=c("synonymous", "stop_gained", "missense", "codon_deletion", "frameshift", "intron", "splice_donor", "splice_acceptor", "UTR"))
density_frame$slim_consequence<- fct_rev(density_frame$slim_consequence)



df <- bap1_sge_v5_unique_density_frame[ -c(7,56:76,80,83:88,91:97,99:103,109,110,131:135,145:200) ]

#stops excel converting to dates
df$INTRON<-paste0("'",df$INTRON,"'")
df$EXON<-paste0("'",df$EXON,"'")

write.csv(df,"./waters_bap1_sge_dataset.csv", row.names = FALSE)

