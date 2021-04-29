#Whilelist barcodes and to plot Cummulative distribution plots
#===========================================================<>
  #Whitelisting
    #criteria: Mean => Drop barcodes with frequency '<' mean
  #Cummulative Distribution
    #ecdf plots
##Check bash qsub sripts for:
    #@Barcode and UMI extraction 
    #@barcode-based read frequency 
    #@Salmon indexing and Alevin alignment 
#===========================================================<>

gc() #for memory management
rm(list=ls()) #Clear environment
setwd('/projectnb/bf528/users/lava_lamp/project_4/esaake_pr4/p4/')

install.packages("BiocManager")
library("dplyr")
library("ggplot2")
library("ggpubr")
library("gplots")
library("tidyverse")
library("samtools")

#Metadata
metapath<-"/projectnb/bf528/users/lava_lamp/project_4/esaake_pr4/p4/p4_sample_metadata.txt"
dirp<"/projectnb/bf528/users/lava_lamp/project_4/esaake_pr4/p4"
metadata<-read_csv(metapath)
spec(metadata)#Metadata fields and type

#Sample codes are under column "Run"
SRR51years <- subset(metadata, metadata$AGE==51 & metadata$sex=="female", select=c("Run","AGE","sex"))
#write_csv(SRR51years, "/projectnb/bf528/users/lava_lamp/project_4/esaake_pr4/p4/shortlistedsamples51.csv")



#================================================>
#==================STEPS========================<:
#------------------------------------------------>
  #1. Sequence QC
    #*GC content/Read quality/phred quality score
      #*apply fastqc or multiqc
  #2. Alignment + QC
    #*After assessing fastq files to be of high quality then you do alignment
      #* RNA-> use STAR
      #* NON-RNA-> bwa/bowtie
    #QC:
      #*RSeQC or multiqc
  #3. Quantification
      #*STAR + htseq-count
      #*Salmon ***
  #4. Count Matrix Normalization
      #*Within cell normalization (divide column by column total)
      #*within dataset(divide by total number of reads)
      #*All methods from bulk also apply
        #**CPM/FPKM/DESeq2 etx

#*There are 13 samples Use only the SRR files associated with the 51 year old female donor for further analysis.
#*
#**************************************************************************
#==========================================================================
#=>=>=>Barcode extraction and processing -> awk (qsub file)
#=>=>=>Barcode frequency count------------> awk (qsun file)
#=>=>=>Whitelisting-----------------------> This R
#=>=>=>Indexing and aligning -------------> grep/awk/sed and Salmon Alevin
#==========================================================================



###BARCODE WHITELISTING
##======================================

#Reading barcodes and their frequecy counts
SRR3879604_counts<- read.csv('/projectnb/bf528/users/lava_lamp/project_4/esaake_pr4/p4/results/count_SRR3879604.csv', col.names=c('barcode', 'count'), sep='')
SRR3879605_counts<- read.csv('/projectnb/bf528/users/lava_lamp/project_4/esaake_pr4/p4/results/count_SRR3879605.csv', col.names=c('barcode', 'count'), sep='')
SRR3879606_counts<- read.csv('/projectnb/bf528/users/lava_lamp/project_4/esaake_pr4/p4/results/count_SRR3879606.csv', col.names=c('barcode', 'count'), sep='')



#Cleaning data => Removing null values (nas) from our data 
SRR3879604_counts <- SRR3879604_counts %>% filter(!is.na(count))
SRR3879605_counts <- SRR3879605_counts %>% filter(!is.na(count))
SRR3879606_counts <- SRR3879606_counts %>% filter(!is.na(count))


#Computing the mean for each sample group
SRR3879604_mean <- mean(SRR3879604_counts$count)
SRR3879605_mean <- mean(SRR3879605_counts$count)
SRR3879606_mean <- mean(SRR3879606_counts$count)

#sorting
SRR3879604_counts <-SRR3879604_counts %>% arrange(desc(count))
SRR3879605_counts<- SRR3879605_counts %>% arrange(desc(count))
SRR3879606_counts <- SRR3879606_counts %>% arrange(desc(count))

#filtering out those values greater than the mean
SRR3879604_whitelisted <- SRR3879604_counts %>% filter(SRR3879604_counts$count>SRR3879604_mean)
SRR3879605_whitelisted <- SRR3879605_counts %>% filter(SRR3879605_counts$count>SRR3879605_mean)
SRR3879606_whitelisted <- SRR3879606_counts %>% filter(SRR3879606_counts$count>SRR3879606_mean)

#Saving data version for saving- it avoids general lost of data when somthing goes wrong with this dataframe
SRR3879604_whitelisted<-data.frame(SRR3879604_whitelisted$barcode)
SRR3879605_whitelisted<-data.frame(SRR3879605_whitelisted$barcode)
SRR3879606_whitelisted<-data.frame(SRR3879606_whitelisted$barcode)

#Plese remind ME:
  #Use write_csv instead of write.csv ... Rememeber to load 'tidyverse'
  #Takes out the unnecessary string formating
  #No rownames
#write_csv(SRR3879604_whitelisted, '/projectnb/bf528/users/lava_lamp/project_4/esaake_pr4/p4/salmon/SRR3879604_whitelist.txt', col_names = FALSE)
#write_csv(SRR3879605_whitelisted, '/projectnb/bf528/users/lava_lamp/project_4/esaake_pr4/p4/salmon/SRR3879605_whitelist.txt', col_names = FALSE)
#write_csv(SRR3879606_whitelisted, '/projectnb/bf528/users/lava_lamp/project_4/esaake_pr4/p4/salmon/SRR3879606_whitelist.txt', col_names = FALSE)


#Data distribution Plots
#commented to avoid rewrite
#==========================================>
#png('./plots/commulative_plot_SRR3879604.png',  width = 900, height = 600 )
plot(ecdf(SRR3879604_counts$count), cex=0,  main="Commulative Distribution Plot- SRR3879604", xlab="Count", ylab="percent")
#dev.off()

#png('./plots/commulative_plot_SRR3879605.png',  width = 900, height = 600 )
plot(ecdf(SRR3879605_counts$count), cex=0, main="Commulative Distribution Plot- SRR3879605", xlab="Count", ylab="percent")
#dev.off()


#png('./plots/commulative_plot_SRR3879606.png',  width = 900, height = 600 )
plot(ecdf(SRR3879606_counts$count), cex=0,  main="Commulative Distribution Plot- SRR3879606", xlab="Count", ylab="percent")
#dev.off()



