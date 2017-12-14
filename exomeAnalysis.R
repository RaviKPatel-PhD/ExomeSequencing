
#Analysis of Exome data

#in Galaxy I got the VCF (see the word doc for note on how I got there)
# converted the vcf to tabular in Galaxy using VCFtools
#I moved the data to the working directory, downloaded directly from Galaxy
setwd("C:/Users/Ravi Patel/Desktop/Thomas Smithgall- Thesis Lab/Experiments/0046 2017-10-05 Exome data/Bam and VCF files/R workspace and tabular files")
#only need the dplyr package
library(dplyr)  
#Loading the data and eliminating the "noise" (DP>4 means read depth >4)
#also adding an additional column called "sample" which just contains the cell line that the data comes from
#\t=tab in the read.csv notation. I had to specify this because the empty spaces in the tabular files are being skipped when read using read.table 

ExomeSamples<-c("MV411", "MV411R2", "MV411R3", "MOLM13", "MOLM13R3", "MOLM14", "MOLM14R1", "MOLM14R2", "MOLM14R3", "MOLM14ac220R3")

Files<-c("Galaxy336-[VCFtoTab-delimited__on_data_165].tabular", "Galaxy337-[VCFtoTab-delimited__on_data_182].tabular", 
         "Galaxy338-[VCFtoTab-delimited__on_data_199].tabular", "Galaxy339-[VCFtoTab-delimited__on_data_216].tabular", 
         "Galaxy340-[VCFtoTab-delimited__on_data_233].tabular", "Galaxy341-[VCFtoTab-delimited__on_data_250].tabular",
         "Galaxy342-[VCFtoTab-delimited__on_data_267].tabular", "Galaxy343-[VCFtoTab-delimited__on_data_284].tabular", 
         "Galaxy344-[VCFtoTab-delimited__on_data_301].tabular", "Galaxy345-[VCFtoTab-delimited__on_data_318].tabular")

loadData<-function(File, Samp){
  Data<-read.table(File, header = T, fill=T, na.strings="", sep="\t") 
  Data$DP<-as.numeric(as.character(Data$DP))
  Data<-Data[Data$DP>4,]
  Data<-mutate(Data, sample=Samp)
  return(Data)
}

loadDataset<-function(files){
  listDF<-list()
  for(i in 1:length(files)){
    listDF[[i]]<-loadData(File = files[i],i)
  }
  return(listDF)
}

parsExomeData<-loadDataset(Files)


#subset only the missense mutations

SubsetEFFbyTxt<-function(files, text){
  SubsetData<-list()
  for(i in 1:length(files)){
    SubsetData[[i]]<-subset(files[[i]], grepl(pattern = text, x=files[[i]]$EFF))
  }
  return(SubsetData)
}

SubExomeData<-SubsetEFFbyTxt(parsExomeData, "missense")


#Renaming the elements of the list and pulling out individual data frames instead of keeping everything in the list

MV411<-SubExomeData[[1]]
MV411R2<-SubExomeData[[2]]
MV411R3<-SubExomeData[[3]]
MOLM13<-SubExomeData[[4]]
MOLM13R3<-SubExomeData[[5]]
MOLM14<-SubExomeData[[6]]
MOLM14R1<-SubExomeData[[7]]
MOLM14R2<-SubExomeData[[8]]
MOLM14R3<-SubExomeData[[9]]
MOLM14ac220R3<-SubExomeData[[10]]


#Subsetting such that mutations that are in the the resistant cell line and not in the parent cell line and vice versa

MV411minR2<-subset(MV411, !(POS %in% MV411R2$POS))
MV411R2unique<-subset(MV411R2, !(POS %in% MV411$POS))
MV411minR3<-subset(MV411, !(POS %in% MV411R3$POS))
MV411R3unique<-subset(MV411R3, !(POS %in% MV411$POS))

MOLM13minR3<-subset(MOLM13, !(POS %in% MOLM13R3$POS))
MOLM13R3unique<-subset(MOLM13R3, !(POS %in% MOLM13$POS))

MOLM14minR1<-subset(MOLM14, !(POS %in% MOLM14R1$POS))
MOLM14R1unique<-subset(MOLM14R1, !(POS %in% MOLM14$POS))
MOLM14minR2<-subset(MOLM14, !(POS %in% MOLM14R2$POS))
MOLM14R2unique<-subset(MOLM14R2, !(POS %in% MOLM14$POS))
MOLM14minR3<-subset(MOLM14, !(POS %in% MOLM14R3$POS))
MOLM14R3unique<-subset(MOLM14R3, !(POS %in% MOLM14$POS))
MOLM14minac220R3<-subset(MOLM14, !(POS %in% MOLM14ac220R3$POS))
MOLM14ac220R3unique<-subset(MOLM14ac220R3, !(POS %in% MOLM14$POS))

#still have >1000 obs per data frame
#maybe use ingenuity variant analysis?

#testing if MV411 R2 is actually MOLM14 or MOLM13 cell line
MOLM13minMV411R2<-subset(MOLM13, !(POS %in% MV411R2$POS))
MV411R2ifMOLM13unique<-subset(MV411R2, !(POS %in% MOLM13$POS))
MOLM14minMV411R2<-subset(MOLM14, !(POS %in% MV411R2$POS))
MV411R2ifMOLM14unique<-subset(MV411R2, !(POS %in% MOLM14$POS))

