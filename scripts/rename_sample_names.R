#!/usr/bin/env Rscript

# library(tidyverse)
# library(readxl)
# library(ComplexHeatmap)
# library(PCAtools)
# library(ggplot2)
# library(colorRamp2)
# library(ggalt)
# date_var <- as.Date(date(), format = "%a %b %d %H:%M:%S %Y")

args <- commandArgs(TRUE)
print(args[1])

# Data load ---------------------------------------------------------------
filepath <- file.path(args[1])
# filepath <- file.path("D:/Terra-Informatix/Rutendo-CHARM/histoSampleList_edt.csv")
outDir <- dirname(filepath)
filename <- paste0(sub('\\..*$', '',  basename(filepath)),"_new.csv")
outpath <- file.path(paste(outDir,filename,sep = "/"))


df1 <- read.csv(filepath)

# Rename sample ids -------------------------------------------------------

len <- nrow(df1)
cols <- ncol(df1)

nudf <- data.frame(matrix(NA,nrow = len, ncol = 4))
renameDF <- data.frame(matrix(NA,nrow = len, ncol = 4))


for (i in (1:len)){ #print(i)}
  oldID <- df1[i,1]
  nuID <- gsub("_","-",oldID) #str_replace_all(oldID,"_","-")
  
  
  read1 <- df1[i,2]
  r1 <- gsub("-","_",read1)
  repID <- gsub("-","_",oldID)
  rd1 <- gsub(repID,nuID,r1) #str_replace(read1,oldID,nuID)
  
  read2 <- df1[i,3]
  r2 <- gsub("-","_",read2)
  rd2 <- gsub(repID,nuID,r2) #str_replace(read2,oldID,nuID)
    
  
  nudf[i,1] <- oldID
  nudf[i,2] <- nuID
  nudf[i,3] <- rd1
  nudf[i,4] <- rd2
  
  
  renameDF[i,1] <- read1
  renameDF[i,2] <- rd1
  renameDF[i,3] <- read2
  renameDF[i,4] <- rd2
  
}

names(nudf) <- c("old_id","sampleID","read1","read2")
# head(nudf)

# Write new sample sheet to file ------------------------------------------


write.csv(nudf, file = outpath,row.names=F,quote = F)


# Rename fastq files using new ID -----------------------------------------
names(renameDF) <- c("read1","newread1","read2","newread2")
# head(renameDF)
write.csv(renameDF, file = "rename_reads_test.csv",row.names=F,quote = F)

for(j in (1:nrow(renameDF))){ #print(j)}
  
  # get path
  dir1 <- dirname(renameDF[j,1])
  if(dir.exists(dir1)){
    setwd(dir = dir1)
    # read1
    old1=basename(renameDF[j,1])
    new1=basename(renameDF[j,2])
    file.rename(old1,new1)
    # read2
    old2=basename(renameDF[j,3])
    new2=basename(renameDF[j,4])
    file.rename(old2,new2)
  }
  
  
  
  
}


# End ---------------------------------------------------------------------


