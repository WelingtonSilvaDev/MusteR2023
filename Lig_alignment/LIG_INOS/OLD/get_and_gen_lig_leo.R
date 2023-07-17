library(tidyverse)

if(1){
  filename = "LIGS/leo-ligs-v1.csv"
  lig.table = read_delim(file=filename,delim=",")
  scriptname = "gapin-pdb2align.R"
  sep =","
  for (i in 1:dim(lig.table)[1]){
  #for (i in 1:1){
    parn = paste0(lig.table[i,2],sep,lig.table[i,3],sep,lig.table[i,4],sep,lig.table[i,5])
    command = paste("Rscript",scriptname,lig.table$pdb_title[i],parn)
    print(paste("RUNNING",command))
    try(system(command))
    #browser()
  }
}