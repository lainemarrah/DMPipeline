library(dplyr)
library(ggplot2)
library(stringr)
library(biomaRt)

#first arg: dir with DMFinder result (summary.sh generated) files for each cancer
#second arg: name of text file with cancer names (in alphabetical order), should be in same dir as result files
#third arg (optional): name of gene of interest
gene = NA
args = c("/Users/lainemarrah/Desktop/LiLab/DMFinder_Results", "cancers.txt", "AVIL")
#args <- commandArgs(trailingOnly = TRUE)

#downloading and preparing args
dir = args[1]
files = list.files(dir)
files = files[-which(files==args[2])]
cancerfile = paste0(dir,"/",args[2])
cancers = data.frame(read.table(cancerfile, sep="\n"))
rm(cancerfile)
gene = args[3]
  
dmlines = function(df){
  df_lines = df %>% filter(!(str_detect(V1, '^DM([1-9]|_)')))
  df_lines = df_lines[,1]
  df_out = data.frame(df, rep(NA, length(df)))
  df_out = df_out %>% filter(!(str_detect(V1, '^DM_index')))
  colnames(df_out) = c("DMs", "Cell_Line")
  
  #assigning cell line to each DM
  for (i in 1:nrow(df_out)){
    if (df_out[i,1] %in% df_lines){
      line=df_out[i,1]
    } else if (df_out[i,1] %in% df_lines==FALSE){
      df_out[i,2] = line
    }
  }
  df_out = df_out %>% filter((str_detect(DMs, '^DM([1-9]|_)')))
  
  if(!is.na(gene)){
    df_out = cbind(df_out, rep(NA, nrow(df_out)))
    colnames(df_out)[3] = gene
    #getting gene coordinates
    m = useMart('ensembl', dataset='hsapiens_gene_ensembl', "https://feb2014.archive.ensembl.org") #hg19
    martdf = getBM(mart=m, attributes=c('hgnc_symbol', 'description', 'chromosome_name',
                                        'start_position', 'end_position', 'strand', 'ensembl_gene_id'), 
                   filters='hgnc_symbol', values=gene)
    startpos = martdf[1,4]
    endpos = martdf[1,5]
    rm(m)
    rm(martdf)
    #return if no DMs   
    if (nrow(df_out)==0){
      return(df_out)
    }
    #determining which DMs contain the gene if it exists
    for (i in 1:nrow(df_out)){
      dmvec = df_out$DMs
      dmvec = substr(dmvec, 11, nchar(dmvec))
      dmlist = strsplit(dmvec,",")
      for (i in 1:length(dmlist)){
        suppressWarnings(dmlist[[i]] <- as.numeric(dmlist[[i]]))
        if (dmlist[[i]][1]<=startpos & dmlist[[i]][2]>=endpos){
          df_out[i,3] = TRUE
        } else {
          df_out[i,3] = FALSE
        }
      }
    }
  }
  return(df_out)
}
  
dmcount = function(df){
  df_lines = df %>% filter(!(str_detect(V1, '^DM([1-9]|_)')))
  df_lines = df_lines[,1]
  lines_summary = data.frame(df_lines, rep(NA, length(df_lines)))
  colnames(lines_summary)=c("Cell_Lines", "DM_Count")
  dl = dmlines(df)
  if (nrow(dl)!=0){ #confirming DM presence
    #counting DMs for each line present
    for(i in 2:nrow(dl)){
      if(dl[i,2]!=dl[i-1,2]){
        num = substr(dl[i-1,1], 3, 3)
        line = dl[i-1, 2]
        lines_summary[(grep(line,df_lines)),2] = num
      } else if(i==nrow(dl)){
        num = substr(dl[i,1], 3, 3)
        line = dl[i,2]
        lines_summary[(grep(line,df_lines)),2] = num
      }
    }
  }
  
  #adding gene column if necessary
  if(!is.na(gene)){
    lines_summary = cbind(lines_summary, rep(NA, nrow(lines_summary)))
    colnames(lines_summary)[3] = gene
    if(nrow(dl)!=0){
      for(i in 1:nrow(dl)){
        if(dl[i,3]==TRUE){lines_summary[which(lines_summary$Cell_Line==dl[i,2]),3]=TRUE}
      }
    }
  }
  
  #entering lines with 0 DMs and 0 AVIL
  for(i in 1:nrow(lines_summary)){
    if(is.na(lines_summary[i,2])){lines_summary[i,2] = 0}
    if(!is.na(gene)){
      if(is.na(lines_summary[i,3])){lines_summary[i,3] = FALSE}
    }
  }
  return(lines_summary)
}

#calls dmcount
dmsummary = function(df, cancer) {
  #adding relevant cancer column
  df_out = dmcount(df)
  Cancer = rep(cancer, nrow(df_out))
  df_out = cbind(df_out, Cancer)
  rm(Cancer)
  return(df_out)
}

#repeating summary function for each result file in directory
if(is.na(gene)){
  all_out = data.frame(NA, NA, NA)
  colnames(all_out) = c("Cell_Lines", "DM_Count", "Cancer")
  for (i in 1:length(files)){
    print(cancers[i,1])
    df = read.table(paste0(dir,"/",files[i]), sep="\t")
    df_out = dmsummary(df, cancers[i,1])
    all_out = rbind(all_out, df_out)
  }
} else {
  all_out = data.frame(NA, NA, NA, NA)
  colnames(all_out) = c("Cell_Lines", "DM_Count", gene, "Cancer")
  for (i in 1:length(files)){
    print(cancers[i,1])
    df = read.table(paste0(dir,"/",files[i]), sep="\t")
    df_out = dmsummary(df, cancers[i,1])
    all_out = rbind(all_out, df_out)
  }
}
all_out = all_out %>% na.omit()

#write output tab-separated file to same dir, optionally can change separator in sep argument
#write.table(all_out, paste0(dir,"/summary.txt"), sep="\t")
