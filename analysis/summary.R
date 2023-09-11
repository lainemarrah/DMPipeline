library(dplyr)
library(ggplot2)
library(stringr)
library(biomaRt)

#first arg: dir with DMFinder result (summary.sh generated) files for each cancer
#second arg: name of text file with cancer names (in alphabetical order), should be in same dir as result files
#third arg (optional): name of gene of interest
gene = NA
args = commandArgs(trailingOnly = TRUE)

#downloading and preparing args
dir = args[1]
files = list.files(dir)
files = files[-which(files==args[2])]
cancerfile = paste0(dir,"/",args[2])
cancers = data.frame(read.table(cancerfile, sep="\n"))
rm(cancerfile)
gene = args[3]

dmcount = function(df){
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
  
  lines_summary = data.frame(df_lines, rep(NA, length(df_lines)))
  colnames(lines_summary)=c("Cell_Lines", "DM_Count")
  
  #no dms predicted, exit function
  if (nrow(df_out)==0){
    lines_summary[,2] = 0
    if(!is.na(gene)){
      lines_summary = data.frame(lines_summary, rep(FALSE, nrow(lines_summary)))
      colnames(lines_summary)[3] = gene
    }
    return(lines_summary)
  }
  
  extras = lines_summary%>%filter(Cell_Lines%nin%df_out$Cell_Line)
  lines_summary = lines_summary%>%filter(Cell_Lines%in%df_out$Cell_Line)
  
  #counting DMs per cell line
  if (nrow(df_out)!=0){ #confirming DM presence
    for(i in 2:nrow(df_out)){
      if(df_out[i,2]!=df_out[i-1,2]){
        num = substr(df_out[i-1,1], 3, 3)
        line = df_out[i-1, 2]
        lines_summary[which(lines_summary$Cell_Lines==line),2] = num
      } else if(i==nrow(df_out)){
        num = substr(df_out[i,1], 3, 3)
        line = df_out[i,2]
        lines_summary[which(lines_summary$Cell_Lines==line),2] = num
      }
      #if(df_out[i,3]==TRUE){lines_summary[which(lines_summary$Cell_Line==df_out[i,2]),3]=TRUE}
    }
  }
  lines_summary = rbind(lines_summary, extras)
  
  if(!is.na(gene)){
    #add column for gene presence
    lines_summary = cbind(lines_summary, rep(NA, nrow(lines_summary)))
    colnames(lines_summary)[3] = gene
    
    #getting gene coordinates
    m = useMart('ensembl', dataset='hsapiens_gene_ensembl', "https://feb2014.archive.ensembl.org") #hg19
    martdf = getBM(mart=m, attributes=c('hgnc_symbol', 'description', 'chromosome_name',
                                        'start_position', 'end_position', 'strand', 'ensembl_gene_id'), 
                   filters='hgnc_symbol', values=gene)
    startpos = martdf[1,4]
    endpos = martdf[1,5]
    
    #finding gene in DMs
    lineswdm = unique(df_out$Cell_Line)
    for (j in 1:length(lineswdm)){ #iterating so i can separate cell lines
      celllinedms = df_out%>%filter(Cell_Line==lineswdm[j]) 
      dmnumber = substr(celllinedms[nrow(celllinedms),1],3,3)
      infovec = c()
      for (k in 1:dmnumber){ #start and finish coordinates of individual dm
        specificdm = celllinedms[k==substr(celllinedms[,1],3,3),]
        specificdmlist = strsplit(specificdm$DMs,",")
        mincoords = rep(NA, length(specificdmlist))
        maxcoords = rep(NA, length(specificdmlist))
        for (l in 1:length(specificdmlist)){
          mincoords[l] = suppressWarnings(as.numeric(specificdmlist[[l]][3]))
          maxcoords[l] = suppressWarnings(as.numeric(specificdmlist[[l]][4]))
        }
        start = min(mincoords)
        finish = max(maxcoords)
        if (min(mincoords)<=startpos & max(maxcoords)>=endpos){ #checking for avil coords
          lines_summary[which(lines_summary$Cell_Lines==lineswdm[j]),3] = TRUE
        } 
      } 
    }
    
    for(i in 1:nrow(lines_summary)){
      if(is.na(lines_summary[i,3])){lines_summary[i,3] = FALSE}
    }
  }

  #entering lines with 0 DMs 
  for(i in 1:nrow(lines_summary)){
    if(is.na(lines_summary[i,2])){lines_summary[i,2] = 0}
  }
  return(lines_summary)
}

#repeating summary function for each result file in directory
if(is.na(gene)){
  all_out = data.frame(NA, NA, NA)
  colnames(all_out) = c("Cell_Lines", "DM_Count", "Cancer")
  for (i in 1:(length(files))){
    df = read.table(paste0(dir,"/",files[i]), sep="\t")
    canceroutput = dmcount(df)
    Cancer = rep(cancers[i,1], nrow(canceroutput))
    canceroutput = cbind(canceroutput, Cancer)
    all_out = rbind(all_out, canceroutput)
  }
} else {
  all_out = data.frame(NA, NA, NA, NA)
  colnames(all_out) = c("Cell_Lines", "DM_Count", gene, "Cancer")
  for (i in 1:(length(files))){
    df = read.table(paste0(dir,"/",files[i]), sep="\t")
    canceroutput = dmcount(df)
    Cancer = rep(cancers[i,1], nrow(canceroutput))
    canceroutput = cbind(canceroutput, Cancer)
    all_out = rbind(all_out, canceroutput)
  }
}
all_out = all_out %>% na.omit()

#write output tab-separated file to same dir, optionally can change separator in sep argument
#write.table(all_out, paste0(dir,"/summary.txt"), sep="\t")
