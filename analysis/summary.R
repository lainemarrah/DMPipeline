library(dplyr)
library(ggplot2)
library(stringr)

#first arg: dir with DMFinder result (summary.sh generated) files
#second arg: name of text file with cancer names (in alphabetical order), should be in same dir as result files
args <- commandArgs(trailingOnly = TRUE)

#downloading and preparing args
dir = args[1]
files = list.files(dir)
files = files[-which(files==args[2])]
cancerfile = paste0(dir,"/",args[2])
cancers = data.frame(read.table(cancerfile, sep="\n"))
rm(cancerfile)

avillines = function(df){
  df_lines = df %>% filter(!(str_detect(V1, '^DM([1-9]|_)')))
  df_lines = df_lines[,1]
  df_out = data.frame(df, rep(NA, length(df)), rep(NA, length(df)))
  df_out = df_out %>% filter(!(str_detect(V1, '^DM_index')))
 
  #assigning cell line to each DM
  for (i in 1:nrow(df_out)){
    if (df_out[i,1] %in% df_lines){
      line=df_out[i,1]
    } else if (df_out[i,1] %in% df_lines==FALSE){
      df_out[i,2] = line
    }
  }
  df_out = df_out %>% filter((str_detect(V1, '^DM([1-9]|_)')))
  colnames(df_out) = c("DMs", "Cell_Line", "AVIL")
 
  #determining which DMs contain AVIL
  for (i in 1:nrow(df_out)){
    dmvec = df_out$DMs
    dmvec = substr(dmvec, 11, nchar(dmvec))
    dmlist = strsplit(dmvec,",")
    for (i in 1:length(dmlist)){
      suppressWarnings(dmlist[[i]] <- as.numeric(dmlist[[i]]))
      if (dmlist[[i]][1]<=58191159 & dmlist[[i]][2]>=58212517){
        df_out[i,3] = TRUE
      } else {
        df_out[i,3] = FALSE
      }
    }
  }
  rm(line)
  rm(dmvec)
  rm(dmlist)
  return(df_out)
}

dmcount = function(df){
  df_lines = df %>% filter(!(str_detect(V1, '^DM([1-9]|_)')))
  df_lines = df_lines[,1]
  lines_summary = data.frame(df_lines, rep(NA, length(df_lines)), rep(NA, length(df_lines)))
  colnames(lines_summary)=c("Cell_Lines", "DM_Count", "AVIL")
  al = avillines(df)
  
  #counting DMs for each line present
  for(i in 2:nrow(al)){
    if(al[i,2]!=al[i-1,2]){
      num = substr(al[i-1,1], 3, 3)
      line = al[i-1, 2]
      lines_summary[(grep(line,df_lines)),2] = num
    } else if(i==nrow(al)){
      num = substr(al[i,1], 3, 3)
      line = al[i,2]
      lines_summary[(grep(line,df_lines)),2] = num
    }
    if(al[i,3]==TRUE){lines_summary[which(lines_summary$Cell_Line==al[i,2]),3]=TRUE}
  }
  
  #entering lines with 0 DMs and 0 AVIL
  for(i in 1:nrow(lines_summary)){
    if(is.na(lines_summary[i,2])){lines_summary[i,2] = 0}
    if(is.na(lines_summary[i,3])){lines_summary[i,3] = FALSE}
  }  
  
  rm(df_lines)
  rm(line)
  rm(al)
  rm(num)
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
all_out = data.frame(NA, NA, NA, NA)
colnames(all_out) = c("Cell_Lines", "DM_Count", "AVIL", "Cancer")
for (i in 1:length(files)){
  df = read.table(paste0(dir,"/",files[i]), sep="\t")
  df_out = dmsummary(df, cancers[i,1])
  all_out = rbind(all_out, df_out)
}
all_out = all_out %>% na.omit()

#write output tab-separated file to same dir, optionally can change separator in sep argument
write.table(all_out, paste0(dir,"/summary.txt"), sep="\t")

