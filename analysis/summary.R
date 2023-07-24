library(dplyr)
library(ggplot2)

avillines = function(df){
  df = df %>% filter(!(str_detect(V1, '^DM_index')))
  df_lines = df %>% filter(!(str_detect(V1, '^DM')))
  df_lines = df_lines[,1]
  df_out = data.frame(df, rep(NA, length(df)), rep(NA, length(df)))
  
  for (i in 1:nrow(df)){
    if (df[i,1] %in% df_lines){
      line=df[i,1]
    } else if (df[i,1] %in% df_lines==FALSE){
      df_out[i,2] = line
    }
  }
  df_out = df_out %>% filter((str_detect(V1, '^DM')))
  colnames(df_out) = c("DMs", "Cell_Line", "AVIL")
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
  rm(df_lines)
  rm(line)
  rm(dmvec)
  rm(dmlist)
  return(df_out)
}

#calls avillines
dmcount = function(df){
  df = df %>% filter(!(str_detect(V1, '^DM_index')))
  df_lines = df %>% filter(!(str_detect(V1, '^DM')))
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


#gbm = read.table("CCLE-GBM-results.csv", sep="\t")
#gbm_out = dmsummary(gbm, "Glioblastoma")
