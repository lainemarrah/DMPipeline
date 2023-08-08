#do this after going through summary.R in the same env, which should create all_out dataframe
#can comment out ggsave() lines if you don't want to download the graph files
#arg: desired directory for output graphs
library(dplyr)
library(ggplot2)

#first arg: dir for graphs to go into
gene = NA
args <- commandArgs(trailingOnly = TRUE)
setwd(args[1])
if(!is.na(colnames(all_out)[4])){
  gene = colnames(all_out)[3]
}

#showing frequency of 0, 1, 2, or 3 DM cell lines for each cancer type
ggplot(all_out, aes(x=Cancer, fill=DM_Count))+geom_bar()+ylab("Cell Line Count")+scale_x_discrete(guide = guide_axis(angle = 90))
ggsave("bar_dm_count_summary.png")

#same thing but showing proportions instead of counts
ggplot(all_out, aes(x=Cancer, fill=DM_Count))+geom_bar(position="fill")+ylab("Cell Line Proportion")+scale_x_discrete(guide = guide_axis(angle = 90))
ggsave("bar_dm_proportion_summary.png")

#showing mean DMs per cell line for each cancer
ggplot(all_out, aes(x=Cancer, y=DM_Count, group=Cancer))+geom_bar(stat = "summary", fun = "mean")+ylab("Mean DMs per Cell Line")+scale_x_discrete(guide = guide_axis(angle = 90))
ggsave("bar_cancer_average_dm.png")

#pie chart of dm count
tblc = data.frame(table(all_out$DM_Count))
colnames(tblc)[1] = "DM count"
ggplot(tblc, aes(x="", y=Freq, fill=`DM count`)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  theme_void()
ggsave("pie_dm_count_summary.png")

metrics = function(){
  print(paste0("Total cell lines checked: ",nrow(all_out)))
  print(paste0("Cell lines with 1+ DMs: ",nrow(filter(all_out, DM_Count>0))))
  print(paste0("Percent of cell lines with 1+ DMs: ",(nrow(filter(all_out, DM_Count>0))/nrow(all_out))*100, "%"))
}

if(!is.na(gene)){
  #showing frequency of gene presence in DM, good after previous to compare
  ggplot(all_out, aes(x=Cancer, fill=all_out[,3]))+geom_bar()+labs(y="Cell Line Count", fill=gene)+scale_x_discrete(guide = guide_axis(angle = 90))
  ggsave("bar_gene_count_summary.png")
  
  #gene lines to screenshot
  print(all_out[which(all_out[,3]==TRUE),])
  
  metrics = function(){
    print(paste0("Total cell lines checked: ",nrow(all_out)))
    print(paste0("Cell lines with 1+ DMs: ",nrow(filter(all_out, DM_Count>0))))
    print(paste0("Percent of cell lines with 1+ DMs: ",(nrow(filter(all_out, DM_Count>0))/nrow(all_out))*100, "%"))
    print(paste0("Cell lines with DM containing ",gene,": ",nrow(all_out[which(all_out[,3]==TRUE),])))
  }
}

#quick printed summary
metrics()
