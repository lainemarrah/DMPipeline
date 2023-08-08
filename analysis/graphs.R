#do this after going through summary.R in the same env, which should create all_out dataframe
#arg: desired directory for output graphs
args <- commandArgs(trailingOnly = TRUE)
setwd(args[1])

#showing frequency of 0, 1, 2, or 3 DM cell lines for each cancer type
ggplot(all_out, aes(x=Cancer, fill=DM_Count))+geom_bar()+ylab("Cell Line Count")+scale_x_discrete(guide = guide_axis(angle = 90))
ggsave("bar_dm_count_summary.png")

#showing frequency of AVIL presence in DM, good after previous to compare
ggplot(all_out, aes(x=Cancer, fill=AVIL))+geom_bar()+ylab("Cell Line Count")+scale_x_discrete(guide = guide_axis(angle = 90))
ggsave("bar_avil_count_summary.png")
#same but removing cell lines with 0 DMs
#filter(all_out, DM_Count!=0) %>% ggplot(aes(x=Cancer, fill=AVIL))+geom_bar()+ylab("Cell Line With >0 DMs Count ")+scale_x_discrete(guide = guide_axis(angle = 90))

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

#quick metrics
metrics = function(){
  print(paste0("Total cell lines checked: ",nrow(all_out)))
  print(paste0("Cell lines with 1+ DMs: ",nrow(filter(all_out, DM_Count>0))))
  print(paste0("Percent of cell lines with 1+ DMs: ",(nrow(filter(all_out, DM_Count>0))/nrow(all_out))*100, "%"))
  print(paste0("Cell lines with an AVIL DM: ",nrow(filter(all_out, AVIL==TRUE))))
}
metrics()
