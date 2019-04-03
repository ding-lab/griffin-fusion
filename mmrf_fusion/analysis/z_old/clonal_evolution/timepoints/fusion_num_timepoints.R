compare_df <- read.table("/Users/lijunyao/Documents/PhD_LAB/compare_fusions_timepoint.txt",header=TRUE)
sample_ID <- unique(compare_df$mmrf)

#get how many timepoints for each sample
sample_df <- read.table("/Users/lijunyao/Documents/PhD_LAB/sample_list.807.txt",header=FALSE)
tp_list = as.vector(table(sample_df$V1)[table(sample_df$V1) > 1])

df <- read.table("/Users/lijunyao/Documents/PhD_LAB/fusion_df.txt",header = TRUE)
nrows=dim(unique(subset(df,df$has_secondary==1)[,c(1,4)]))[1]
out_df_header <- c("mmrf","time_point","fusions_num")
out_df <- as.data.frame(matrix(NA,nrows,length(out_df_header)))
names(out_df) <- out_df_header
tps <- c("time_1","time_2","time_3","time_4")
n <- 0
for (i in 1:length(sample_ID)) {
  x=tp_list[i]
  for (j in 1:x){
    n <- n+1
    sample_subset <- subset(compare_df,mmrf==sample_ID[i])[,tps[j]]
    sample_subset <- sample_subset[sample_subset!=0 & sample_subset!="NA" ]
    fusion_num <- length(sample_subset)
    out_df[n,] <- c(as.character(sample_ID[i]),j,as.vector(fusion_num))
  }
}
write.table(out_df,file="/Users/lijunyao/Documents/PhD_LAB/fusion_comparison_plot.txt",row.names = FALSE,sep = "\t",quote = FALSE)

plot_df <- read.table("/Users/lijunyao/Documents/PhD_LAB/fusion_comparison_plot.txt",header = TRUE)

jpeg(file = "/Users/lijunyao/Documents/PhD_LAB/fusion_timepoint_comparison.jpeg")
#plot(plot_df$time_point,plot_df$fusions_num,xlab="time points",ylab="number of fusions")
for (z in 1:length(sample_ID)) {
  sample_plot <- subset(plot_df,mmrf==sample_ID[z])
  if (z==1) {
    
    plot(sample_plot[,"time_point"],sample_plot[,"fusions_num"],xlab="time points",ylab="number of fusions",type="b",xlim=c(1,4),ylim=c(0,25))
  }
  else {
    lines(sample_plot[,"time_point"],sample_plot[,"fusions_num"],col=z,type="b")
  }
}
dev.off()