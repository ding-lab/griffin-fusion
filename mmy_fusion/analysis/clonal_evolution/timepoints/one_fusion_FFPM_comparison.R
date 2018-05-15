compare_df <- read.table("/Users/lijunyao/Documents/PhD_LAB/compare_fusions_timepoint.txt",header=TRUE)
#change the fusion of interest here
fusion_subset <- subset(compare_df,fusion=="IGH@--WHSC1")
fusion_sample_ID <- unique(fusion_subset[,"mmrf"])

#determine the number of rows in the output dataframe
df <- read.table("/Users/lijunyao/Documents/PhD_LAB/fusion_df.txt",header = TRUE)
row_n <- 0
for (m in 1:length(sample_fusion_ID)) {
  fusion_l <- length(grep(sample_fusion_ID[m],unique(subset(df,df$has_secondary==1)[,c(1,4)])[,1]))
  row_n <- row_n + fusion_l
}

sample_header <- c("mmrf","time_point","FFPM")
one_sample_df <- as.data.frame(matrix(NA,row_n,length(sample_header)))
names(one_sample_df) <- sample_header
tps <- c("time_1","time_2","time_3","time_4")

n <- 0
for (i in 1:length(fusion_sample_ID)){
  sample_subset <- subset(fusion_subset,mmrf==fusion_sample_ID[i])
  unlist(sample_subset, use.names=FALSE)
  NA_sample_subset <- sample_subset[sample_subset=="NA"]
  fusion_num <- (4-length(NA_sample_subset))
  for (j in 1:fusion_num){
    n <- n+1
    one_sample_df[n,] <- c(as.character(fusion_sample_ID[i]),j,sample_subset[2+j])
  }
}

#one_sample_df
jpeg(file = "/Users/lijunyao/Documents/PhD_LAB/one_fusion_timepoint_comparison.jpeg")
#plot(plot_df$time_point,plot_df$fusions_num,xlab="time points",ylab="number of fusions")
for (z in 1:length(fusion_sample_ID)) {
  sample_plot <- subset(one_sample_df,mmrf==fusion_sample_ID[z])
  if (z==1) {
    plot(sample_plot[,"time_point"],sample_plot[,"FFPM"],xlab="time points",ylab="FFPM",type="b",xlim=c(1,4),ylim=c(0,10))
  }
  else {
    lines(sample_plot[,"time_point"],sample_plot[,"FFPM"],col=z,type="b")
  }
}
dev.off()