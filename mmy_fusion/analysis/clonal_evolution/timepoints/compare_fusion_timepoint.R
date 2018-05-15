#from 807 samples, get the sample list and #time_points
sample_df <- read.table("/Users/lijunyao/Documents/PhD_LAB/sample_list.807.txt",header=FALSE)
tp_sample_lists <- unique(sample_df[duplicated(sample_df$V1),]$V1)
tp_list = c()
for (y in 1:length(tp_sample_lists)) {
  tp <- dim((subset(sample_df,V1==tp_sample_lists[y])))[1]
  tp_list <-c (tp_list,tp)
}
#creat a matrix with NA
df <- read.table("/Users/lijunyao/Documents/PhD_LAB/fusion_df.txt",header = TRUE)
header <- c("mmrf","fusion","time_1","time_2","time_3","time_4")
df_sec <- subset(df,df$has_secondary==1)
nrows=dim(unique(df_sec[,c(1,3)]))[1]
output_df <- as.data.frame(matrix(NA,nrows,length(header)))
names(output_df) <- header
#loop through 
n <- 0 #initialize row num
for(i in 1:length(tp_sample_lists)) {
  sample_name<-as.character(tp_sample_lists[i])
  sample_subset <- subset(df_sec,mmrf==sample_name)
  fusion_list <- unique(sample_subset$fusion)
  for (j in 1:length(fusion_list)){
    n <- n+1
    fusion_name <- as.character(fusion_list[j])
    rep_subv<-subset(sample_subset,fusion==fusion_list[j])$sample_number
    rep <- as.numeric(c(1:tp_list[i]) %in% rep_subv)
    #replace the [0,1,1] with [0,FFPM_time_point2,FFPM_time_point3]
    for (z in 1:length(rep)){
      if (rep[z]==1) {
        rep[z] <- subset(sample_subset,fusion==fusion_list[j] & sample_number==z)$FFPM
      }
    }
    row <- c(sample_name,fusion_name,rep)
    output_df[n,c(1:(2+tp_list[i]))]<- row
  }
}
#if a fusion is dectected in a certain time point, the value is raw FFPM instead of 1.
write.table(output_df,quote=FALSE,row.names = FALSE,sep="\t",file="/Users/lijunyao/Documents/PhD_LAB/compare_fusions_timepoint.txt")