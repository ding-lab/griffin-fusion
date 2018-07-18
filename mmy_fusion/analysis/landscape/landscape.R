#library(ggplot2)

n <- 10 #number of top fusions to consider

top_fusions_names <- names(sort(table(primary_df$fusion), decreasing=TRUE)[1:n])
top_fusions_count <- as.vector(sort(table(primary_df$fusion), decreasing=TRUE)[1:n])
top_fusions_rate <- vector()
for(name in top_fusions_names){
  val <- subset(primary_df, fusion==name)$n_discordant
  rate <- mean(val >= 3 , na.rm=T)
  top_fusions_rate <- c(top_fusions_rate, rate)
}

top_callers <- names(sort(table(primary_df$Callers), decreasing=TRUE))
top_callers_count <- as.vector(sort(table(primary_df$Callers), decreasing=TRUE))
top_callers_validation_rate <- vector()
for(name in top_callers){
  val <- subset(primary_df, Callers==name)$n_discordant
  rate <- mean(val >= 3 , na.rm=T)
  top_callers_validation_rate <- c(top_callers_validation_rate, rate)
}

hyperdiploid_status <- unique(primary_df[,c("mmrf", "seqfish_Hyperdiploidy")])

#plots

##### top n fusions #####
#update coverge to show N of samples with WGS
plot_df <- data.frame(fusion=top_fusions_names, n_samples=top_fusions_count, validation_rate=format(top_fusions_rate,digits=1))
plot_df$fusion <- factor(plot_df$fusion, levels = rev(top_fusions_names))
ggplot(plot_df, aes(x=fusion, weight=n_samples, label=validation_rate)) + geom_bar() + coord_flip() + labs(x="Fusion name", y="Number of samples", title=paste0("Top ", n," Recurrent Fusions")) + theme_bw(base_size=20) + geom_text(aes(x = fusion, y = n_samples), hjust=1, color="white", size=5)

##### callers #####
#update coverge to show N of samples with WGS
plot_df <- data.frame(caller=top_callers, fusions=top_callers_count, validation_rate=format(top_callers_validation_rate,digits=1))
plot_df$caller <- factor(plot_df$caller, levels = rev(top_callers))
ggplot(plot_df, aes(x=caller, weight=fusions, label=validation_rate)) + geom_bar() + coord_flip() + labs(x="Fusion Caller Combination", y="Number of Fusions", title="Number of Fusions by Caller") + theme_bw(base_size=20) + geom_text(aes(x = caller, y = fusions), hjust=1, color="white", size=5)

##### distribution of fusions #####
n_fusions = vector()
for(sample in unique(samples[,1])){
  n_fusions <- c(n_fusions, dim(subset(primary_df, mmrf==sample))[1])
}
plot_df <- data.frame(srr = unique(samples[,1]), n_fusions = n_fusions)
med <- median(plot_df$n_fusions)
ggplot(plot_df, aes(x=n_fusions)) + geom_histogram(breaks=seq(min(plot_df[,2])-0.5,max(plot_df[,2])+0.5,1)) + labs(x="Number of Fusions", y="Number of Samples", title="Number of Fusions per Sample") + theme_bw(base_size=20) + geom_vline(xintercept=med, lty=2, color="white", size=1) + geom_text(label=paste0("Median=",med), x=med, y=50, angle=90, vjust=-1, color="white", size=5)

##### coverage by validation #####
plot_df <- subset(primary_df, !is.na(n_discordant))[,c("depthA","depthB","n_discordant")]
plot_df <- data.frame(plot_df, mean_cov=(plot_df$depthA + plot_df$depthB)/2, val=(plot_df$n_discordant>0))
plot_df$mean_cov[is.na(plot_df$mean_cov)] <- 0
ggplot(plot_df, aes(x=val, y=mean_cov)) + geom_violin(size=2) + theme_bw(base_size=20) + labs(x="Fusion Validated by WGS", y="WGS Coverage at Fusion", title="Fusion Validation and WGS Coverage") + geom_boxplot(width=.1, outlier.alpha = 0)

##### samples per fusion #####
plot_df <- data.frame(table(primary_df$fusion))
ggplot(plot_df, aes(x=Freq)) + geom_histogram(breaks=seq(min(plot_df[,2])-0.5,max(plot_df[,2])+0.5,1)) + theme_bw(base_size=20) + labs(x="Number of Samples", y="Number of Fusions", title="Number of Samples per Fusion")
ggplot(subset(plot_df, Freq>1), aes(x=Freq)) + geom_histogram(breaks=seq(min(plot_df[,2])-0.5,max(plot_df[,2])+0.5,1)) + theme_bw(base_size=20) + labs(x="Number of Samples", y="Number of Fusions", title="Number of Samples per Fusion")

