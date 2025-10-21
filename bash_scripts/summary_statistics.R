# File name 
reference <- snakemake@wildcards[[1]]
isolate <- snakemake@wildcards[[2]]
depth_file <- read.csv(snakemake@input[[1]])

calculate_overall_depth <- function(depth_file,isolate,reference){
  isolate <- isolate
  reference <- reference 
  mean <-  sum(depth_file$Total_Depth) / nrow(depth_file)
  median <- median(depth_file$Total_Depth)  
  non_zero <- subset(depth_file$Total_Depth, depth_file$Total_Depth > 0)
  mean_non_zero <-  sum(non_zero) / length(non_zero)
  median_non_zero <- median(non_zero)  
  results <- cbind.data.frame(isolate,reference,mean,median,mean_non_zero,median_non_zero)
  return(results)
}

summary_statistics <- calculate_overall_depth(depth_file = depth_file, isolate = isolate, reference = reference)

write.table(summary_statistics,file = snakemake@output[[1]], sep = "\t", quote = FALSE, row.names = FALSE, na = "")