library("readr")
pfam_table<-read_tsv("data/11.MEBs/CP9_out.tsv")
pfam_table_df<-as.data.frame(pfam_table)
pfam_table_df_2 <- pfam_table_df[,-1]
rownames(pfam_table_df_2) <- pfam_table_df[,1]

pfam_only_table_df<-pfam_table_df_2[,75:18003]

write.table(pfam_only_table_df, file = "data/11.MEBs/CP9_pfam.tsv", sep = "\t", quote = FALSE, row.names = T)
