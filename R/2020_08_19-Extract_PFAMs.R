library(vroom)

########################################################### Read input data ###########################################################

pan_matrix<-vroom("/home/mirnis/Documentos/02.CP9/data/14.PFAM/Pfam_Clustering-20200722T194919Z-001/Pfam_Clustering/CP9_pfam_2.tsv", 
                  delim="\t")
clusters<-vroom("/home/mirnis/Documentos/02.CP9/data/14.PFAM/PFAM_clusters", delim="\t")
#dplyr::right_join, join the tables based on one column
pfam_matrix_clusters<-dplyr::right_join(pan_matrix, clusters, by = "Genomes")

############################################################ Function to get the pfams ################################################

get_unique_pfams<-function(df_pfam, groupOfInterest){
  gruop_1<-dplyr::filter(df_pfam, Groups == groupOfInterest)
  allGroupNames <- unique(df_pfam$Groups)
  otherGroups <- allGroupNames[which(allGroupNames  != groupOfInterest )]
  The_rest <- dplyr::filter(df_pfam, Groups %in%  otherGroups)
  uniquePFAMs <- c()
  for(i in 2:(length(gruop_1)-2)){
    if(isTRUE(sum(gruop_1[,i]) > 0 & sum(The_rest[,i]) == 0)){
      uniquePFAMs <- c(uniquePFAMs, colnames(gruop_1[i]))
    }
  }
  return(uniquePFAMs)
}

########################################################## Ejecutar la funcion ##########################################################

Pfam_group1<-get_unique_pfams(pfam_matrix_clusters, "Group 1")
Pfam_group2<-get_unique_pfams(pfam_matrix_clusters, "Group 2")
Pfam_group3<-get_unique_pfams(pfam_matrix_clusters, "Group 3")
Pfam_group4<-get_unique_pfams(pfam_matrix_clusters, "Group 4")
Pfam_group5<-get_unique_pfams(pfam_matrix_clusters, "Group 5")
Pfam_group6<-get_unique_pfams(pfam_matrix_clusters, "Group 6")
Pfam_group7<-get_unique_pfams(pfam_matrix_clusters, "Group 7")

save(pan_matrix, clusters, pfam_matrix_clusters, Pfam_group1, Pfam_group2, Pfam_group3, Pfam_group4, 
     Pfam_group5, Pfam_group6, Pfam_group7,   file="/home/mirnis/Documentos/02.CP9/data/00.Rdata/pfam.RData")



##########################################################################################################################################
############################################################ Cosas aprendidas ############################################################
##########################################################################################################################################
############################################################### Subsetting ###############################################################

#filter take a data frame, and made a subset based on a column
gruop_1<-dplyr::filter(pfam_matrix_clusters, Groups == "Group 1")
The_rest<-dplyr::filter(pfam_matrix_clusters, 
                        Groups == otherGroups) 
#isTRUE, evalua las dos condiciones que dan valores booleanos, si son verdad ambas condiciones, pasa a la siguiente.
uniquePFAMs <- c()
for(i in 2:(length(gruop_1)-2)){
  if(isTRUE(sum(gruop_1[,i]) > 0 & sum(The_rest[,i]) == 0)){
    #print(colnames(gruop_1[i]) )
    uniquePFAMs <- c(uniquePFAMs, colnames(gruop_1[i]))
  }
}

save(pan_matrix, clusters, pfam_matrix_clusters, gruop_1, The_rest, file="/home/mirnis/Documentos/02.CP9/pfam.RData")
load("/home/mirnis/Documentos/02.CP9/data/00.Rdata/pfam.RData")

######################################################## Getting groups ##############################################################

groupOfInterest <- c("Group 1")
allGroupNames <- unique(pfam_matrix_clusters$Groups)
otherGroups <- allGroupNames[which(allGroupNames  != groupOfInterest )]
 
#%in%  inside the other object

TheRest2 <- dplyr::filter(pfam_matrix_clusters, Groups %in%  otherGroups)

groupOfInterest %in% allGroupNames
groupOfInterest %in% otherGroups
allGroupNames %in% groupOfInterest
