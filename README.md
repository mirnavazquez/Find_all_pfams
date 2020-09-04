Identify all possible PFAMs in your genome and extract usefull information
================

##### Looking for the PFAM domains.

  - Install [MEBs](https://github.com/valdeanda/mebs).

  - Get into your **MEBs** directory

  - Within the MEBs directory, create a directory call **pfam** within
    the cycles dir.

  - Move the following files to your directory cycles/pfam:
    [**entropies.tab**](https://github.com/mirnavazquez/Find_all_pfams/blob/master/data/),
    **my\_Pfam.pfam.hmm** and
    [**pfam2kegg.tab**](https://github.com/mirnavazquez/Find_all_pfams/blob/master/data/).

  - Move into the config directory and add the next line into the end of
    the **config.txt**.

<!-- end list -->

``` bash
pfam    cycles/pfam/    cycles/pfam/pfam2kegg.tab   1   1   1   1   1   1   1   1
```

  - Run MEBs.

<!-- end list -->

``` bash
perl mebs.pl -input /path/to/the/genomes/ -type genomic -comp > out_file.tsv
```

  - Explore the output file.
  - Make sure that all the PFAMs where calculated correctly for all the
    genomes. To do that we can check that all the genomes got an **ok**
    message. In the case a genome is missing you may want to run again
    MEBs.

<!-- end list -->

``` bash
tail /path/to/the/genomes/*pfam.hmmsearch.tab | grep "ok" | wc -l
```

##### Make the a clustering analysis

``` bash
python mebs/scripts/groups_1_mod.py out_file.tsv
```

##### Which PFAMs are present in a certain cluster but not in the rest?

We read the \*.tsv file within R. As well, we load a file that contains
your genome ID and to which group/cluster it corresponds.

``` r
library(vroom)
pfam_matrix<-vroom("out_pfam_2.tsv", delim="\t")
```

``` r
clusters<-vroom("PFAM_clusters", delim="\t")
```

``` r
head(clusters)
```

    ##           Genomes  Groups
    ## 1     AB_03_Bin_3 Group 3
    ## 2 AB_3033_Bin_153 Group 3
    ## 3 AB_3033_Bin_163 Group 3
    ## 4 AB_3033_Bin_184 Group 3
    ## 5  AB_3033_Bin_57 Group 5
    ## 6  AB_3033_Bin_59 Group 3

``` r
pfam_matrix_clusters<-dplyr::right_join(pfam_matrix, clusters, by = "Genomes")
```

Now run the function **get\_unique\_pfams**. This function that takes
the PFAM data frame and parse the information based on any user group.
It sums all the values of each column for a given group and the rest. If
the group’s sum is greater than 0, and if the sum of the rest is equal
to 0, we assume that the PFAM in that column is unique for that group.

``` r
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
```

We use that function to extract the unique PFAMs for each group.

``` r
Pfam_group1<-get_unique_pfams(pfam_matrix_clusters, "Group 1")
Pfam_group2<-get_unique_pfams(pfam_matrix_clusters, "Group 2")
```

##### What does this PFAM groups mean in terms of function?

First, I will map the IDs of raw names to real PFAM names.

``` r
pfam_mapping<-vroom::vroom("pfam_mapping.txt")
Project_pfams<-as.data.frame(colnames(pfam_matrix))
New_Project_pfams<-as.data.frame(Project_pfams[-1,])
colnames(New_CP9_pfams) <- c("pfam_PATHWAY")
New_Project_pfams$pfam_PATHWAY<-as.character(New_Project_pfams$pfam_PATHWAY)
Mapped_pfams<-dplyr::left_join(New_Project_pfams, pfam_mapping, by = "pfam_PATHWAY")
write.table(Mapped_pfams, file="Mapped_pfams.txt", quote = F, col.names = T, row.names = F, sep = "\t")
```

``` r
head(Mapped_pfams)
```

To get more informartion of the PFAMs, I ran this script:

``` bash
#!/bin/bash
#Source:
#https://github.com/Cantalapiedra/pfam_terms
#Author Carlos Cantalapiedra:
# code in bash
# pfam_terms.tab contains a list of PFAM identifiers

cat pfam_names.txt | while read  pfam; do
desc=$(curl http://pfam.xfam.org/family/"$pfam"/desc | head -1);
printf "$pfam\t";
printf "$desc\n";
done 2> /dev/null \
> pfam_terms.desc.tab

# Postprocessing of not found terms (since HTTP request returns always 200, even when the PFAM term was not found and an error is reported in HTML)

cat pfam_terms.desc.tab | sed 's#<\!DOCTYPE.*#NF#' > tmp && mv tmp pfam_terms.desc.tab

# END
```

##### Now the enrichement

``` r
library(dcGOR)
library(tibble)
```

``` r
get_input_for_enrch<-function(group_list, mapping_file){
  group_list_mapped<-as_tibble(group_list)
  colnames(group_list_mapped) <- c("pfam_PATHWAY")
  Mapped_pfams<-dplyr::left_join(group_list_mapped, mapping_file, by = "pfam_PATHWAY")
  return(Mapped_pfams)
  }
```

``` r
Enriched1<-get_input_for_enrch(Pfam_group1, pfam_mapping)
Enriched2<-get_input_for_enrch(Pfam_group2, pfam_mapping)
```

``` r
Pfam <- dcRDataLoader('Pfam')
```

``` r
eoutput1 <-dcEnrichment(Enriched1$PFAM, background = Mapped_CP9_pfams$PFAM, domain = 
"Pfam", ontology="GOMF")
eoutput2 <-dcEnrichment(Enriched2$PFAM, background = Mapped_CP9_pfams$PFAM, domain = 
"Pfam", ontology="GOMF")
```

