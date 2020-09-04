Identify all possible PFAMs in your genome and extract usefull
information
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

``` r
library(vroom)
pfam_matrix<-vroom("out_pfam_2.tsv", delim="\t")
clusters<-vroom("PFAM_clusters", delim="\t")
pfam_matrix_clusters<-dplyr::right_join(pfam_matrix, clusters, by = "Genomes")
```

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
