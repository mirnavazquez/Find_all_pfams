Identify all possible PFAMs in your genome and extract usefull
information
================

### Looking for the PFAM domains.

  - Install [MEBs](https://github.com/valdeanda/mebs).

  - Get into your **MEBs** directory

  - Within the MEBs directory, create a directory call **pfam** within
    the cycles dir.

  - Move the following files to your directory cycles/pfam:
    **entropies.tab**, **my\_Pfam.pfam.hmm** and **pfam2kegg.tab**.

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
