###################################
Pairwise disease similarity: pairdc
###################################

This program creates input data for HN2V.
It creates two files:

1. g2d_associations_test_7_2017.tsv
This file contains only gene/disease links from date 2 and should be used for testing
7_2017 is "Date 2", i.e., the date for which we want to predict gene/disease links.


2. "pairwise_disease_similarity.tsv"
This file contains disease/disease links from date 2 as well as gene/disease links from date 1
and should be used for training. Not that the first file only contains the "new" gene/disease
links from date2 and does not have the "old" links from date 1.

To run the app. ::

    java -jar pairdc.jar --hpo hp.obo \\
        -a phenotype.hpoa \\
        --geneinfo Homo_sapiens_gene_info.gz \\
        --mim2genemedgen mim2gene_medgen \\
        --date1 6/2014 \\
        --date2 7/2017

Note that it is easiest to use the LIRICAL download --overwrite command to get
new versions of the input files (they get written to the data directory in
the main LIRICAL directory).