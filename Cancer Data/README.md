This directory is designated for storing the cancer data files required for running GenMorw, and the directory format must adhere to the following specifications:

The storage directory for breast invasive carcinoma (BRCA) is as follows:

    -Cancer Data
	    -GDC_BRCA
		    -TCGA-BRCA.gistic.tsv.gz
		    -TCGA-BRCA.htseq_fpkm.tsv.gz
		    -TCGA-BRCA.methylation450.tsv.gz
		    -TCGA-BRCA.mirna.tsv.gz
		    -TCGA-BRCA.mutect2_snv.tsv.gz
The directory structure for multiple cancers is as follows:

    -Cancer Data
	    >GDC_BLCA
	    >GDC_BRCA
	    >...
