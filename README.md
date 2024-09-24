# GenMorw
GenMorw is a computational tool designed for identifying **cancer driver genes** by prioritizing genes within **single patient**, **cancer cohort** and **pancancer** levels. It utilizes essential Single Nucleotide Variants (SNV) data, along with optional Copy Number Variants (CNV), DNA methylation data, gene expression, and miRNA expression data to construct extensive gene-patient heterogeneity networks. These networks are then analyzed using a Random Walk with Restart model to predict relationships between genes and patients in the network. The integration of diverse data types ensures a comprehensive approach to driver gene identification, leveraging the complex interactions within the cellular environment.
![flowchart.png](https://github.com/JLiuing/GenMorw/blob/master/flowchart.png)

###  GenMorw overview
GenMorw segments patients into distinct groups based on features derived from multiple data types, constructing heterogeneous networks that connect these patient groups with their mutated genes using various PPI networks. By employing the Random Walk with Restart (RWR) model, patient-gene associations are captured, allowing GenMorw to integrate driver gene predictions across multiple datasets. The final prioritization of mutated genes is achieved at both the individual and cohort levels through a cumulative rank voting method. Additionally, GenMorw offers extended functionalities, including the **construction of the GenMorw-network**, **model interpretation**, and **the assessment and validation of potential gene-drug interactions**. 

## How to use GenMorw
### Required environment configuration
GenMorw is designed to operate within a **Windows 11/10** environment and leverages multi-threading to enhance processing efficiency. It is recommended to use a computing setup equipped with a minimum of **10 cores** and **40 GB** of RAM to handle the intensive computational tasks effectively. This specification ensures optimal performance and stability during the execution of the algorithm's complex processes. Below are the **R packages** required to run the algorithm:

 - R ( >= 4.2.0, < 4.3.0)
 - BiocManager ( = 1.30.20)
 - tidyverse ( = 2.0.0)
 - data.table ( = 1.14.8)
 - Seurat ( = 4.3.0)
 - Matrix ( = 1.6-1.1)
 - readr ( = 2.4.1)
 - foreach ( = 1.5.2)
 - doParallel ( = 1.0.17)
 ### Optional packages for additional function
The packages mentioned above are essential for the process of predicting driver genes. Additional packages are required for extended functionalities such as the construction of the GenMorw-network, identification of its strongly connected components and clique structures, model interpretation, and the assessment and validation of potential gene-drug interactions.
 - rentrez ( = 1.2.3)
 - openxlsx ( = 4.2.5.2)
 - ggplot2 ( = 3.5.0)
 - ggrepel ( = 0.9.5)
 - PharmacoGx ( = 3.2.0)
 - igraph ( = 2.0.3)
 ### 1. Prioritization of predicted individual patient and cancer cohort driver gene
Once the environment is properly configured, users need to input various cancer data types to complete the prediction of driver genes. GenMorw accepts inputs in the format used by [XenaBrowser](https://xenabrowser.net/datapages/) for GDC TCGA data. Specifically, it requires gene-level Copy Number Variants (CNV), DNA methylation data, MuTect2 Variant Aggregation for somatic mutations, gene expression data (FPKM), and miRNA expression data. Data should be organized in the directory structure `Cancer Data/GDC_[cancer]/...` and stored within the `Cancer Data` folder for processing. The storage directory for breast invasive carcinoma (BRCA) is as follows:

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
GenMorw is capable of predicting driver genes for 33 GDC TCGA or more cancer types, provided the data is formatted consistently. The `main.R` program offers methods for predicting driver gene prioritization both cancer cohorts and individual patients. After placing the desired cancer data into the `Cancer Data` folder, the `main.R` program will automatically recognize (or user can manually specify) all cancers and sequentially predict driver genes for each.

Additionally, GenMorw comes with 12 predefined Protein-Protein Interaction (PPI) networks used for constructing heterogeneous networks, which are stored in the `Data/network_rds` directory. Users can also add their own networks to this folder. The `main.R` program utilizes all the contents in the `Data/network_rds` directory to construct heterogeneous networks, or allows for manual definition of these networks. The format for these network files is as follows:

		TP53		PIK3CA		TTN	...
	TP53	0		1		1
	PIK3CA	1		0		0
	TTN	1		0		0
	...			
After the program completes execution, all results are stored in the `Results/Data/[cancer]` directory. For example, the prioritization of driver genes for BRCA is stored in `Results/Data/BRCA/final_symbol_sort_patientFirst.txt`, and the prioritization of driver genes for individual cancer patients within BRCA is stored in `Results/Data/BRCA/single_patients_pred.rds`.
 ### 2. Prioritization of predicted pancancer driver gene
Before executing this functionality, it is necessary to complete the driver gene prediction for at least two different cancers. GenMorw is capable of performing pancancer analysis for multiple cancer types, not limited to just 33. Specifically, the `pancancer.R` script integrates the results from all cancers in the `Results/Data` directory and consolidates them into a pancancer prioritization of driver genes. Alternatively, users can freely specify combinations of cancers they are interested in analyzing. The prediction results are stored in the file `Results/Data/pancancer_prediction.txt`.
 ### 3. Construction of the GenMorw-network
During execution, GenMorw thoroughly explores the interactions between genes, incorporating multiple types of information such as gene mutations, methylation, gene and miRNA expression, and known gene interaction relationships. Before utilizing this functionality, the `retain_mmPPI` parameter in the `main.R` program should be set to `True`. This setting enables the program to automatically retain the interactions in the stabilized heterogeneous network, which will be stored in the `Models/Data/[cancer]/mmPPI` folder. It is also advisable to allocate **40-80 GB** of disk space for each cancer analysis.
After completing the prediction for an individual cancer, users can utilize the `GenMorw-network construction.R` script to specify which cancer's GenMorw-network to construct. This network consolidates gene interaction information from all heterogeneous networks while retaining high-confidence nodes and edges. Research has shown that the strongly connected components and cliques within this network capture many predicted and known driver genes. These strongly connected components and cliques are stored in the `Results/Data/[cancer]/[cancer]_GenMorw-network_strongly_connected_components.rds` and `Results/Data/[cancer]/[cancer]_GenMorw-network_cliques.rds`, respectively.
 ### 4. Model interpretation
In the pan-cancer prioritization of driver genes, the top-ranked genes typically command significant attention. Our primary interest lies in discerning the specific cancers that propel these genes to the forefront of the rankings, as well as elucidating the degree of contribution each cancer exerts on these premier genes. This understanding will provide critical insights into which cancers merit intensified scrutiny in relation to the top-ranked driver genes.
The `Explain pancancer.R` script calculates the contribution percentage of each  cancer under current study to the every gene of pancancer prioritization. The results are stored in `Results/Data/Cancer_contribution.rds`.
Similarly, the `Explain single cancer.R` script computes the contribution percentage of different PPI networks to the prioritization of the cancer cohort during the algorithm's execution. The results are stored in `Results/Data/[cancer]/Explaination_PPI_to_gene.rds`.
 ### 5. Estimate and validation of potential gene-drug interactions
Another supplementary function of GenMorw is its capability to assess the potential relationships and significance between the predicted driver genes and 202 different drugs, utilizing drug sensitivity data from the [GDSC2](https://zenodo.org/records/5787145) and [CCLE](https://zenodo.org/records/3905462) datasets provided by `PharmacoGx`. The `Anti-cancer drug identification.R` script requires users to download and specify the paths to these datasets to complete the analysis. Subsequently, it automatically searches for studies based on `"gene symbol AND drug"` across all databases on [Entrez](https://www.ncbi.nlm.nih.gov/), recording these as evidence. The results are saved in `Results/Data/[cancer]/GDSC2_drug_record_with_evidence.xlsx` and `Results/Data/[cancer]/CCLE_drug_record_with_evidence.xlsx`, along with visualizations.
