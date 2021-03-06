### Principal Component Analysis (PCA) and Canonical Correlation Analysis (CCA)

* 1pseudobulk.R: Code for generating pseudobulk gene expression data.

* 2-1pca_concat_features.R: Code for running PCA and CCA by concatenating cluster-specific highly variable genes (HVGs). Global analysis (Fig 1e, 1f and Extended Data Fig 4a) and CD8 overall analysis (Extended Data Fig18c) was based on this code.

* 2-2pca_individual_celltype.R: Code for PCA and CCA within a cell cluster. A 'MANA.combined' cluster was generated by combining combining 4 MANA enriched cluster highlighted in Extended Data Fig1a. The combined MANA-enriched cluster along with the remaining 10 non-MANA enriched clusters were evaluated separately in Extended Data Fig18b and d.

* 3permute_cancor_pvalue.R: Generate permutation-based p-value by permuting sample label 10000 times.
