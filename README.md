# ILMS

If you have any questions, please contact [wangzr27@mail.sysu.edu.cn](mailto:wangzr27@mail.sysu.edu.cn);  lren@mail3.sysu.edu.cn; linrx6@mail2.sysu.edu.cn.



## Package Requirements

```R
###Model Construction
library(survival)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(CoxBoost)
library(survivalsvm)
library(BART)
library(compareC)

###Survival Analysis Related Packages
library(survival)
library(survminer)
library(timeROC)
library(pROC)
library(rms)
library(powerSurvEpi)

##Differential Analysis and Pathway Enrichment Packages
library(DESeq2)
library(limma)
library(GSVA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)
library(pathview)
library(topGO)

###Multi-omics Integration Analysis Packages
library(MOVICS)
library(IOBR)
library(maftools)

###Visualization Packages
library(ggbreak)
library(ggwaffle)
library(reshape)
library(ggsci)
library(ggplot2)
library(VennDiagram)
library(ConsensusClusterPlus)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(export)
library(forestplot)
library(patchwork)
library(ggExtra)

###Other Utility Packages
library(tidyverse)
library(dplyr)
library(tidyr)
library(tibble)
library(data.table)
library(miscTools)



sessionInfo()
```



```R
sessionInfo()
R version 4.2.2 (2022-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 22631)

Matrix products: default

locale:
[1] LC_COLLATE=Chinese (Simplified)_China.utf8 
[2] LC_CTYPE=Chinese (Simplified)_China.utf8   
[3] LC_MONETARY=Chinese (Simplified)_China.utf8
[4] LC_NUMERIC=C                               
[5] LC_TIME=Chinese (Simplified)_China.utf8    

attached base packages:
[1] stats4    grid      stats     graphics  grDevices utils     datasets  methods  
[9] base     

other attached packages:
 [1] maftools_2.14.0             IOBR_0.99.9                
 [3] tidyHeatmap_1.8.1           MOVICS_0.99.17             
 [5] topGO_2.50.0                GO.db_3.16.0               
 [7] graph_1.76.0                pathview_1.38.0            
 [9] DOSE_3.24.2                 enrichplot_1.18.4          
[11] org.Hs.eg.db_3.16.0         AnnotationDbi_1.60.2       
[13] clusterProfiler_4.6.2       GSVA_1.46.0                
[15] limma_3.54.2                DESeq2_1.38.3              
[17] SummarizedExperiment_1.28.0 Biobase_2.58.0             
[19] MatrixGenerics_1.10.0       matrixStats_1.0.0          
[21] GenomicRanges_1.50.2        GenomeInfoDb_1.34.9        
[23] IRanges_2.32.0              S4Vectors_0.36.2           
[25] BiocGenerics_0.44.0         miscTools_0.6-28           
[27] data.table_1.14.8           forcats_1.0.0              
[29] stringr_1.5.0               dplyr_1.1.3                
[31] purrr_1.0.2                 readr_2.1.4                
[33] tidyr_1.3.0                 tibble_3.2.1               
[35] tidyverse_1.3.2             ggExtra_0.10.0             
[37] patchwork_1.1.3             forestplot_3.1.1           
[39] abind_1.4-5                 checkmate_2.1.0            
[41] export_0.3.0                circlize_0.4.15            
[43] ComplexHeatmap_2.14.0       ConsensusClusterPlus_1.62.0
[45] VennDiagram_1.7.3           futile.logger_1.4.3        
[47] ggsci_3.0.0                 reshape_0.8.9              
[49] ggwaffle_0.2.5              ggbreak_0.1.2              
[51] rms_6.3-0                   SparseM_1.81               
[53] Hmisc_4.7-1                 Formula_1.2-4              
[55] lattice_0.20-45             pROC_1.18.4                
[57] timeROC_0.4                 survminer_0.4.9            
[59] ggpubr_0.6.0                ggplot2_3.4.4              
[61] compareC_1.3.2              BART_2.9.4                 
[63] nnet_7.3-18                 nlme_3.1-160               
[65] survivalsvm_0.0.5           CoxBoost_1.5               
[67] prodlim_2023.08.28          gbm_2.1.8.1                
[69] superpc_1.12                plsRcox_1.7.7              
[71] glmnet_4.1-8                Matrix_1.6-5               
[73] randomForestSRC_3.2.2       survival_3.4-0             
[75] powerSurvEpi_0.1.3         

loaded via a namespace (and not attached):
  [1] KEGGREST_1.38.0             genefilter_1.80.3          
  [3] spam_2.9-1                  locfit_1.5-9.8             
  [5] corrplot_0.92               labelled_2.12.0            
  [7] vctrs_0.6.4                 utf8_1.2.4                 
  [9] blob_1.2.4                  risksetROC_1.0.4.1         
 [11] withr_2.5.1                 foreign_0.8-83             
 [13] gdtools_0.3.3               registry_0.5-1             
 [15] uuid_1.1-1                  readxl_1.4.1               
 [17] lifecycle_1.0.3             cellranger_1.1.0           
 [19] munsell_0.5.0               ragg_1.2.6                 
 [21] aricode_1.0.3               devEMF_4.4-1               
 [23] ScaledMatrix_1.6.0          ggalluvial_0.12.5          
 [25] codetools_0.2-18            preprocessCore_1.60.2      
 [27] rmeta_3.0                   annotate_1.76.0            
 [29] parallelly_1.36.0           fs_1.5.2                   
 [31] fastmatch_1.1-4             ellipse_0.4.3              
 [33] impute_1.72.3               stringi_1.7.8              
 [35] stargazer_5.2.3             polyclip_1.10-6            
 [37] rhdf5filters_1.10.1         yulab.utils_0.1.0          
 [39] sandwich_3.0-2              cluster_2.1.4              
 [41] ggraph_2.1.0                ape_5.7-1                  
 [43] pkgconfig_2.0.3             sparseMatrixStats_1.10.0   
 [45] googledrive_2.0.0           lubridate_1.9.3            
 [47] timechange_0.2.0            rvg_0.3.3                  
 [49] httr_1.4.7                  flextable_0.9.2            
 [51] IntNMF_1.2.0                igraph_1.5.1               
 [53] treeio_1.22.0               GetoptLong_1.0.5           
 [55] modeltools_0.2-23           beachmat_2.14.2            
 [57] mitools_2.4                 graphlayouts_1.0.1         
 [59] haven_2.5.3                 ggfun_0.1.3                
 [61] gson_0.1.0                  htmltools_0.5.7            
 [63] miniUI_0.1.1.1              viridisLite_0.4.2          
 [65] NMF_0.26                    pillar_1.9.0               
 [67] later_1.3.0                 httpcode_0.3.0             
 [69] glue_1.6.2                  DBI_1.1.3                  
 [71] BiocParallel_1.32.6         plyr_1.8.9                 
 [73] foreach_1.5.2               dotCall64_1.0-2            
 [75] timereg_2.0.4               gtable_0.3.4               
 [77] GOSemSim_2.24.0             survivalROC_1.0.3          
 [79] pls_2.8-2                   rsvd_1.0.5                 
 [81] caTools_1.18.2              GlobalOptions_0.1.2        
 [83] latticeExtra_0.6-30         fastmap_1.1.0              
 [85] broom_1.0.5                 rARPACK_0.11-0             
 [87] ggpp_0.5.4                  promises_1.2.0.1           
 [89] FNN_1.1.3.2                 textshaping_0.3.7          
 [91] hms_1.1.3                   ggforce_0.4.1              
 [93] askpass_1.2.0               png_0.1-8                  
 [95] survey_4.2-1                clue_0.3-65                
 [97] ggtree_3.6.2                tableone_0.13.2            
 [99] lazyeval_0.2.2              crayon_1.5.2               
[101] gridBase_0.4-7              reprex_2.0.2               
[103] boot_1.3-28                 tidyselect_1.2.0           
[105] xfun_0.40                   BiocSingular_1.14.0        
[107] kernlab_0.9-32              splines_4.2.2              
[109] PINSPlus_2.0.6              mogsa_1.32.0               
[111] rappdirs_0.3.3              bit64_4.0.5                
[113] rngtools_1.5.2              lambda.r_1.2.4             
[115] modelr_0.1.10               jstable_1.1.2              
[117] jpeg_0.1-10                 fields_14.1                
[119] MatrixModels_0.5-2          ggsignif_0.6.4             
[121] data.tree_1.0.0             permute_0.9-7              
[123] quantreg_5.97               htmlTable_2.4.1            
[125] pamr_1.56.1                 xtable_1.8-4               
[127] googlesheets4_1.0.1         cachem_1.0.6               
[129] DelayedArray_0.24.0         vegan_2.6-4                
[131] systemfonts_1.0.5           mime_0.12                  
[133] rjson_0.2.21                aplot_0.2.2                
[135] ggrepel_0.9.4               rstatix_0.7.2              
[137] SuppDists_1.1-9.7           numDeriv_2016.8-1.1        
[139] quadprog_1.5-8              tools_4.2.2                
[141] cli_3.6.1                   magrittr_2.0.3             
[143] rgl_0.110.2                 proxy_0.4-27               
[145] future.apply_1.11.0         ggplotify_0.1.2            
[147] DelayedMatrixStats_1.20.0   mixOmics_6.22.0            
[149] assertthat_0.2.1            officer_0.6.3              
[151] qvalue_2.30.0               fgsea_1.24.0               
[153] lpSolve_5.6.17              DNAcopy_1.72.3             
[155] sna_2.7                     HDF5Array_1.26.0           
[157] fontquiver_0.2.1            survcomp_1.48.0            
[159] mgcv_1.8-41                 survMisc_0.5.6             
[161] gfonts_0.2.0                tweenr_2.0.2               
[163] Rgraphviz_2.42.0            InterSIM_2.2.0             
[165] entropy_1.3.1               multcomp_1.4-20            
[167] zlibbioc_1.44.0             zip_2.3.0                  
[169] survAUC_1.2-0               coxme_2.2-18.1             
[171] shadowtext_0.1.2            tzdb_0.4.0                 
[173] geneplotter_1.76.0          DiagrammeR_1.0.10          
[175] fansi_1.0.5                 tidygraph_1.2.3            
[177] GSEABase_1.60.0             TH.data_1.1-1              
[179] KernSmooth_2.23-20          backports_1.4.1            
[181] XVector_0.38.0              bootstrap_2019.6           
[183] interp_1.1-3                farver_2.1.1               
[185] bit_4.0.5                   gplots_3.1.3               
[187] openxlsx_4.2.5.1            shiny_1.7.3                
[189] KMsurv_0.1-5                scatterpie_0.2.1           
[191] maps_3.4.1                  futile.options_1.0.1       
[193] dendextend_1.16.0           downloader_0.4             
[195] CIMLR_1.0.0                 KEGGgraph_1.58.3           
[197] minqa_1.2.6                 viridis_0.6.4              
[199] rstudioapi_0.14             iterators_1.0.14           
[201] Rhdf5lib_1.20.0             shape_1.4.6                
[203] gtools_3.9.4                network_1.18.0             
[205] listenv_0.9.0               reshape2_1.4.4             
[207] rhdf5_2.42.1                gargle_1.2.1               
[209] flexclust_1.4-1             generics_0.1.3             
[211] colorspace_2.1-0            oompaData_3.1.3            
[213] coca_1.1.0                  graphite_1.44.0            
[215] base64enc_0.1-3             XML_3.99-0.14              
[217] e1071_1.7-13                dbplyr_2.2.1               
[219] RColorBrewer_1.1-3          GenomeInfoDbData_1.2.9     
[221] Biostrings_2.66.0           evaluate_0.22              
[223] memoise_2.0.1               coda_0.19-4                
[225] knitr_1.44                  doParallel_1.0.17          
[227] httpuv_1.6.6                fontLiberation_0.1.0       
[229] ClassDiscovery_3.4.0        class_7.3-20               
[231] irlba_2.3.5.1               lars_1.3                   
[233] Rcpp_1.0.9                  polynom_1.4-1              
[235] openssl_2.1.1               formatR_1.14               
[237] BiocManager_1.30.22         bipartite_2.18             
[239] iClusterPlus_1.34.3         jsonlite_1.8.7             
[241] km.ci_0.5-6                 fontBitstreamVera_0.1.1    
[243] clusterRepro_0.9            RSpectra_0.16-1            
[245] CMScaller_2.0.1             digest_0.6.30              
[247] polspline_1.1.20            oompaBase_3.2.9            
[249] cowplot_1.1.1               bitops_1.0-7               
[251] RSQLite_2.3.1               crul_1.4.0                 
[253] rmarkdown_2.18              globals_0.16.2             
[255] compiler_4.2.2              zoo_1.8-12                 
[257] ridge_3.3                   limSolve_1.5.6             
[259] carData_3.0-5               pec_2022.05.04             
[261] pracma_2.4.2                gridGraphics_0.5-1         
[263] rlang_1.1.1                 nloptr_2.0.3               
[265] SingleCellExperiment_1.20.1 heatmap.plus_1.3           
[267] lava_1.7.2.1                sva_3.46.0                 
[269] bdsmatrix_1.3-6             rvest_1.0.3                
[271] svd_0.5.4.1                 visNetwork_2.1.2           
[273] future_1.33.0               mvtnorm_1.2-3              
[275] htmlwidgets_1.6.3           geepack_1.3.9              
[277] curl_4.3.3                  parallel_4.2.2             
[279] plsRglm_1.5.1               SNFtool_2.2                
[281] edgeR_3.40.2                corpcor_1.6.10             
[283] scales_1.2.1                RcppParallel_5.1.7         
[285] lme4_1.1-34                 HDO.db_0.99.1              
[287] deldir_1.0-6                gridExtra_2.3              
[289] RCurl_1.98-1.12             car_3.1-2                  
[291] MASS_7.3-58.1               ellipsis_0.3.2             
[293] tidytree_0.4.5              xml2_1.3.5                 
[295] rpart_4.1.19                R6_2.5.1                   
[297] mclust_6.0.0                statnet.common_4.7.0       
```

