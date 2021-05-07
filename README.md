## Retrieving Human Variant Alleles from Ensembl Using API by Mira Sohn

### Aims

#### - This workflow is aimed at retrieving human variant alleles from Ensembl by using the [Ensembl REST API (GET overlap/id/:id)](https://rest.ensembl.org/documentation/info/overlap_id)

#### - Reference: [The Ensembl Variant Effect Predictor (McLaren et al., 2016)](https://pubmed.ncbi.nlm.nih.gov/27268795)

### Conda environment

#### - This analysis was performed under [conda](https://conda.io/projects/conda/en/latest/index.html) environment (see the [recipe](https://github.com/Mira0507/TRF_demo/blob/master/conda_r.yml))

### Loading packages

    library(tidyverse)
    library(httr) 
    library(jsonlite) 

### Importing gene ids

#### - This demonstration uses ensembl human gene id

    # Importing my genes to a vector 
    geneid <- scan("geneid.txt", character(), quote="")

    # Exploring the imported object 
    geneid[1:10]

    ##  [1] "ENSG00000125148" "ENSG00000101361" "ENSG00000162676" "ENSG00000164687"
    ##  [5] "ENSG00000127951" "ENSG00000108691" "ENSG00000171791" "ENSG00000079459"
    ##  [9] "ENSG00000069011" "ENSG00000130164"

    class(geneid)

    ## [1] "character"

    length(geneid)

    ## [1] 11

### Retrieving data

    # Assign the server for API 
    server <- "https://rest.ensembl.org"


    variants.df <- data.frame()

    for (i in 1:length(geneid)) {

        # Retrieve data using ensembl gene ids 
        ext <- paste0("/overlap/id/",
                      geneid[i],
                      "?feature=variation")

        r <- GET(paste(server, ext, sep = ""), 
                 content_type("application/json"))
         
        stop_for_status(r)

        # Extract the variant ids 
        df <- data.frame(t(sapply(content(r),c))) %>% 
            gather("Variable", "Variant") %>% 
            dplyr::filter(str_detect(Variant, "rs")) %>%
            mutate(GENEID=geneid[i]) %>%
            dplyr::select(-Variable) 

        # Combine the retrieved variants to variants.df by rbinding
        variants.df <- rbind(variants.df, df)

    }

    # Ensure that any duplicated rows are removed
    variants.df <- variants.df[!duplicated(variants.df),] %>%
        mutate(Variant=unlist(Variant))

    # Explore the output data frame
    head(variants.df)

    ##        Variant          GENEID
    ## 1  rs570458452 ENSG00000125148
    ## 2 rs1195013231 ENSG00000125148
    ## 3  rs779703210 ENSG00000125148
    ## 4 rs1597058734 ENSG00000125148
    ## 5  rs777588985 ENSG00000125148
    ## 6  rs761890309 ENSG00000125148

    dim(variants.df)

    ## [1] 93759     2

    glimpse(variants.df)

    ## Rows: 93,759
    ## Columns: 2
    ## $ Variant <chr> "rs570458452", "rs1195013231", "rs779703210", "rs1597058734", …
    ## $ GENEID  <chr> "ENSG00000125148", "ENSG00000125148", "ENSG00000125148", "ENSG…

    # Count the number of variant alleles per gene id 
    variants.sum <- variants.df %>%
        group_by(GENEID) %>%
        summarize(Counts=n())

    # Explore the output data frame
    head(variants.sum)

    ## # A tibble: 6 x 2
    ##   GENEID          Counts
    ##   <chr>            <int>
    ## 1 ENSG00000069011   2265
    ## 2 ENSG00000079459  20404
    ## 3 ENSG00000101361   2288
    ## 4 ENSG00000108691    514
    ## 5 ENSG00000125148    294
    ## 6 ENSG00000127951   1403

    dim(variants.sum)

    ## [1] 11  2

    glimpse(variants.sum)

    ## Rows: 11
    ## Columns: 2
    ## $ GENEID <chr> "ENSG00000069011", "ENSG00000079459", "ENSG00000101361", "ENSG0…
    ## $ Counts <int> 2265, 20404, 2288, 514, 294, 1403, 13418, 3628, 1103, 45343, 30…

    # Save the output as csv files
    write.csv(variants.df, "total_variants.csv")
    write.csv(variants.sum, "counted_variants.csv")

### Session Info

    sessionInfo()

    ## R version 4.0.3 (2020-10-10)
    ## Platform: x86_64-conda-linux-gnu (64-bit)
    ## Running under: Ubuntu 20.04.2 LTS
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /home/mira/miniconda3/envs/snakemake_r/lib/libopenblasp-r0.3.15.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] jsonlite_1.7.2  httr_1.4.2      forcats_0.5.1   stringr_1.4.0  
    ##  [5] dplyr_1.0.5     purrr_0.3.4     readr_1.4.0     tidyr_1.1.3    
    ##  [9] tibble_3.1.1    ggplot2_3.3.3   tidyverse_1.3.1
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.1.1  xfun_0.20         haven_2.4.1       colorspace_2.0-1 
    ##  [5] vctrs_0.3.8       generics_0.1.0    htmltools_0.5.1.1 yaml_2.2.1       
    ##  [9] utf8_1.2.1        rlang_0.4.11      pillar_1.6.0      glue_1.4.2       
    ## [13] withr_2.4.2       DBI_1.1.1         dbplyr_2.1.1      modelr_0.1.8     
    ## [17] readxl_1.3.1      lifecycle_1.0.0   munsell_0.5.0     gtable_0.3.0     
    ## [21] cellranger_1.1.0  rvest_1.0.0       evaluate_0.14     knitr_1.31       
    ## [25] ps_1.6.0          curl_4.3.1        fansi_0.4.2       broom_0.7.6      
    ## [29] Rcpp_1.0.6        scales_1.1.1      backports_1.2.1   fs_1.5.0         
    ## [33] hms_1.0.0         digest_0.6.27     stringi_1.5.3     grid_4.0.3       
    ## [37] cli_2.5.0         tools_4.0.3       magrittr_2.0.1    crayon_1.4.1     
    ## [41] pkgconfig_2.0.3   ellipsis_0.3.2    xml2_1.3.2        reprex_2.0.0     
    ## [45] lubridate_1.7.10  assertthat_0.2.1  rmarkdown_2.7     rstudioapi_0.13  
    ## [49] R6_2.5.0          compiler_4.0.3
