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

    ##  [1] "ENSG00000101846" "ENSG00000167552" "ENSG00000005075" "ENSG00000228451"
    ##  [5] "ENSG00000260747" "ENSG00000271425" "ENSG00000204389" "ENSG00000242960"
    ##  [9] "ENSG00000263001" "ENSG00000274840"

    class(geneid)

    ## [1] "character"

    length(geneid)

    ## [1] 50

### Retrieving data

    # Assign the server for API 
    server <- "https://rest.ensembl.org"


    # Create an empty vector storing the number of variant alleles per gene id 
    variants.df <- data.frame()


    for (i in 1:length(geneid)) {

        # Retrieve data using ensembl gene ids 
        ext <- paste0("/overlap/id/",
                      geneid[i],
                      "?feature=variation")

        r <- GET(paste(server, ext, sep = ""), content_type("application/json"))

         
        stop_for_status(r)

        # Extract the variant ids 
        df <- data.frame(t(sapply(content(r),c))) %>% 
            gather("Variable", "Value") 

        # If no variant is found:
        if (ncol(df) == 0) {

            # gene id with no variant is skipped
            df.filtered <- data.frame(GENEID=c(), 
                                      Variants=c()) 

            print(paste0("Number of variant alleles in ", 
                         geneid[i], 
                         ": 0")) 

        } else {

            # Save variant ids filtered as a data frame
            df <- df %>% 
                dplyr::filter(str_detect(Value, "rs")) %>% 
                mutate(Value=unlist(Value)) 

            # Create a vector storing unique variant ids
            var <- unique(df$Value)

            # Create a data frame storing gene id and corresponding unique variant ids 
            df.filtered <- data.frame(GENEID=rep(geneid[i], 
                                                 length(var)),
                                      Variants=var)

            print(paste0("Number of variant alleles in ", 
                         geneid[i], 
                         ": ", 
                         length(var)))

        }

        # Combine the variant allele data frame to the variant.df data frame
        variants.df <- rbind(variants.df, df.filtered)


    }

    ## [1] "Number of variant alleles in ENSG00000101846: 108520"
    ## [1] "Number of variant alleles in ENSG00000167552: 1195"
    ## [1] "Number of variant alleles in ENSG00000005075: 1921"
    ## [1] "Number of variant alleles in ENSG00000228451: 1697"
    ## [1] "Number of variant alleles in ENSG00000260747: 128"
    ## [1] "Number of variant alleles in ENSG00000271425: 42095"
    ## [1] "Number of variant alleles in ENSG00000204389: 590"
    ## [1] "Number of variant alleles in ENSG00000242960: 143"
    ## [1] "Number of variant alleles in ENSG00000263001: 24187"
    ## [1] "Number of variant alleles in ENSG00000274840: 6501"
    ## [1] "Number of variant alleles in ENSG00000219507: 123"
    ## [1] "Number of variant alleles in ENSG00000214391: 472"
    ## [1] "Number of variant alleles in ENSG00000255303: 2476"
    ## [1] "Number of variant alleles in ENSG00000267736: 158"
    ## [1] "Number of variant alleles in ENSG00000140839: 4279"
    ## [1] "Number of variant alleles in ENSG00000228532: 70"
    ## [1] "Number of variant alleles in ENSG00000234545: 6937"
    ## [1] "Number of variant alleles in ENSG00000244161: 1892"
    ## [1] "Number of variant alleles in ENSG00000215784: 2860"
    ## [1] "Number of variant alleles in ENSG00000140459: 6816"
    ## [1] "Number of variant alleles in ENSG00000271254: 0"
    ## [1] "Number of variant alleles in ENSG00000225813: 413"
    ## [1] "Number of variant alleles in ENSG00000163283: 1634"
    ## [1] "Number of variant alleles in ENSG00000174469: 578486"
    ## [1] "Number of variant alleles in ENSG00000254635: 2798"
    ## [1] "Number of variant alleles in ENSG00000235314: 842"
    ## [1] "Number of variant alleles in ENSG00000262188: 1756"
    ## [1] "Number of variant alleles in ENSG00000214425: 9159"
    ## [1] "Number of variant alleles in ENSG00000268916: 576"
    ## [1] "Number of variant alleles in ENSG00000237510: 1238"
    ## [1] "Number of variant alleles in ENSG00000196139: 20709"
    ## [1] "Number of variant alleles in ENSG00000273472: 344"
    ## [1] "Number of variant alleles in ENSG00000165949: 2856"
    ## [1] "Number of variant alleles in ENSG00000235559: 83"
    ## [1] "Number of variant alleles in ENSG00000214518: 242"
    ## [1] "Number of variant alleles in ENSG00000228218: 250"
    ## [1] "Number of variant alleles in ENSG00000233476: 336"
    ## [1] "Number of variant alleles in ENSG00000223722: 101"
    ## [1] "Number of variant alleles in ENSG00000177855: 131"
    ## [1] "Number of variant alleles in ENSG00000224546: 447"
    ## [1] "Number of variant alleles in ENSG00000259341: 2929"
    ## [1] "Number of variant alleles in ENSG00000272189: 409"
    ## [1] "Number of variant alleles in ENSG00000284727: 21139"
    ## [1] "Number of variant alleles in ENSG00000037042: 2093"
    ## [1] "Number of variant alleles in ENSG00000153902: 5084"
    ## [1] "Number of variant alleles in ENSG00000276550: 19162"
    ## [1] "Number of variant alleles in ENSG00000146707: 5179"
    ## [1] "Number of variant alleles in ENSG00000166268: 31376"
    ## [1] "Number of variant alleles in ENSG00000243678: 1642"
    ## [1] "Number of variant alleles in ENSG00000266869: 18557"

    # Ensure that any duplicated rows are removed
    variants.df <- variants.df[!duplicated(variants.df),] 

    # Explore the output data frame
    head(variants.df)

    ##            GENEID     Variants
    ## 1 ENSG00000101846 rs1602854972
    ## 2 ENSG00000101846 rs1422412432
    ## 3 ENSG00000101846 rs1015958816
    ## 4 ENSG00000101846  rs935459767
    ## 5 ENSG00000101846 rs1490518983
    ## 6 ENSG00000101846 rs1268004492

    dim(variants.df)

    ## [1] 943031      2

    glimpse(variants.df)

    ## Rows: 943,031
    ## Columns: 2
    ## $ GENEID   <chr> "ENSG00000101846", "ENSG00000101846", "ENSG00000101846", "ENS…
    ## $ Variants <chr> "rs1602854972", "rs1422412432", "rs1015958816", "rs935459767"…

    # Count the number of variant alleles per gene id 
    variants.sum <- variants.df %>%
        group_by(GENEID) %>%
        summarize(Counts=n()) 

    # Explore the output data frame
    head(variants.sum)

    ## # A tibble: 6 x 2
    ##   GENEID          Counts
    ##   <chr>            <int>
    ## 1 ENSG00000005075   1921
    ## 2 ENSG00000037042   2093
    ## 3 ENSG00000101846 108520
    ## 4 ENSG00000140459   6816
    ## 5 ENSG00000140839   4279
    ## 6 ENSG00000146707   5179

    dim(variants.sum)

    ## [1] 49  2

    glimpse(variants.sum)

    ## Rows: 49
    ## Columns: 2
    ## $ GENEID <chr> "ENSG00000005075", "ENSG00000037042", "ENSG00000101846", "ENSG0…
    ## $ Counts <int> 1921, 2093, 108520, 6816, 4279, 5179, 5084, 1634, 2856, 31376, …

    # Save the output as csv files
    write.csv(variants.df, "total_variants.csv")
    write.csv(variants.sum, "counted_variants.csv")


    ## NOTE: The genes missing in the csv files have zero variant!! 

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
    ##  [5] dplyr_1.0.6     purrr_0.3.4     readr_1.4.0     tidyr_1.1.3    
    ##  [9] tibble_3.1.1    ggplot2_3.3.3   tidyverse_1.3.1
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] bslib_0.2.4       tidyselect_1.1.1  xfun_0.20         haven_2.4.1      
    ##  [5] colorspace_2.0-1  vctrs_0.3.8       generics_0.1.0    htmltools_0.5.1.1
    ##  [9] yaml_2.2.1        utf8_1.2.1        rlang_0.4.11      jquerylib_0.1.4  
    ## [13] pillar_1.6.0      glue_1.4.2        withr_2.4.2       DBI_1.1.1        
    ## [17] dbplyr_2.1.1      modelr_0.1.8      readxl_1.3.1      lifecycle_1.0.0  
    ## [21] munsell_0.5.0     gtable_0.3.0      cellranger_1.1.0  rvest_1.0.0      
    ## [25] evaluate_0.14     knitr_1.31        ps_1.6.0          curl_4.3.1       
    ## [29] fansi_0.4.2       broom_0.7.6       Rcpp_1.0.6        scales_1.1.1     
    ## [33] backports_1.2.1   fs_1.5.0          hms_1.0.0         digest_0.6.27    
    ## [37] stringi_1.5.3     grid_4.0.3        cli_2.5.0         tools_4.0.3      
    ## [41] sass_0.3.1        magrittr_2.0.1    crayon_1.4.1      pkgconfig_2.0.3  
    ## [45] ellipsis_0.3.2    xml2_1.3.2        reprex_2.0.0      lubridate_1.7.10 
    ## [49] assertthat_0.2.1  rmarkdown_2.7     rstudioapi_0.13   R6_2.5.0         
    ## [53] compiler_4.0.3
