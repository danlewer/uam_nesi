# Analysis of demographic changes among people who inject illicit drugs in the UK

## Repository contents
1. R code for processing raw [UAM](https://github.com/danlewer/uam_nesi/blob/main/r_code/process_raw_uam.R) and [NESI](https://github.com/danlewer/uam_nesi/blob/main/r_code/process_raw_nesi.R) data. These scripts require individual-level data that is not available in this repository to ensure participants' confidentiality.
2. Processed non-identifiable summary data:
    * [Histograms](https://github.com/danlewer/uam_nesi/tree/main/histogram) of duration of injecting by survey year.
    * [Quantiles of age, duration, and age at initiation](https://github.com/danlewer/uam_nesi/tree/main/quantiles), both observed and modelled.
3. R code for processing summary tables (these scripts read data directly from Github):
    * [Analysis of distribution of age, duration, and age at initiation](https://github.com/danlewer/uam_nesi/blob/main/r_code/plots_and_summaries.R)
    * [Modelled estimates of numbers starting injecting each year](https://github.com/danlewer/uam_nesi/blob/main/r_code/model_new_injectors.R)
