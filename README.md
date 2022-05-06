# Analysis of demographic changes among people who inject illicit drugs in the UK

## Introduction
This repository includes code and data associated with the following study:

Lewer D, Croxford S, Desai M, Emanuel E, Hope VD, McAuley A, Phipps E, Tweed EJ. The characteristics of people who inject drugs in the United Kingdom: changes in age, duration, and incidence of injecting, 1980–2019, using evidence from repeated cross-sectional surveys. Addiction. 2022. https://doi.org/10.1111/add.15911

## Repository contents
1. R code for processing raw [UAM](https://github.com/danlewer/uam_nesi/blob/main/r_code/process_raw_uam.R) and [NESI](https://github.com/danlewer/uam_nesi/blob/main/r_code/process_raw_nesi.R) data. These scripts require individual-level data that is not available in this repository to ensure participants' confidentiality.
2. Processed non-identifiable summary data:
    * [Histograms](https://github.com/danlewer/uam_nesi/tree/main/histogram) of duration of injecting by survey year.
    * [Quantiles of age, duration, and age at initiation](https://github.com/danlewer/uam_nesi/tree/main/quantiles), both observed and modelled.
3. R code for processing summary tables (these scripts read data directly from Github):
    * [Analysis of distribution of age, duration, and age at initiation](https://github.com/danlewer/uam_nesi/blob/main/r_code/plots_and_summaries.R)
    * [Modelled estimates of numbers starting injecting each year](https://github.com/danlewer/uam_nesi/blob/main/r_code/model_new_injectors.R)
    * [Additional sensitivity analyses](https://github.com/danlewer/uam_nesi/blob/main/r_code/additional_sensitivity_analyses.R) that were done during peer review but not included in the published article.

## Abstract

**Background and aims**. Mortality and drug treatment data suggest that people who inject drugs are getting older. We aimed to describe changes in the characteristics of people injecting drugs in the United Kingdom (UK).  

**Design**. Repeat cross-sectional surveys and modelling.  

**Setting**. Low-threshold services such as needle and syringe programmes.  

**Participants**. 79,900 people who recently injected psychoactive drugs in the UK, recruited as part of the Unlinked Anonymous Monitoring Survey (England, Wales, Northern Ireland, 1990-2019) and Needle Exchange Surveillance Initiative (Scotland, 2008-2019).  

**Measurements**. Age of current injectors, age at first injection, duration of injecting (each 1990-2019), and estimates of new people who started injecting (1980-2019). 

**Findings**. In England, Wales, and Northern Ireland between 1990 and 2019, the median age of current injectors increased from 27 (IQR 24-31) to 40 (IQR 34-46); median age at first injection increased from 22 (IQR 19-25) to 33 (IQR 28-39); and median years of injecting increased from 7 (IQR 3-11) to 18 (IQR 9-23). Values in Scotland were similar after 2008. The estimated number that started injecting annually in England increased from 5,470 (95% prediction interval 3,120-6,940) in 1980 to a peak of 10,270 (95% PrI 8,980-12,780) in 1998, and then decreased to 2,420 (95% PrI 1,320-5,580) in 2019. The number in Scotland followed a similar pattern, increasing from 1,220 (95% PrI 740-2,430) in 1980 to a peak of 3,080 (95% PrI 2,160-3,350) in 1998, then decreased to a 270 (95% PrI 130-600) in 2018. The timing of the peak differed between regions, with earlier peaks in London and the North West of England. 

**Conclusions**: In the UK, large cohorts started injecting drugs in the 1980s and 1990s and many still inject today. Relatively few people started in more recent years. This has led to changes in the population, including an older average age and longer injecting histories.
