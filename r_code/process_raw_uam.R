# This script processes UAM data. It requires raw UAM data.
# 1. Imputes missing age/duration, etc.
# 2. Produces the chart in supplementary material showing proportion that injected in the past year
# 3. Drops ineligible participants
# 4. Produces descriptive table
# 5. Reports quantiles of age, duration, and age at initiation
# 6. Estimates modelled quantiles
# 7. Estimates pairwise ratios for modelling of number injecting for the first time each year

setwd("H:/uam_age_trends")

# -- libraries --

library(data.table) # for reading and manipulating data
library(quantreg) # for quantile regression
library(RColorBrewer) # for plot colours
library(devEMF) # for enhanced metafile device (graphics)

# do quantile regression for interquartile range and predict confidence intervals

qm <- function(d, outcome = 'age', exposure = 'year', q = c(0.25, 0.5, 0.75), pyears = 1990:2019, sampleVal = F) {
  f <- as.formula(paste0(outcome, '~poly(', exposure, ',2)'))
  modelVals <- rbindlist(lapply(q, function(x) {
    m <- rq(f, data = d, tau = x)
    nd <- data.table(exposure = pyears)
    names(nd)[1] <- exposure
    cbind(nd, predict(m, newdata = nd, interval = 'confidence'), tau = x)
  }))
  if (sampleVal) {
    sampleVals <- d[get(exposure) %in% pyears, .(sample = quantile(get(outcome), probs = q)), by = exposure]
    sampleVals$tau <- rep(q, length(pyears))
    return(sampleVals[modelVals, on = c('tau', exposure)])
  } else {
    return(modelVals)
  }
}

#  :::::::::
#  read data
#  .........

d <- fread("uam_age_27jan2022.csv")
d[, id := .I]

#  :::::::::::::::::::::::::::::::::::::::::::::::::::
#  impute year started, age, duration, and age started
#  ...................................................

take_first <- function (...) {
  x <- cbind(...)
  i <- max.col(!is.na(x), ties.method = 'first')
  x[cbind(seq_len(length(i)), i)]
}

# using year, fyear, fage, injdur, age

# age

d[, age2 := fage + injdur]
d[, ageI := take_first(age, age2)]

# first year

d[, fyear2 := year - injdur]
d[, fyear3 := year - (age - fage)]
d[, fyearI := take_first(fyear, fyear2, fyear3)]

# duration

d[, injdur2 := year - fyearI]
d[, injdur3 := ageI - fage]
d[, injdurI := take_first(injdur, injdur2, injdur3)]

# age started

d[, fage2 := ageI - injdurI]
d[, fageI := take_first(fage, fage2)]

#  ::::::::::::::::::::::::::::::::
#  impute injected in the last year
#  ................................

d$injy1[d$injy1 == ''] <- NA
d[, injy2 := (year - lastyr < 2) | (ageI - lage < 2)]
d[, injy2 := fifelse(injy2, 'Yes', 'No')]
d[, injyI := take_first(injy1, injy2)]

d[, table(year, is.na(injyI))] # missing 1990-1992

#  plot of recent injecting
#  ........................

d[, inyM := fifelse(is.na(injyI), 'missing', injyI)]
d[, inyM := factor(inyM, c('Yes', 'No', 'missing'))]

recent_injecting <- dcast(d, year ~ inyM, value.var = 'year', fun.aggregate = length)
# fwrite(recent_injecting, 'recent_injecting.csv')

#  :::::::::::::::
#  drop ineligible
#  ...............

nrow(d) # 91181
d[!(injyI == 'Yes' | year %in% c(1990:1992)), .N] # 17014

# did not inject in the past year

d <- d[injyI == 'Yes' | year %in% c(1990:1992)]
nrow(d) # 69947

# missing age, duration, age started, or year started

d[, missing := rowSums(is.na(cbind(ageI, injdurI, fyearI, fageI))) > 0]
d[, sum(missing)] # 3274
d <- d[missing == F]
nrow(d) # 66673

# implausible values

d[(injdurI < 0 | ageI < 10 | fyearI > year), .N] # 127
d <- d[!(injdurI < 0 | ageI < 10 | fyearI > year)]
nrow(d) # 66655

#  ::::::::::::
#  missing data
#  ............

d[d == ''] <- NA

colSums(is.na(d)) # number missing
round(colSums(is.na(d)) / nrow(d) * 100, 1) # % missing

missing_data <- sapply(split.data.frame(is.na(d), f = d$year), colSums)
sample_size <- d[, .N, year]
round(t(missing_data) / sample_size[, N] * 100, 0)

#  :::::::::::::::::::::
#  table of participants
#  .....................

d[, total := 1]
year_lims <- seq(1990, 2020, 5)
d[, yearGroup := findInterval(year, year_lims)]
d[, yearGroup := factor(yearGroup, seq_along(year_lims), year_lims)]

dv <- function (v) {
  r <- d[, .N, get(v)][order(get)]
  r[, pc := N / sum(N) * 100]
  r[, pc := format(round(pc, 1), digits = 1, nsmall = 1)]
  r[, N := prettyNum(N, big.mark = ',')]
  r[, N := paste0(N, '(', pc, ')')]
  r[, N := gsub('\\(', ' (', gsub(' ', '', N))]
  r[, pc := NULL]
  names(r)[1] <- 'level'
  cbind(variable = c(v, rep('', nrow(r)-1)), r)
}

desctab <- rbindlist(lapply(c('total', 'gen', 'govoff', 'yearGroup'), dv))
fwrite(desctab, 'descriptive_table_1april2022.csv')

# proportion injecting heroin

heroin <- dcast(d, year ~ heroin1, fun.aggregate = length)
names(heroin) <- c('year', 'missing', 'no', 'yes')
heroin <- heroin[, c('year', 'yes', 'no', 'missing')]
heroin[, number_not_missing := yes + no]
heroin[, pc_heroin := yes / number_not_missing * 100]
fwrite(heroin, 'pc_using_heroin_1april2022.csv')

#  :::::
#  drugs
#  .....

table(d$heroin1, useNA = 'always')
#   No   Yes  <NA> 
#   1801 19989 44755 
prop.table(table(d$heroin1))
# No         Yes 
# 0.08265259 0.91734741

#  :::::::::::::::::::
#  age of participants
#  ...................

d_age <- d[!is.na(govoff) & govoff != 'Northern Ireland' & !is.na(age), c('govoff', 'year', 'ageI')]
d_age <- c(split(d_age, f = d_age$govoff), England = list(d_age))

# model and actual values

age_reg <- lapply(d_age, qm, outcome = 'ageI')
for (i in seq_along(age_reg)) {age_reg[[i]]$region <- names(age_reg)[i]}
age_act <- lapply(d_age, function(x) x[, .(md = quantile(ageI, probs = 0.5), lower = quantile(ageI, probs = 0.25), upper = quantile(ageI, probs = 0.75)), year])
for (i in seq_along(age_act)) {age_act[[i]]$region <- names(age_act)[i]}

age_reg <- rbindlist(age_reg)
age_act <- rbindlist(age_act)

#  :::::::::::::::::
#  age of initiation
#  .................

ndy <- 1980:2019

d_ini <- d[, c('govoff', 'fyearI', 'fageI', 'year')]
d_ini <- d_ini[fyearI %in% ndy]
d_ini <- d_ini[(year - fyearI) <= 3]
d_ini <- c(split(d_ini, f = d_ini$govoff), England = list(d_ini))

# model and actual values

ini_reg <- lapply(d_ini, qm, exposure = 'fyearI', outcome = 'fageI', pyears = 1990:2019)
for (i in seq_along(ini_reg)) {ini_reg[[i]]$region <- names(ini_reg)[i]}
ini_act <- lapply(d_ini, function(x) x[, .(md = quantile(fageI, probs = 0.5), lower = quantile(fageI, probs = 0.25), upper = quantile(fageI, probs = 0.75)), fyearI])
for (i in seq_along(ini_act)) {ini_act[[i]]$region <- names(ini_act)[i]}

ini_reg <- rbindlist(ini_reg)
ini_act <- rbindlist(ini_act)

# including all participants

ndy <- 1980:2019

d_ini2 <- d[, c('govoff', 'fyearI', 'fageI', 'year')]
d_ini2 <- d_ini2[fyearI %in% ndy]
d_ini2 <- c(split(d_ini2, f = d_ini2$govoff), England = list(d_ini2))

# model and actual values

ini_reg2 <- lapply(d_ini2, qm, exposure = 'fyearI', outcome = 'fageI', pyears = 1990:2019)
for (i in seq_along(ini_reg2)) {ini_reg2[[i]]$region <- names(ini_reg2)[i]}
ini_act2 <- lapply(d_ini2, function(x) x[, .(md = quantile(fageI, probs = 0.5), lower = quantile(fageI, probs = 0.25), upper = quantile(fageI, probs = 0.75)), fyearI])
for (i in seq_along(ini_act2)) {ini_act2[[i]]$region <- names(ini_act2)[i]}

ini_reg2 <- rbindlist(ini_reg2)
ini_act2 <- rbindlist(ini_act2)

#  :::::::::::::::::::::
#  duration of injecting
#  .....................

ndy <- 1990:2019

# -- percentage by year --    
    
id <- d[, .(n = .N), c('year', 'injdurI')]
id <- id[CJ(year = ndy, injdurI = 0:100), on = c('year', 'injdurI')]
id$n[is.na(id$n)] <- 0L
id <- id[, .(N = sum(n)), year][id, on = 'year']
id[, pc := n / N]
id <- id[injdurI <= 40 & year >= ndy[1]]
fwrite(id[, c('year', 'injdurI', 'pc')], 'uam_histogram_dur_1april2022.csv')

# -- percentage by initiation year --

iy <- d[, .(n = .N), c('year', 'fyearI')]
iy <- iy[CJ(year = ndy, fyearI = 1948:2019), on = c('year', 'fyearI')]
iy$n[is.na(iy$n)] <- 0L
iy <- iy[, .(N = sum(n)), year][iy, on = 'year']
iy[, pc := n / N]
iy <- iy[fyearI >= 1970 & fyearI <= year & year >= ndy[1]]
fwrite(iy[, c('year', 'fyearI', 'pc')], 'uam_histogram_yearInit_1april2022.csv')

# -- duration by region --

d_dur <- d[, c('year', 'govoff', 'injdurI')]
d_dur <- c(split(d_dur, f = d_dur$govoff), England = list(d_dur))

dur_reg <- lapply(d_dur, qm, outcome = 'injdurI')
for (i in seq_along(dur_reg)) {dur_reg[[i]]$region <- names(dur_reg)[i]}
dur_act <- lapply(d_dur, function(x) x[, .(md = quantile(injdurI, probs = 0.5), lower = quantile(injdurI, probs = 0.25), upper = quantile(injdurI, probs = 0.75)), year])
for (i in seq_along(dur_act)) {dur_act[[i]]$region <- names(dur_act)[i]}

dur_reg <- rbindlist(dur_reg)
dur_act <- rbindlist(dur_act)

#  ::::::::::::::::::
#  combine for saving
#  ..................

# for age of initiation, year is year of initiation (not survey year)

age_reg[, parameter := 'age']
ini_reg[, parameter := 'age_init']
ini_reg2[, parameter := 'age_init_whole_sample']
setnames(ini_reg, 'fyearI', 'year')
setnames(ini_reg2, 'fyearI', 'year')
dur_reg[, parameter := 'dur']
reg_results <- rbind(age_reg, ini_reg, ini_reg2, dur_reg)

age_act[, parameter := 'age']
ini_act[, parameter := 'age_init']
ini_act2[, parameter := 'age_init_whole_sample']
setnames(ini_act, 'fyearI', 'year')
setnames(ini_act2, 'fyearI', 'year')
dur_act[, parameter := 'dur']
act_results <- rbind(age_act, ini_act, ini_act2, dur_act)

fwrite(reg_results, 'uam_modelled_quantiles_1april2022.csv')
fwrite(act_results, 'uam_actual_quantiles_1april2022.csv')

#  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#  sensitivity excluding those reporting previous participation
#  ............................................................

prevPartSens <- function (d, outcome = 'ageI', exposure = 'year', label = 'age') {
  a <- qm(d = d, outcome = outcome, exposure = exposure, q = 0.5, sample = T)
  a$sens <- 'all'
  b <- qm(d = d[spit != 'Yes'], outcome = outcome, exposure = exposure, q = 0.5, sample = T)
  b$sens <- 'new'
  d <- qm(d = d[spit == 'Yes'], outcome = outcome, exposure = exposure, q = 0.5, sample = T)
  d$sens <- 'prev'
  p <- rbind(a, b, d)
  p$label <- label
  p
}

ageSpitSens <- prevPartSens(d = d, outcome = 'ageI', exposure = 'year', label = 'age')
durSpitSens <- prevPartSens(d = d, outcome = 'injdurI', exposure = 'year', label = 'dur')
iniSpitSens <- prevPartSens(d = d[(year - fyearI) <= 3], outcome = 'fageI', exposure = 'fyearI', label = 'ini')
names(iniSpitSens)[1] <- 'year'
save(ageSpitSens, durSpitSens, iniSpitSens, file = 'previousParticipationSensitivities_1april2022.Rdata')

#  ::::::::::::::::::::::::::::::::::::::::::::
#  pairwise ratios for modeling of cohort sizes
#  ............................................

pairwise_bs_variance <- function(decay = 15, B = 10, region = unique(d$govoff), minSurveyYear = 1990, name = NA, printProgress = T) {
  
  if (printProgress) {print(name)}
  
  bsd <- d[!is.na(fyearI) & !is.na(year) & (govoff %in% region) & (year >= minSurveyYear), c('year', 'fyearI')] # bootstrap data
  
  ratio_points <- rbindlist(lapply(1993:2019, function (z) {
    x <- bsd[year == z & fyearI >= 1980, .N, fyearI]
    y <- CJ(y1 = x$fyearI, y2 = x$fyearI)
    y <- y[y2 > y1]
    y <- data.table(y1 = x$fyearI, n1 = x$N)[y, on = 'y1']
    y <- data.table(y2 = x$fyearI, n2 = x$N)[y, on = 'y2']
    y[, half := exp(-(y2-y1)/decay)] # exponential decay function
    # y[, half := (1-decay)^(y2-y1)] # previous function where 1/15 stop per year
    y[, n1_rescale := n1/half]
    y[, ratio := n2 / n1_rescale]
    y[, syear := z]
    y[, c('syear', 'y1', 'y2', 'ratio')]
  }))
  
  bs_var <- function(z, B) { # z = survey year
    x <- bsd[year == z & fyearI >= 1980]
    n <- nrow(x)
    x <- x[sample(n, n * B, replace = T)]
    x[, b := rep(seq_len(B), each = n)]
    x <- x[, .N, c('fyearI', 'b')]
    yrs <- unique(x$fyearI)
    y <- CJ(y1 = yrs, y2 = yrs)
    y <- y[y2 > y1]
    n2 <- nrow(y)
    y <- rbindlist(rep(list(y), B))
    y[, b := rep(seq_len(B), each = n2)]
    y <- data.table(y1 = x$fyearI, n1 = x$N, b = x$b)[y, on = c('y1', 'b')]
    y <- data.table(y2 = x$fyearI, n2 = x$N, b = x$b)[y, on = c('y2', 'b')]
    y <- y[!is.na(n1) & !is.na(n2)]
    y[, half := exp(-(y2-y1)/decay)] # exponential decay function
    # y[, half := (1-decay)^(y2-y1)] # previous function where 1/15 stop per year
    y[, n1_rescale := n1/half]
    y[, ratio := n2 / n1_rescale]
    y[, logratio := log(ratio)]
    y[, .(syear = z, v = var(logratio)), c('y1', 'y2')]
  }
  
  bs_res <- lapply(1993:2019, bs_var, B = B)
  bs_res <- rbindlist(bs_res)
  
  bs_res[ratio_points, on = c('syear', 'y1', 'y2')]
  
}

# input values
# first regions, then England, then England 2008, then decay values

decay_values <- 10:20
regions <- d[!is.na(govoff), unique(govoff)]
regions_england <- regions[!(regions %in% c('Wales', 'Northern Ireland'))]

regions_input <- c(as.list(regions_england), rep(list(regions_england), length(decay_values) + 2))
decay_input <- c(rep(15, length(regions_england)+2), decay_values)
minYear_input <- c(rep(1990, length(regions_england) + 1), 2008, rep(1990, length(decay_values))) 

names_input <- c(regions_england, 'England', 'England2008', decay_values)

# bootstrap pairwise ratios and variances

set.seed(44)
bs <- mapply(pairwise_bs_variance,
             decay = decay_input,
             region = regions_input,
             minSurveyYear = minYear_input,
             name = names_input,
             B = 1000L,
             SIMPLIFY = F)
for(i in seq_along(bs)) {bs[[i]]$decay <- decay_input[i]}
for(i in seq_along(bs)) {bs[[i]]$minYear <- minYear_input[i]}
for(i in seq_along(bs)) {bs[[i]]$region <- names_input[i]}
bs <- rbindlist(bs)

#  do inverse variance weighting
#  .............................

bs[, iv := 1/v]
bs <- bs[, .(total_iv = sum(iv)), c('y1', 'y2', 'region')][bs, on = c('y1', 'y2', 'region')]
bs[, iwv := iv / total_iv]
bs[, logratio := log(ratio)]

# the variance of the combined effect is defined as the reciprocal of the sum of the weights

weighted_ratio <- bs[, .(logratio = sum(logratio * iwv), var_log_ratio = 1/sum(iv)), c('y1', 'y2', 'region', 'decay')]
weighted_ratio[, ratio := exp(logratio)]
weighted_ratio[, ratio_lci := exp(logratio - qnorm(0.975) * sqrt(var_log_ratio))]
weighted_ratio[, ratio_uci := exp(logratio + qnorm(0.975) * sqrt(var_log_ratio))]

fwrite(weighted_ratio, 'uam_pairwise_ratios_1april2022.csv')


#  :::::::::::::::::::::::::::::::::::::::::::::::
#  sensitivity analysis of changing the start year
#  ...............................................

start_year_input <- 1990:2019

set.seed(619)
bs <- mapply(pairwise_bs_variance,
             decay = 15,
             minSurveyYear = start_year_input,
             region = list(regions_england),
             name = start_year_input,
             B = 1000L,
             SIMPLIFY = F)
for(i in seq_along(bs)) {bs[[i]]$startYear <- start_year_input[i]}
bs <- rbindlist(bs)

#  do inverse variance weighting
#  .............................

bs[, iv := 1/v]
bs <- bs[, .(total_iv = sum(iv)), c('y1', 'y2', 'startYear')][bs, on = c('y1', 'y2', 'startYear')]
bs[, iwv := iv / total_iv]
bs[, logratio := log(ratio)]

# the variance of the combined effect is defined as the reciprocal of the sum of the weights

weighted_ratio <- bs[, .(logratio = sum(logratio * iwv), var_log_ratio = 1/sum(iv)), c('y1', 'y2', 'startYear')]
weighted_ratio[, ratio := exp(logratio)]
weighted_ratio[, ratio_lci := exp(logratio - qnorm(0.975) * sqrt(var_log_ratio))]
weighted_ratio[, ratio_uci := exp(logratio + qnorm(0.975) * sqrt(var_log_ratio))]

fwrite(weighted_ratio, 'uam_startYear_sensitivities_1april2022.csv')
