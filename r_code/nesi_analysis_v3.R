# This script processes NESI data. It requires raw NESI data.
# 1. Drops ineligible participants
# 2. Produces descriptive table
# 3. Reports quantiles of age, duration, and age at initiation
# 4. Estimates modelled quantiles
# 5. Estimates pairwise ratios for modelling of number injecting for the first time each year

setwd("//freddy/dept/PHIBCS/PHI/SHPN-PA/01 SHBBV/Surveillance/NESI/03 Data/Dan")

# -- libraries --

library(data.table) # for reading and manipulating data
library(quantreg) # for quantile regression
library(lubridate) # for date processing

# -- named functions --

# do quantile regression for interquartile range and predict confidence intervals
qm <- function(d, outcome = 'age', exposure = 'year', q = c(0.25, 0.5, 0.75), pyears = 1990:2019) {
  f <- as.formula(paste0(outcome, '~poly(', exposure, ',2)'))
  rbindlist(lapply(q, function(x) {
    m <- rq(f, data = d, tau = x)
    nd <- data.table(exposure = pyears)
    names(nd)[1] <- exposure
    cbind(nd, predict(m, newdata = nd, interval = 'confidence'), tau = x)
  }))
}

# -- read data --

d <- fread("data_exports/nesi_18sept2021_v3.csv", select = c('Date', 'Healthboard', 'Q2', 'Q3', 'Q15', 'Q18', 'Q19heroin', 'Q19cocaine', 'Q19crack', 'Q19HandCtogether', 'Q19Speed', 'Q19Benzo', 'Q19bodybuilding', 'Q19other', 'Q19_legal_high', 'time_since_last_injection', 'Age', 'length_inj_adj_np'))

d[, id := .I]
setnames(d, c('Q2', 'Q3', 'Q15', 'Q18'), c('sex', 'dob', 'age_init', 'inject_6_mth'))

# -- replace 9999, 8888, and 7777 with NA --

d[d == 9999] <- NA
d[d == 8888] <- NA
d[d == 7777] <- NA
replace_Q19 <- function (x) replace(x, is.na(x), 2L)
Q19s <- names(d)[grepl('Q19', names(d))]
d[, (Q19s) := lapply(.SD, replace_Q19), .SDcols = Q19s]

# -- health boards --

d[, area := as.character(Healthboard)]
d$area[d$area == 801] <- 'GCC'
d$area[d$area %in% c('802', '807')] <- 'Lothian_and_Tayside'
d$area[d$area %in% c(803:806, 808:811)] <- 'Other'

# -- process dates --

date_cols <- c('Date', 'dob')
d[, (date_cols) := lapply(.SD, mdy), .SDcols = date_cols]
d[, yr := year(Date)]
d[, syage := floor(Age)]
d[, date_init := Date - length_inj_adj_np * 365.25]
d[, year_init := year(date_init)]

# -- format other variables --

d[, sex := factor(sex, 1:2, c('male', 'female'))]

#  :::::::::::::::::
#  remove ineligible
#  .................

# -- remove 2020 --

d <- d[yr != 2020]

# -- remove those who did not inject drugs in the past 6 months --

nrow(d) # 17591
d[!(is.na(inject_6_mth) | inject_6_mth == 1), .N] # 3892
d <- d[is.na(inject_6_mth) | inject_6_mth == 1]
nrow(d) # 13699

# -- injecting bodybuilding drugs only --

illicit_vars <- c('Q19heroin', 'Q19cocaine', 'Q19crack', 'Q19HandCtogether', 'Q19Speed', 'Q19Benzo')
d[, illicit := rowSums(d[, illicit_vars, with = F] == 1, na.rm = T) > 0]
d[(Q19bodybuilding == 1 & illicit == F), .N] # 259
d <- d[!(Q19bodybuilding == 1 & illicit == F)]
nrow(d) # 13440

# -- missing variables -- 

d[, missing := is.na(Date) | is.na(Age) | is.na(length_inj_adj_np)]
d[missing == T, .N] # 85
d <- d[missing == F]
nrow(d) # 13355

#  ::::::::::::::
#  describe cases
#  ..............

d[, total := 1]
year_lims <- seq(1990, 2020, 5)
d[, yearGroup := findInterval(yr, year_lims)]
d[, yearGroup := factor(yearGroup, seq_along(year_lims), year_lims)]
d[, crackOrCoke := Q19crack == 1 | Q19cocaine == 1 | Q19HandCtogether == 1]
d[, heroin := Q19heroin == 1 | Q19HandCtogether == 1]
d[, other := Q19Speed == 1 | Q19Benzo == 1 | Q19bodybuilding == 1 | Q19other == 1 | Q19_legal_high == 1]

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

desctab <- rbindlist(lapply(c('total', 'sex', 'yearGroup', 'area', 'heroin', 'crackOrCoke', 'other'), dv))
fwrite(desctab, 'results/nesi_desc_2oct2021.csv')

#  ::::::::::::::::::::::::::::::::::::
#  split data for quantile calculations
#  ....................................

dq <- c(split(d, f = d$area), Scotland = list(d))

#  :::::::::::::::::::
#  age of participants
#  ...................

age_reg <- lapply(dq, qm, outcome = 'Age', exposure = 'yr', pyears = c(2008:2020))
for (i in seq_along(age_reg)) {age_reg[[i]]$region <- names(age_reg)[i]}
age_act <- lapply(dq, function(x) x[, .(md = quantile(Age, probs = 0.5), lower = quantile(Age, probs = 0.25), upper = quantile(Age, probs = 0.75)), yr])
for (i in seq_along(age_act)) {age_act[[i]]$region <- names(age_act)[i]}

age_reg <- rbindlist(age_reg)
age_act <- rbindlist(age_act)

#  :::::::::::::::::
#  age of initiation
#  .................

ndy <- 1980:2019

d_ini <- d[year_init %in% ndy]
d_ini <- d_ini[(yr - year_init) <= 3]
d_ini <- c(split(d_ini, f = d_ini$area), Scotland = list(d_ini))

# model and actual values

ini_reg <- lapply(d_ini, qm, exposure = 'year_init', outcome = 'age_init', pyears = 2005:2019)
for (i in seq_along(ini_reg)) {ini_reg[[i]]$region <- names(ini_reg)[i]}
ini_act <- lapply(d_ini, function(x) x[, .(md = quantile(age_init, probs = 0.5), lower = quantile(age_init, probs = 0.25), upper = quantile(age_init, probs = 0.75)), year_init])
for (i in seq_along(ini_act)) {ini_act[[i]]$region <- names(ini_act)[i]}

ini_reg <- rbindlist(ini_reg)
ini_act <- rbindlist(ini_act)

#  :::::::::::::::::::::
#  duration of injecting
#  .....................

dur_reg <- lapply(dq, qm, outcome = 'length_inj_adj_np', exposure = 'yr', pyears = 2008:2019)
for (i in seq_along(dur_reg)) {dur_reg[[i]]$region <- names(dur_reg)[i]}
dur_act <- lapply(dq, function(x) x[, .(md = quantile(length_inj_adj_np, probs = 0.5), lower = quantile(length_inj_adj_np, probs = 0.25), upper = quantile(length_inj_adj_np, probs = 0.75)), yr])
for (i in seq_along(dur_act)) {dur_act[[i]]$region <- names(dur_act)[i]}

dur_reg <- rbindlist(dur_reg)
dur_act <- rbindlist(dur_act)

#  ::::::::::::::::::
#  combine for saving
#  ..................

# for age of initiation, year is year of initiation (not survey year)

age_reg[, parameter := 'age']
ini_reg[, parameter := 'age_init']
setnames(ini_reg, 'year_init', 'yr')
dur_reg[, parameter := 'dur']
reg_results <- rbind(age_reg, ini_reg, dur_reg)

age_act[, parameter := 'age']
ini_act[, parameter := 'age_init']
setnames(ini_act, 'year_init', 'yr')
dur_act[, parameter := 'dur']
act_results <- rbind(age_act, ini_act, dur_act)

fwrite(reg_results, 'results/nesi_modelled_quantiles_2oct2021.csv')
fwrite(act_results, 'results/nesi_actual_quantiles_2oct2021.csv')

#  :::::::::::::::::::::
#  histogram of duration
#  .....................

ndy <- 2008:2019

# -- percentage by year --    

d[, dur_int := round(length_inj_adj_np, 0)]
id <- d[, .(n = .N), c('yr', 'dur_int')]
id <- id[CJ(yr = ndy, dur_int = 0:100), on = c('yr', 'dur_int')]
id$n[is.na(id$n)] <- 0L
id <- id[, .(N = sum(n)), yr][id, on = 'yr']
id[, pc := n / N]
id <- id[dur_int <= 40 & yr >= ndy[1]]
fwrite(id[, c('yr', 'dur_int', 'pc')], 'results/nesi_histogram_dur_2oct2021.csv')

# -- percentage by initiation year --

iy <- d[, .(n = .N), c('yr', 'year_init')]
iy <- iy[CJ(yr = ndy, year_init = 1948:2019), on = c('yr', 'year_init')]
iy$n[is.na(iy$n)] <- 0L
iy <- iy[, .(N = sum(n)), yr][iy, on = 'yr']
iy[, pc := n / N]
iy <- iy[year_init >= 1970 & year_init <= yr & yr >= ndy[1]]
fwrite(iy[, c('yr', 'year_init', 'pc')], 'results/nesi_histogram_yearInit_2oct2021.csv')

#  ::::::::::::::::::::::::::::::::::::::::::::
#  pairwise ratios for modeling of cohort sizes
#  ............................................

pairwise_bs_variance <- function(decay = 1/15, B = 10, region = unique(d$govoff), name = NA, printProgress = T) {
  
  if (printProgress) {print(name)}
  
  bsd <- d[(area %in% region), c('yr', 'year_init')] # bootstrap data
  
  ratio_points <- rbindlist(lapply(1993:2019, function (z) {
    x <- bsd[yr == z & year_init >= 1980, .N, year_init]
    y <- CJ(y1 = x$year_init, y2 = x$year_init)
    y <- y[y2 > y1]
    y <- data.table(y1 = x$year_init, n1 = x$N)[y, on = 'y1']
    y <- data.table(y2 = x$year_init, n2 = x$N)[y, on = 'y2']
    y[, half := (1-decay)^(y2-y1)]
    y[, n1_rescale := n1/half]
    y[, ratio := n2 / n1_rescale]
    y[, syear := z]
    y[, c('syear', 'y1', 'y2', 'ratio')]
  }))
  
  bs_var <- function(z, B) { # z = survey year
    x <- bsd[yr == z & year_init >= 1980]
    n <- nrow(x)
    x <- x[sample(n, n * B, replace = T)]
    x[, b := rep(seq_len(B), each = n)]
    x <- x[, .N, c('year_init', 'b')]
    yrs <- unique(x$year_init)
    y <- CJ(y1 = yrs, y2 = yrs)
    y <- y[y2 > y1]
    n2 <- nrow(y)
    y <- rbindlist(rep(list(y), B))
    y[, b := rep(seq_len(B), each = n2)]
    y <- data.table(y1 = x$year_init, n1 = x$N, b = x$b)[y, on = c('y1', 'b')]
    y <- data.table(y2 = x$year_init, n2 = x$N, b = x$b)[y, on = c('y2', 'b')]
    y <- y[!is.na(n1) & !is.na(n2)]
    y[, half := (1-decay)^(y2-y1)]
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
# first regions, then Scotland, then decay values

decay_values <- 1/(10:20)
regions <- d[!is.na(area), unique(area)]
regions_input <- c(as.list(regions), rep(list(regions), length(decay_values) + 1))
decay_input <- c(rep(1/15, length(regions)+1), decay_values)
names_input <- c(regions, 'Scotland', paste0('1/', 10:20))

# bootstrap pairwise ratios and variances

set.seed(12)
bs <- mapply(pairwise_bs_variance,
             decay = decay_input,
             region = regions_input,
             name = names_input,
             B = 1000L,
             SIMPLIFY = F)
for(i in seq_along(bs)) {bs[[i]]$decay <- 1/decay_input[i]}
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

fwrite(weighted_ratio, 'results/nesi_pairwise_ratios_2oct2021.csv')