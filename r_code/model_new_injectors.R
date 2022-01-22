options(scipen = 999)

library(data.table)
library(RColorBrewer)
library(devEMF)

add.alpha <- function(cols, alpha) rgb(t(col2rgb(cols)/255), alpha = alpha)

#  read pairwise ratios of numbers of new initiators
#  .................................................

ratios_uam <- read.csv('https://raw.githubusercontent.com/danlewer/uam/main/ratios_of_new_initiators/uam_pairwise_ratios_6oct2021.csv', stringsAsFactors = F)
ratios_uam_decay <- read.csv('https://raw.githubusercontent.com/danlewer/uam/main/ratios_of_new_initiators/uam_decay_sensitivities_20jan2022.csv', stringsAsFactors = F)
ratios_uam_startYear <- read.csv('https://raw.githubusercontent.com/danlewer/uam/main/ratios_of_new_initiators/uam_startYear_sensitivities_21jan2022.csv', stringsAsFactors = F)
ratios_nesi <- read.csv('https://raw.githubusercontent.com/danlewer/uam/main/ratios_of_new_initiators/nesi_pairwise_ratios_2oct2021.csv', stringsAsFactors = F)
setDT(ratios_uam); setDT(ratios_nesi); setDT(ratios_uam_decay); setDT(ratios_uam_startYear)

ratios_uam[, survey := 'uam']
ratios_nesi[, survey := 'nesi']

ratios <- rbind(ratios_uam, ratios_nesi)
ratios_england <- ratios[region == 'England' & decay == 15]

# /// compare original ratios with decay sensitivities
compare_ratios <- ratios_uam[region == 'England' & decay == 15, .(y1 = y1, y2 = y2, r1 = ratio)][ratios_uam_decay[, c('y1', 'y2', 'decay', 'ratio')], on = c('y1', 'y2')]
compare_ratios[, dif := abs(r1 - ratio)]
compare_ratios[, mean(dif), decay] # v. similar for 'constant'

#  heatmap of pairwise ratios
#  ..........................

cols_red <- colorRampPalette(brewer.pal(9, 'Reds'))(100)
cols_blue <- rev(colorRampPalette(brewer.pal(9, 'Blues'))(100))
pos_range <- c(0, max(ratios_england$logratio))
neg_range <- c(0, min(ratios_england$logratio))
pos_interval <- seq(pos_range[1], pos_range[2], length.out = 100)
neg_interval <- seq(neg_range[2], neg_range[1], length.out = 100)
interval <- c(neg_interval, pos_interval)
ratios_england[, interval := findInterval(logratio, interval)]
ratios_england[, col := c(cols_blue, cols_red)[interval]]

RRs <- c(0.2, 0.3, 0.5, 0.75, 1, 1.25, 1.5)
ticks <- sapply(log(RRs), function(x) which.min(abs(x - interval)))
ys_height <- 2020-1981
total_rr <- -neg_range[2] + pos_range[2]
neg_height <- ys_height * (-neg_range[2] / total_rr)

ys_neg <- seq(1981, 1981 + neg_height, length.out = 101)
ys_pos <- seq(1981 + neg_height, 2020, length.out = 101)

emf('new_initiator_ratio_heatmap_20jan2022.emf', height = 7, width = 9, family = 'Candara')

par(xpd = NA, mar = c(4, 4, 2, 8))
plot(1, type = 'n', xlim = c(1980, 2019), ylim = c(1980, 2020), axes = F, xlab = NA, ylab = NA)
with(ratios_england, rect(y1, y2, y1+1, y2+1, col = col))
axis(1, 1980:2019, label = F, pos = 1981, tck = -0.01)
axis(1, 1980 + 0:7 * 5, label = F, pos = 1981, tck = -0.03)
text(1980 + 0:7 * 5 + 0.5, y = 1980, label = 1980 + 0:7 * 5, adj = 1, srt = 90)
axis(2, 1981:2020, label = F, pos = 1980, tck = -0.01)
axis(2, 1985 + 0:7 * 5, label = F, pos = 1980, tck = -0.03)
text(1979, 1985 + 0:6 * 5 + 0.5, label = 1985 + 0:6 * 5, adj = 1)
rect(1980, 1981, 2019, 2020)
title(xlab = 'Year 1', line = 2)
title(ylab = 'Year 2', line = 2)
rect(2020, ys_neg[-length(ys_neg)], 2021, ys_neg[-1], col = cols_blue, border = NA)
rect(2020, ys_pos[-length(ys_pos)], 2021, ys_pos[-1], col = cols_red, border = NA)
rect(2020, 1981, 2021, 2020)
segments(2021, c(ys_neg, ys_pos)[ticks], 2021.5)
text(2022, c(ys_neg, ys_pos)[ticks], RRs, adj = 0)
text(2020, 2022, 'Ratio of new intiators\nYear 2 / Year 1', adj = 0)

dev.off()

#  do monte-carlo model to create estimates of ratios vs. 1980
#  ...........................................................

mc <- function(weighted_data, mcB = 1000, printProgress = T) { # monte-carlo function. mcB = number of simulations
  if (printProgress) {print(weighted_data$group[1])}
  sapply(seq_len(mcB), function(x) {
    ty <- sample(1980:2019, 1)
    y <- weighted_data[(y1 == ty | y2 == ty)]
    # standard deviation is square root of variance
    y[, sampled_log_ratio := rnorm(.N, mean = logratio, sd = sqrt(var_log_ratio))]
    y[, sampled_ratio := exp(sampled_log_ratio)]
    y[, sampled_ratio2 := fifelse(y2 == ty, 1/sampled_ratio, sampled_ratio)]
    y[, yr := y1 + y2 - ty]
    sr <- y$sampled_ratio2[match(1980:2019, y$yr)]
    sr[1980:2019 == ty] <- 1
    sr <- sr/sr[1]
    return(sr)
  })
}

ratios[, group := paste0(survey, '_', region)]
set.seed(14)
mc_results <- lapply(split(ratios, ratios$group), mc)
mc_estimates <- lapply(mc_results, function (x) apply(x, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = T))

# decay sensitivity estimates

set.seed(77)
mc_results_decay <- lapply(split(ratios_uam_decay, ratios_uam_decay$decay), mc)
mc_estimates_decay <- lapply(mc_results_decay, function (x) apply(x, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = T))

# start year sensitivity estimates

set.seed(811)
mc_results_startYear <- lapply(split(ratios_uam_startYear, ratios_uam_startYear$startYear), mc)
mc_estimates_startYear <- lapply(mc_results_startYear, function (x) apply(x, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = T))

#  plot of main scenario
#  .....................

cols <- brewer.pal(3, 'Set1')[1:2]
cols2 <- add.alpha(cols, 0.2)
rrs <- c(0.1, 0.2, 0.4, 0.7, 1, 2, 3, 5)

emf('main_scenario_ratio_1980_20jan2022.emf', height = 5, width = 8, family = 'Candara')

par(mar = c(4, 4, 0, 11), xpd = NA)
plot(1, type = 'n', xlim = c(1980, 2019), ylim = range(log(rrs)), axes = F, xlab = NA, ylab = NA)
segments(1980, 0, 2019)

pd <- log(mc_estimates$uam_England)
polygon(c(1980:2019, 2019:1980), c(pd[1,], rev(pd[3,])), col = cols2[2], border = NA)
lines(1980:2019, pd[2,], col = cols[2])

pd <- log(mc_estimates$nesi_Scotland)
polygon(c(1980:2019, 2019:1980), c(pd[1,], rev(pd[3,])), col = cols2[1], border = NA)
lines(1980:2019, pd[2,], col = cols[1])

axis(1, 1980:2019, pos = log(rrs[1]), labels = F, tck = -0.01)
axis(1, 1980 + 5 * 0:7, pos = log(rrs[1]), tck = -0.03)
axis(2, log(rrs), rrs, pos = 1980, las = 2)
rect(1980, log(rrs[1]), 2019, log(max(rrs)))
title(xlab = 'Year', line = 2.5)
title(ylab = 'Ratio of people injecting\nfor the first time, vs. 1980', line = 2)

ys <- seq(-0.6, 0.6, length.out = 6)
rect(2020, ys[1], 2022, ys[3], col = cols2[1], border = NA)
rect(2020, ys[4], 2022, ys[6], col = cols2[2], border = NA)
segments(2020, ys[c(2, 5)], x1 = 2022, col = cols)
text(2022.5, ys[c(2, 5)], c('NESI (Scotland)', 'UAM (England, Wales &\nNorthern Ireland'), adj = 0)

dev.off()

#  plot of sensitivity analyses on lambda
#  ......................................

cols <- colorRampPalette(c('blue', 'red'))(length(10:20))
rrs <- c(0.1, 0.2, 0.4, 0.7, 1, 2, 3)
uam_sens <- mc_estimates[paste0('uam_1/', 10:20)]
nesi_sens <- mc_estimates[paste0('nesi_1/', 10:20)]
uam_final_values <- sapply(uam_sens, function(x) x[2,40])
nesi_final_values <- sapply(nesi_sens, function(x) x[2,40])

emf('lambda_sensitivty_analysis_20jan2022.emf', height = 10, width = 7, family = 'Candara')

par(mfrow = c(2, 1), mar = c(4, 4, 0, 6), xpd = NA)

# uam

plot(1, type = 'n', xlim = c(1980, 2019), ylim = range(log(rrs)), axes = F, xlab = NA, ylab = NA)
segments(1980, 0, 2019)
for (i in seq_along(10:20)) {
  pd <- log(uam_sens[[i]])
  lines(1980:2019, pd[2,], col = cols[i])
}
axis(1, 1980:2019, pos = log(rrs[1]), labels = F, tck = -0.01)
axis(1, 1980 + 5 * 0:7, pos = log(rrs[1]), tck = -0.03)
axis(2, log(rrs), rrs, pos = 1980, las = 2)
rect(1980, log(rrs[1]), 2019, log(max(rrs)))
title(ylab = 'Ratio of people injecting\nfor the first time, vs. 1980', line = 2)
text(2020, log(uam_final_values), 10:20, col = cols, adj = 0, cex = 0.7)
text(2020, max(log(uam_final_values)) + 0.2, 'Average duration\nof injecting', adj = 0, cex = 0.7)
text(1980.5, log(2.5), 'UAM (England, Wales and\nNorthern Ireland)', adj = 0)

# nesi

plot(1, type = 'n', xlim = c(1980, 2019), ylim = range(log(rrs)), axes = F, xlab = NA, ylab = NA)
segments(1980, 0, 2019)
for (i in seq_along(10:20)) {
  pd <- log(nesi_sens[[i]])
  lines(1980:2019, pd[2,], col = cols[i])
}
axis(1, 1980:2019, pos = log(rrs[1]), labels = F, tck = -0.01)
axis(1, 1980 + 5 * 0:7, pos = log(rrs[1]), tck = -0.03)
axis(2, log(rrs), rrs, pos = 1980, las = 2)
rect(1980, log(rrs[1]), 2019, log(max(rrs)))
title(xlab = 'Year', line = 2.5)
title(ylab = 'Ratio of people injecting\nfor the first time, vs. 1980', line = 2)
text(2020, log(nesi_final_values), 10:20, col = cols, adj = 0, cex = 0.7)
text(2020, max(log(nesi_final_values)) + 0.2, 'Average duration\nof injecting', adj = 0, cex = 0.7)
text(1980.5, log(2.5), 'NESI (Scotland)', adj = 0)

dev.off()

# sensitivity on lambda changing over time

emf('lambda_sensitivty_analysis_2_20jan2022.emf', height = 5, width = 7, family = 'Candara')

par(mar = c(4, 4, 1, 10), xpd = NA)

cols <- colorRampPalette(c('blue', 'red'))(length(1:3))
plot(1, type = 'n', xlim = c(1980, 2019), ylim = range(log(rrs)), axes = F, xlab = NA, ylab = NA)
segments(1980, 0, 2019)
for (i in seq_along(1:3)) {
  pd <- log(mc_estimates_decay[[i]])
  lines(1980:2019, pd[2,], col = cols[i])
}
axis(1, 1980:2019, pos = log(rrs[1]), labels = F, tck = -0.01)
axis(1, 1980 + 5 * 0:7, pos = log(rrs[1]), tck = -0.03)
axis(2, log(rrs), rrs, pos = 1980, las = 2)
rect(1980, log(rrs[1]), 2019, log(max(rrs)))
title(xlab = 'Year', line = 2.5)
title(ylab = 'Ratio of people injecting\nfor the first time, vs. 1980', line = 2)
ys <- seq(0, 1, length.out = 3)
segments(2020, ys, 2021.5, ys, col = cols)
text(2022, ys, c('Constant at\n1/15 per year', 'Decreasing from 1/10 to\n1/20 per year', 'Increasing from 1/20 to\n1/10 per year'), adj = 0)

dev.off()

#  estimate absolute size of cohorts
#  .................................

# -- existing estimate of number of PWID in 2011, by region. From https://www.drugsandalcohol.ie/21931/1/estimates-of-the-prevalence-of-opiate-use-and-or-crack-cocaine-use-2011-12.pdf #
# -- scotland population in 2006: http://eprints.gla.ac.uk/45433/3/45433.pdf 

targets <- data.table(region = c('East of England', 'East Midlands', 'London', 'North East', 'North West', 'South East', 'South West', 'West Midlands', 'Yorkshire & Humber', 'England', 'GCC', 'Lothian_and_Tayside', 'Other', 'Scotland'),
                      country = rep(c('England', 'Scotland'), c(10, 4)),
                      target = c(6650, 7808, 11351, 6334, 13110, 11047, 10134, 9175, 11692, 87302, 8862, 4516, 10354, 23933),
                      lowerCI = c(5995, 6989, 10711, 5948, 12233, 9635, 9474, 8281, 11024, 85307, 7091, 3536, 7857, 21655),
                      upperCI = c(7386, 8630, 12347, 6770, 14305, 12368, 10958, 10082, 12457, 90353, 11330, 5962, 10517, 27143),
                      year = c(rep(2011, 10), rep(2006, 4)),
                      stringsAsFactors = F)

# -- assess symmetry of confidence intervals --
# log looks slightly better

untransformed_ci_width <- with(targets, cbind(target - lowerCI, upperCI - target))
logged_ci_width <- with(log(targets[, -c('region', 'country')]), cbind(target - lowerCI, upperCI - target))

data.frame(region = targets$region,
           untransformed = apply(untransformed_ci_width, 1, diff) / untransformed_ci_width[,1] * 100,
           logged = apply(logged_ci_width, 1, diff) / logged_ci_width[,1] * 100)

targets$log.target <- log(targets$target)
targets$log.lowerCI <- log(targets$lowerCI)
targets$log.upperCI <- log(targets$upperCI)
targets$log.meanCIwidth <- with(targets, rowMeans(cbind(log.target - log.lowerCI, log.upperCI - log.target)))
targets$log.CIsd <- targets$log.meanCIwidth / qnorm(0.975)

# check results are similar to reported CIs
# works well

cbind(targets, t(mapply(function(...) quantile(exp(rnorm(1e5, ...)), probs = c(0.5, 0.025, 0.975)), 
                        mean = targets$log.target, 
                        sd = targets$log.CIsd)))

# -- function estimating population size each year, given the number starting in 1980, ratios of new starters vs. 1980, and rate of quitting

mm <- function(ratios, new1980 = 10000, decay = 1/15) {
  new15a <- new1980 * ratios
  m <- diag(new15a)
  m <- t(new15a * (1-decay)^(col(m)-row(m))) * (col(m) <= row(m))
  r <- rowSums(m, na.rm = T)
  r[colSums(is.na(m)) > 0] <- NA
  r
}

# -- estimate regional and national cohort sizes with prediction intervals --

find_multiple_cohorts <- function(ratio_matrix, log.target, log.target.sd, targetYear = 2011, test_vals = 0:20 * 1000, printTestVal = T) {
  B <- ncol(ratio_matrix)
  tg <- round(exp(rnorm(B, mean = log.target, sd = log.target.sd)), 0)
  y <- sapply(test_vals, function(x) {
    pop <- apply(ratio_matrix, 2, mm, new1980 = x)
    pop[1980:2011 == targetYear,]
  })
  closest_to_target <- apply(abs(y - tg), 1, which.min)
  closest_to_target[sapply(closest_to_target, length) == 0] <- NA
  closest_to_target <- unlist(closest_to_target)
  vals1980 <- test_vals[closest_to_target]
  new_cohorts <- rep(vals1980, each = length(1980:2019)) * ratio_matrix
  return(new_cohorts)
}

targets2 <- data.table(label = c(paste0('nesi_1/', 10:20),
                                 paste0('nesi_', c('GCC', 'Lothian_and_Tayside', 'Other', 'Scotland')),
                                 paste0('uam_1/', 10:20),
                                 paste0('uam_', c('East of England', 'East Midlands', 'London', 'North East', 'North West', 'South East', 'South West', 'West Midlands', 'Yorkshire & Humber', 'England', 'England2008'))),
                       region = c(rep('Scotland', 11),
                                  'GCC', 'Lothian_and_Tayside', 'Other', 'Scotland',
                                  rep('England', 11),
                                  'East of England', 'East Midlands', 'London', 'North East', 'North West', 'South East', 'South West', 'West Midlands', 'Yorkshire & Humber', 'England', 'England'))
targets2 <- targets[targets2, on = 'region']

# # -- This code takes a long time to run, so the results are saved as an 'Rdata' object and loaded
# set.seed(156)
# mc_cohort <- lapply(targets2$label, function (x) {
#   print(x)
#   find_multiple_cohorts(ratio_matrix = mc_results[[x]],
#                         log.target = targets2[label == x, log.target],
#                         log.target.sd = targets2[label == x, log.CIsd],
#                         targetYear = targets2[label == x, year],
#                         test_vals = 10 * seq_len(2000))
# })
# names(mc_cohort) <- targets2$label
# save(mc_cohort, file = 'mc_cohort_6Oct2021.Rdata')
load(url('https://github.com/danlewer/uam/raw/main/ratios_of_new_initiators/mc_cohort_6Oct2021.Rdata'))

cohort_size <- lapply(mc_cohort, function(x) apply(x, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = T))
cohort_size <- lapply(cohort_size, t)
cohort_size <- lapply(cohort_size, data.frame)
cohort_size <- lapply(cohort_size, `names<-`, value = c('lower', 'point', 'upper'))

#  -- point estimates --

test_vals <- 10 * seq_len(2000)
point_ratios <- lapply(mc_results[targets2$label], function(x) apply(x, 1, median, na.rm = T))
vals1980 <- sapply(targets2$label, function(z) test_vals[which.min(abs(sapply(test_vals, function(x) mm(point_ratios[[z]], new1980 = x)[1980:2019 == 2011]) - targets2[label == z, target]))])
point_estimates <- data.table(cbind(year = 1980:2019, round(mapply(`*`, a = point_ratios, b = vals1980), -1)))
for(i in targets2$label) {
  cohort_size[[i]]$point <- point_estimates[[i]]
}

#  plot of cohort sizes for main article
#  .....................................

cols <- brewer.pal(4, 'Paired')

region_labels_nesi <- c('GCC', 'Lothian_and_Tayside', 'Other')
region_labels_uam <- c('East of England', 'East Midlands', 'London', 'North East', 'North West', 'South East', 'South West', 'West Midlands', 'Yorkshire & Humber')
region_labels <- data.table(label = c(paste0('nesi_', region_labels_nesi), paste0('uam_', region_labels_uam)),
                            tag = c(region_labels_nesi, region_labels_uam),
                            country = c(rep('Scotland', length(region_labels_nesi)), rep('England', length(region_labels_uam))))
region_labels[, tag := gsub('_', ' ', tag)]
region_labels$tag[region_labels$tag == 'GCC'] <- 'Greater Glasgow & Clyde'
region_labels$tag[region_labels$tag == 'Other'] <- 'Other regions of Scotland'
region_labels[, col_line := as.character(factor(country, c('Scotland', 'England'), cols[c(2, 4)]))]
region_labels[, col_shad := as.character(factor(country, c('Scotland', 'England'), cols[c(1, 3)]))]

panels <- cohort_size[region_labels$label]
max_year <- sapply(panels, function(x) which.max(x$point))
max_cohort <- sapply(panels, function(x) max(x$point))
region_order <- names(max_year)[order(max_year, -max_cohort, decreasing = T)]
panels <- panels[region_order]
n_panels <- length(panels)
ymax <- 2500
yoff <- (seq_along(panels) - 1) * ymax
ygrids <- 1:4 * 500

#emf('cohort_sizes_2oct2021.emf', width = 10, height = 9, family = 'Corbel')
cairo_pdf('cohort_sizes_6oct2021.pdf', width = 10, height = 9, family = 'Corbel')

layout(matrix(1:2, ncol = 2), widths = c(4, 2))
par(mar = c(0, 0, 0, 0), oma = c(4, 5, 0, 4), xpd = NA)

ygap <- 0.1

plot(1, type = 'n', xlim = c(1980, 2019), ylim = c(0, 2 + ygap), axes = F, xlab = NA, ylab = NA)

ymax1 <- 15000
segments(c(1990, 2000, 2010), 1 + ygap, y1 = 2 + ygap)
segments(1980, c(5000, 10000) / ymax1 + 1 + ygap, x1 = 2019, lty = 3)
with(cohort_size$uam_England, {
  polygon(c(1980:2019, 2019:1980), c(lower, rev(upper)) / ymax1 + 1 + ygap, col = cols[3], border = NA)
  lines(1980:2019, point / ymax1 + 1 + ygap, col = cols[4], lwd = 2)
})
axis(2, 0:6 * 2500 / ymax1 + 1 + ygap, labels = prettyNum(0:6 * 2500, big.mark = ','), pos = 1980, las = 2)
rect(1980, 1 + ygap, 2019, 2 + ygap)
text(1981, 1.95 + ygap, 'England', adj = 0)

ymax2 <- 3500
segments(c(1990, 2000, 2010), 0, y1 = 1)
segments(1980, c(1000, 2000, 3000) / ymax2, x1 = 2019, lty = 3)
with(cohort_size$nesi_Scotland, {
  polygon(c(1980:2019, 2019:1980), c(lower, rev(upper)) / ymax2, col = cols[1], border = NA)
  lines(1980:2019, point / ymax2, col = cols[2], lwd = 2)
})
axis(2, 0:7 * 500 / ymax2, labels = prettyNum(0:7 * 500, big.mark = ','), pos = 1980, las = 2)
rect(1980, 0, 2019, 1)
text(1981, 0.95, 'Scotland', adj = 0)

axis(1, 1980:2019, pos = 0, labels = F, tck = -0.01)
axis(1, 1980 + 5 * 0:7, pos = 0, tck = -0.03)

plot(1, type = 'n', xlim = c(1980, 2019), ylim = c(0, ymax * n_panels), axes = F, xlab = NA, ylab = NA)
segments(c(1990, 2000, 2010), 0, y1 = ymax * n_panels)
for(i in seq_along(panels)) {
  segments(1980, ygrids + yoff[i], 2019, ygrids + yoff[i], lty = 3)
  pd <- panels[[i]]
  pd <- pd[complete.cases(pd),]
  years <- seq_len(nrow(pd)) + 1979
  with(pd, {
    polygon(c(years, rev(years)), c(lower, rev(upper)) + yoff[i], col = region_labels[label == names(panels)[i], col_shad], border = NA)
    lines(years, point + yoff[i], col = region_labels[label == names(panels)[i], col_line], lwd = 2)
  })
  rect(1980, yoff[i], 2019, yoff[i] + ymax)
  text(1981, yoff[i] + ymax * 0.9, region_labels[label == names(panels)[i], tag], adj = 0, cex = 0.7)
}
axis(1, 1980:2019, pos = 0, labels = F, tck = -0.01)
axis(1, 1980 + 5 * 0:7, pos = 0, tck = -0.03)
text(x = 2019.5, 
     y = rep(seq_len(n_panels) - 1, each = length(ygrids)) * ymax + rep(ygrids, times = n_panels),
     labels = prettyNum(rep(ygrids, times = n_panels), big.mark = ','),
     adj = 0, cex = 0.7)

mtext('Number of people injecting drugs for the first time', outer = T, side = 2, line = 3.5)

dev.off()

#  plot of sensitivites for main article
#  .....................................

uam_sens_abs <- lapply(cohort_size[paste0('uam_1/', 10:20)], function (x) x[,2])
cols <- colorRampPalette(c('blue', 'red'))(length(10:20))

cairo_pdf('Figure5.pdf', width = 9, height = 6, family = 'Candara')

par(mar = c(5, 5, 1, 10), xpd = NA)
plot(1, type = 'n', xlim = c(1980, 2019), ylim = c(0, 14000), axes = F, xlab = NA, ylab = NA)
rect(1980, 0, 2019, 14000)
axis(1, 1980:2019, pos = 0, labels = F, tck = -0.01)
axis(1, 1980 + 5 * 0:7, pos = 0, tck = -0.02)
axis(2, seq(0, 14000, 2000), labels = prettyNum(seq(0, 14000, 2000), big.mark = ','), pos = 1980, las = 2)
lwds <- rep(1, 11)
lwds[6] <- 3
mapply(lines, x = list(1980:2019), y = uam_sens_abs, col = cols, lwd = lwds)
ys <- seq(3500, 10500, length.out = 11)
segments(2020, ys, 2022, ys, col = cols, lwd = lwds)
labs <- paste0('1/', 10:20)
labs[6] <- '1/15 (main results)'
text(2023, ys, labs, adj = 0)
text(2020, max(ys) * 1.15, 'Proportion that stop\ninjecting per year', adj = 0)
title(xlab = 'Year', line = 2)
title(ylab = 'Number of people injecting drugs for the first time')

dev.off()

#  table of cohort sizes for supplementary information
#  ...................................................

roundf <- function (x) {
  y <- round(x, -1)
  paste0(y$point, ' (', y$lower, '-', y$upper, ')')
}

si_absolute <- sapply(cohort_size[c('uam_England', 'nesi_Scotland', region_labels$label)], roundf)
colnames(si_absolute) <- c('England', 'Scotland', region_labels$tag)

# check regions roughly total England/Scotland

compare_england <- cbind(rowSums(sapply(cohort_size[region_labels[country == 'England', label]], function(x) x[,2])),
                         cohort_size$uam_England[,2])
compare_scotland <- cbind(rowSums(sapply(cohort_size[region_labels[country == 'Scotland', label]], function(x) x[,2])),
                          cohort_size$nesi_Scotland[,2])
par(mfrow = c(1, 2))
plot(compare_england[,1], type = 'l', ylim = c(0, 12000), ylab = NA, xlab = NA, main = 'England'); lines(compare_england[,2], col = 'red')
plot(compare_scotland[,1], type = 'l', ylim = c(0, 4000), ylab = NA, xlab = NA, main = 'Scotland'); lines(compare_scotland[,2], col = 'red')

si_absolute <- cbind(year = 1980:2019, as.data.frame.matrix(si_absolute))

fwrite(si_absolute, 'absolute_numbers_6oct2021.csv')

#  comparison of models with and without data before 2008
#  ......................................................

cols <- brewer.pal(3, 'Set1')[1:2]
cols2 <- add.alpha(cols, alpha = 0.3)

emf('sensitivity_2008_6oct2021.emf', height = 5, width = 7, family = 'Corbel')

par(mar = c(4, 5, 0, 7), xpd = NA)
plot(1, type = 'n', xlim = c(1980, 2019), ylim = c(0, 15000), axes = F, xlab = NA, ylab = NA)
with(cohort_size$uam_England, {
  lines(1980:2019, point, col = cols[1], lwd = 3)
  polygon(c(1980:2019, 2019:1980), c(lower, rev(upper)), col = cols2[1], border = NA)
})
with(cohort_size$uam_England2008, {
  lines(1980:2019, point, col = cols[2], lwd = 3)
  polygon(c(1980:2019, 2019:1980), c(lower, rev(upper)), col = cols2[2], border = NA)
})
axis(1, 1980:2019, labels = F, pos = 0, tck = -0.01)
axis(1, seq(1980, 2015, 5), pos = 0, tck = -0.03)
axis(2, seq(0, 15000, 2500), pos = 1980, las = 2)
rect(1980, 0, 2019, 15000)
title(ylab = 'Number of people injecting\ndrugs for the first time')
title(xlab = 'Year of initiation', line = 2.5)
ys <- seq(5000, 10000, length.out = 6)
rect(2020, ys[1], 2022, ys[3], col = cols2[1], border = NA)
rect(2020, ys[4], 2022, ys[6], col = cols2[2], border = NA)
segments(2020, ys[c(2, 5)], x1 = 2022, col = cols, lwd = 2)
text(2022.5, ys[c(2, 5)], c('England\n(all data)', 'England\n(2008-2019)'), adj = 0)

dev.off()

#  :::::::::::::::::::::::::::::::::::::::::
#  plot of senstivities by survey start year
#  .........................................

# -- This code takes a long time to run, so the results are saved as an 'Rdata' object and loaded
# set.seed(332)
# mc_cohort_startYear <- mapply(find_multiple_cohorts,
#                               ratio_matrix = mc_results_startYear,
#                               log.target = targets2[label == 'uam_England', log.target],
#                               log.target.sd = targets2[label == 'uam_England', log.CIsd],
#                               targetYear = targets2[label == 'uam_England', year],
#                               test_vals = list(10 * seq_len(2000)),
#                               SIMPLIFY = F)
# save(mc_cohort_startYear, file = 'mc_cohort_startYear_21jan2022.Rdata')

load(url('https://github.com/danlewer/uam_nesi/blob/main/ratios_of_new_initiators/mc_cohort_startYear_21jan2022.Rdata?raw=true'))

cohort_size_startYear <- lapply(mc_cohort_startYear, function(x) apply(x, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = T))
cohort_size_startYear <- lapply(cohort_size_startYear, t)
cohort_size_startYear <- lapply(cohort_size_startYear, data.frame)
cohort_size_startYear <- lapply(cohort_size_startYear, `names<-`, value = c('lower', 'point', 'upper'))

#  -- point estimates --

test_vals <- 10 * seq_len(2000)
point_ratios <- lapply(mc_results_startYear, function(x) apply(x, 1, median, na.rm = T))
vals1980 <- sapply(as.character(1990:2019), function(z) test_vals[which.min(abs(sapply(test_vals, function(x) mm(point_ratios[[z]], new1980 = x)[1980:2019 == 2011]) - targets2[label == 'uam_England', target]))])
point_estimates <- data.table(cbind(year = 1980:2019, round(mapply(`*`, a = point_ratios, b = vals1980), -1)))
for(i in as.character(1990:2019)) {
  cohort_size_startYear[[i]]$point <- point_estimates[[i]]
}

# -- plot --- 

emf('sensitivity_survey_years_included_21jan2022.emf', height = 8, width = 15, family = 'Candara')

par(mar = c(5, 3, 1, 12), xpd = NA, mfrow = c(1, 2), oma = c(0, 3, 0, 0))

cols <- colorRampPalette(c(brewer.pal(3, 'Set1')[1], 'yellow', brewer.pal(3, 'Set1')[2]))(30)

plot(1, type = 'n', xlim = c(1980, 2019), ylim = c(0, 15000), axes = F, xlab = NA, ylab = NA)
rect(1980, 0, 2019, 15000)
axis(1, 1980:2019, pos = 0, labels = F, tck = -0.01)
axis(1, 1980 + 5 * 0:7, pos = 0, tck = -0.02)
axis(2, 0:5 * 3000, labels = prettyNum(0:5 * 3000, big.mark = ','), pos = 1980, las = 2)
lwds <- rep(0.2, 30)
lwds[c(1, 19)] <- 3
ltys <- rep(1, 30)
ltys[c(1, 19)] <- 1
mapply(lines, x = list(1980:2019), y = lapply(cohort_size_startYear, function (x) x[,2]), col = cols, lwd = lwds, lty = ltys)
ys <- seq(0, 14000, length.out = 30)
segments(2020, ys, 2022, ys, col = cols, lwd = lwds, lty = ltys)
labs <- paste0(1990:2019, ':2019')
labs[1] <- '1990:2019 (main results)'
labs[19] <- '2008:2019 (sensitivity)'
labs[length(labs)] <- '2019'
text(2023, ys, labs, adj = 0)
text(2020, max(ys) * 1.065, 'Survey years\nincluded', adj = 0)
title(xlab = 'Year', line = 2)
title(ylab = 'Number of people injecting drugs for the first time', line = 3.5)
text(1999.5, 15500, 'Point estimates')

cols <- brewer.pal(3, 'Set1')[1:2]
cols2 <- add.alpha(cols, 0.2)

plot(1, type = 'n', xlim = c(1980, 2019), ylim = c(0, 15000), axes = F, xlab = NA, ylab = NA)
rect(1980, 0, 2019, 15000)
axis(1, 1980:2019, pos = 0, labels = F, tck = -0.01)
axis(1, 1980 + 5 * 0:7, pos = 0, tck = -0.02)
axis(2, 0:5 * 3000, labels = prettyNum(0:5 * 3000, big.mark = ','), pos = 1980, las = 2)
polygon(x = c(1980:2019, 2019:1980), y = c(cohort_size_startYear$`1990`$lower, rev(cohort_size_startYear$`1990`$upper)), col = cols2[1], border = NA)
polygon(x = c(1980:2019, 2019:1980), y = c(cohort_size_startYear$`2019`$lower, rev(cohort_size_startYear$`2019`$upper)), col = cols2[2], border = NA)
lines(1980:2019, cohort_size_startYear$`1990`$point, col = cols[1])
lines(1980:2019, cohort_size_startYear$`2019`$point, col = cols[2])
ys <- seq(6000, 9000, length.out = 6)
rect(2020, ys[c(1, 4)], 2022, ys[c(3, 6)], col = cols2, border = NA)
segments(2020, ys[c(2, 5)], 2022, ys[c(2, 5)], col = cols)
text(2022.5, ys[c(2, 5)], c('1990:2019\n(main results)', '2019 only'), adj = 0)
text(2020, max(ys) * 1.1, 'Survey years\nincluded', adj = 0)
title(xlab = 'Year', line = 2)
title(ylab = 'Number of people injecting drugs for the first time', line = 3.5)
text(1999.5, 15500, 'Point estimates and 95% prediction intervals')

dev.off()
