options(scipen = 999)

library(data.table)
library(RColorBrewer)
library(devEMF)

#  read pairwise ratios of numbers of new initiators
#  .................................................

ratios_lambda <- read.csv('https://raw.githubusercontent.com/danlewer/uam/main/ratios_of_new_initiators/by_cessation_rate.csv')
ratios_region <- read.csv('https://raw.githubusercontent.com/danlewer/uam/main/ratios_of_new_initiators/by_region.csv', stringsAsFactors = F)
setDT(ratios_lambda); setDT(ratios_region)

#  heatmap of pairwise ratios
#  ..........................

decay15 <- ratios_lambda[round(abs(decay - 1/15), 6) < 0.00001]
cols_red <- colorRampPalette(brewer.pal(9, 'Reds'))(100)
cols_blue <- rev(colorRampPalette(brewer.pal(9, 'Blues'))(100))
pos_range <- c(0, max(decay15$logratio))
neg_range <- c(0, min(decay15$logratio))
pos_interval <- seq(pos_range[1], pos_range[2], length.out = 100)
neg_interval <- seq(neg_range[2], neg_range[1], length.out = 100)
interval <- c(neg_interval, pos_interval)
decay15[, interval := findInterval(logratio, interval)]
decay15[, col := c(cols_blue, cols_red)[interval]]

RRs <- c(0.2, 0.3, 0.5, 0.75, 1, 1.25, 1.5)
ticks <- sapply(log(RRs), function(x) which.min(abs(x - interval)))
ys_height <- 2020-1981
total_rr <- -neg_range[2] + pos_range[2]
neg_height <- ys_height * (-neg_range[2] / total_rr)

ys_neg <- seq(1981, 1981 + neg_height, length.out = 101)
ys_pos <- seq(1981 + neg_height, 2020, length.out = 101)

    emf('new_initiator_ratio_heatmap_3july2021.emf', height = 7, width = 9, family = 'Corbel')
    
    par(xpd = NA, mar = c(4, 4, 2, 8))
    plot(1, type = 'n', xlim = c(1980, 2019), ylim = c(1980, 2020), axes = F, xlab = NA, ylab = NA)
    with(decay15, rect(y1, y2, y1+1, y2+1, col = col))
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

mc <- function(weighted_data, mcB = 1000, filter = NA) { # monte-carlo function. mcB = number of simulations
  if (is.na(filter[1])) {filter <- T}
  sapply(seq_len(mcB), function(x) {
    ty <- sample(1980:2019, 1)
    y <- weighted_data[(y1 == ty | y2 == ty) & filter]
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

get_estimate_and_intervals <- function (x) apply(x, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = T)

decay_values <- 1/10:20
set.seed(14)
mc_decays <- lapply(decay_values, function(x) mc(ratios_lambda, filter = round(ratios_lambda$decay - x, 5) < 0.0001))
names(mc_decays) <- decay_values
mc_decays_estimates <- lapply(mc_decays, get_estimate_and_intervals)

regions <- unique(ratios_region$region)
set.seed(50)
mc_regions <- lapply(regions, function(x) mc(ratios_region, filter = ratios_region$region == x))
names(mc_regions) <- regions
mc_regions_estimates <- lapply(mc_regions, get_estimate_and_intervals)

#  plot of main scenario
#  .....................

cols <- brewer.pal(4, 'Paired')
rrs <- c(0.2, 0.4, 0.7, 1, 2, 3, 5)

    emf('main_scenario_ratio_1980_3july2021.emf', height = 5, width = 7, family = 'Corbel')
    
    par(mar = c(4, 4, 0, 0))
    plot(1, type = 'n', xlim = c(1980, 2019), ylim = range(log(rrs)), axes = F, xlab = NA, ylab = NA)
    segments(1980, 0, 2019)
    pd <- log(mc_decays_estimates[decay_values == 1/15][[1]])
    polygon(c(1980:2019, 2019:1980), c(pd[1,], rev(pd[3,])), col = cols[1], border = NA)
    lines(1980:2019, pd[2,], col = cols[2])
    axis(1, 1980:2019, pos = log(rrs[1]), labels = F, tck = -0.01)
    axis(1, 1980 + 5 * 0:7, pos = log(rrs[1]), tck = -0.03)
    axis(2, log(rrs), rrs, pos = 1980, las = 2)
    rect(1980, log(rrs[1]), 2019, log(max(rrs)))
    title(xlab = 'Year', line = 2.5)
    title(ylab = 'Ratio of new injectors vs. 1980', line = 2.5)
    
    dev.off()

#  plot of sensitivity analyses on lambda
#  ......................................

cols <- colorRampPalette(c('blue', 'red'))(length(decay_values))
rrs <- c(0.1, 0.2, 0.4, 0.7, 1, 2, 3)
final_values <- sapply(mc_decays_estimates, function(x) x[2,40])

    emf('lambda_sensitivty_analysis_3july2021.emf', height = 5, width = 7, family = 'Corbel')
    
    par(mar = c(4, 4, 0, 6))
    plot(1, type = 'n', xlim = c(1980, 2019), ylim = range(log(rrs)), axes = F, xlab = NA, ylab = NA)
    segments(1980, 0, 2019)
    for (i in seq_along(decay_values)) {
      pd <- log(mc_decays_estimates[decay_values == decay_values[i]][[1]])
      lines(1980:2019, pd[2,], col = cols[i])
    }
    axis(1, 1980:2019, pos = log(rrs[1]), labels = F, tck = -0.01)
    axis(1, 1980 + 5 * 0:7, pos = log(rrs[1]), tck = -0.03)
    axis(2, log(rrs), rrs, pos = 1980, las = 2)
    rect(1980, log(rrs[1]), 2019, log(max(rrs)))
    title(xlab = 'Year', line = 2.5)
    title(ylab = 'Ratio of new injectors vs. 1980', line = 2.5)
    text(2020, log(final_values), 1/decay_values, col = cols, adj = 0, cex = 0.7)
    text(2020, max(log(final_values)) + 0.2, 'Average duration\nof injecting', adj = 0, cex = 0.7)
    
    dev.off()

#  estimate absolute size of cohorts
#  .................................

# -- existing estimate of number of PWID in 2011, by region. From https://www.drugsandalcohol.ie/21931/1/estimates-of-the-prevalence-of-opiate-use-and-or-crack-cocaine-use-2011-12.pdf #

targets <- data.frame(region = c('London', 'East of England', 'South East', 'North West', 'South West', 'Yorkshire & Humber', 'East Midlands', 'West Midlands', 'North East', 'England'),
                      target = c(6650, 7808, 11351, 6334, 13110, 11047, 10134, 9175, 11692, 87302),
                      stringsAsFactors = F)

# -- function estimating population size each year, given the number starting in 1980, ratios of new starters vs. 1980, and rate of quitting

mm <- function(ratios, new1980 = 10000, decay = 1/15) {
  new15a <- new1980 * ratios
  m <- diag(new15a)
  m <- t(new15a * (1-decay)^(col(m)-row(m))) * (col(m) <= row(m))
  rowSums(m)
}

# -- estimate regional and national cohort sizes with prediction intervals --

find_multiple_cohorts <- function(ratio_matrix, target2011, test_vals = 10 * seq_len(2000)) {
  y <- sapply(test_vals, function(x) {
    pop <- apply(ratio_matrix, 2, mm, new1980 = x)
    pop[1980:2011 == 2011,]
  })
  closest_to_target <- apply(abs(y - target2011), 1, which.min)
  closest_to_target[sapply(closest_to_target, length) == 0] <- NA
  closest_to_target <- unlist(closest_to_target)
  vals1980 <- test_vals[closest_to_target]
  new_cohorts <- rep(vals1980, each = length(1980:2019)) * ratio_matrix
  return(new_cohorts)
}

regions2 <- regions[!(regions %in% c('Northern Ireland', 'Wales'))]
mc_regions2 <- mc_regions[names(mc_regions) %in% regions2]
mc_regions2 <- c(mc_regions2,
                 England = list(mc_decays[abs(round(decay_values - 1/15, 5)) < 0.0001][[1]]))

# -- This code takes a long time to run (depending on the computer about an hour), so the results are saved as an 'Rdata' object and loaded
# mc_cohort_regions2 <- lapply(targets$region, function (x) {
#   print(x)
#   find_multiple_cohorts(mc_regions2[[x]], target2011 = targets[targets$region == x, 'target'])
# })
# save(mc_cohort_regions2, file = 'mc_cohort_regions_4july2021.Rdata')
load(url('https://github.com/danlewer/uam/raw/main/ratios_of_new_initiators/mc_cohort_regions_4july2021.Rdata'))

cohort_size <- lapply(mc_cohort_regions2, function(x) apply(x, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = T))
names(cohort_size) <- names(mc_regions2)
cohort_size <- lapply(cohort_size, t)
cohort_size <- lapply(cohort_size, data.frame)
cohort_size <- lapply(cohort_size, `names<-`, value = c('lower', 'point', 'upper'))

#  -- point estimates --

test_vals <- 10 * seq_len(2000)
point_ratios <- lapply(mc_regions2, function(x) apply(x, 1, median, na.rm = T))
vals1980 <- sapply(targets$region, function(z) test_vals[which.min(abs(sapply(test_vals, function(x) mm(point_ratios[[z]], new1980 = x)[1980:2019 == 2011]) - targets[targets$region == z,'target']))])
point_estimates <- data.table(cbind(year = 1980:2019, round(mapply(`*`, a = point_ratios, b = vals1980), -1)))
for(i in targets$region) {
  cohort_size[[i]]$point <- point_estimates[[i]]
}

#  -- plot of cohort sizes for main article --

cols <- brewer.pal(4, 'Paired')

    emf('cohort_sizes.emf', width = 12, height = 9, family = 'Corbel')
    
    layout(matrix(1:2, ncol = 2), widths = c(2, 1))
    par(mar = c(4, 0, 0, 0), oma = c(0, 5, 0, 4), xpd = NA)
    
    ymax <- 15000
    plot(1, type = 'n', xlim = c(1980, 2019), ylim = c(0, ymax), axes = F, xlab = NA, ylab = NA)
    segments(c(1990, 2000, 2010), 0, y1 = ymax)
    segments(1980, c(5000, 10000), x1 = 2019, lty = 3)
    with(cohort_size$England, {
      polygon(c(1980:2019, 2019:1980), c(lower, rev(upper)), col = cols[1], border = NA)
      lines(1980:2019, point, col = cols[2], lwd = 2)
    })
    axis(1, 1980:2019, pos = 0, labels = F, tck = -0.01)
    axis(1, 1980 + 5 * 0:7, pos = 0, tck = -0.03)
    axis(2, 0:6 * 2500, labels = prettyNum(0:6 * 2500, big.mark = ','), pos = 1980, las = 2)
    rect(1980, 0, 2019, ymax)
    text(1981, ymax * 0.975, 'England', adj = 0)
    title(ylab = 'Number of people starting injecting\nillicit drugs')
    
    panels <- cohort_size[names(cohort_size) != 'England']
    n_panels <- length(panels)
    ymax <- 2000
    yoff <- (seq_along(panels) - 1) * ymax
    ygrids <- 1:3 * 500
    
    plot(1, type = 'n', xlim = c(1980, 2019), ylim = c(0, ymax * n_panels), axes = F, xlab = NA, ylab = NA)
    segments(c(1990, 2000, 2010), 0, y1 = ymax * n_panels)
    for(i in seq_along(panels)) {
      segments(1980, ygrids + yoff[i], 2019, ygrids + yoff[i], lty = 3)
      with(panels[[i]], {
        polygon(c(1980:2019, 2019:1980), c(lower, rev(upper)) + yoff[i], col = cols[1], border = NA)
        lines(1980:2019, point + yoff[i], col = cols[2], lwd = 2)
      })
      rect(1980, yoff[i], 2019, yoff[i] + ymax)
      text(1981, yoff[i] + ymax * 0.9, names(panels)[i], adj = 0, cex = 0.7)
    }
    axis(1, 1980:2019, pos = 0, labels = F, tck = -0.01)
    axis(1, 1980 + 5 * 0:7, pos = 0, tck = -0.03)
    text(x = 2019.5, 
         y = rep(seq_len(n_panels) - 1, each = length(ygrids)) * ymax + rep(ygrids, times = n_panels),
         labels = prettyNum(rep(ygrids, times = n_panels), big.mark = ','),
         adj = 0, cex = 0.7)
    
    dev.off()

# -- table of cohort sizes for supplementary information --

si_absolute <- sapply(cohort_size, function(x) x[,2])

# check regions roughly total England
plot(1980:2019, si_absolute[, 'England'], type = 'l')
lines(1980:2019, rowSums(si_absolute[, -ncol(si_absolute)]), col = 'red')

si_absolute <- cbind(year = 1980:2019, as.data.frame.matrix(si_absolute))

fwrite(si_absolute, 'absolute_number_new_injectors.csv')
