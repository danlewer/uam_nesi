library(data.table)
library(RColorBrewer) # for plot colours
library(devEMF)

# -- named functions --

#  add transparency to colours
add.alpha <- function(cols, alpha) rgb(t(col2rgb(cols)/255), alpha = alpha)

# empty plot
empty_plot <- function(...) plot(1, type = 'n', xlim = c(0, 10), ylim = c(0, 10), axes = F, xlab = NA, ylab = NA, ...)

# panel plot of quantiles
panel_plots <- function (models = reg, actuals = act, REGIONS = regions, parm, YLIM, ndy, COL = brewer.pal(3, 'Pastel1')[1], xlabs = c(13:15), ylabs = c(1, 5, 9, 13), YLAB = 'Age of participants', XLAB = 'Survey year') {
  layout(m = matrix(1:16, ncol = 4, byrow = T))
  par(mar = c(0, 0, 0, 0), oma = c(5, 5, 0, 0), xpd = NA)
  for (i in REGIONS) {
    plot(1, type = 'n', xlim = range(ndy), ylim = YLIM, axes = F, xlab = NA, ylab = NA)
    gridy <- YLIM[1]:YLIM[2]
    gridy <- gridy[gridy %% 5 == 0]
    gridx <- ndy[ndy %% 5 == 0]
    segments(min(ndy), gridy, max(ndy), col = 'grey70', lwd = 0.5)
    segments(gridx, YLIM[1], y1 = YLIM[2], col = 'grey70', lwd = 0.5)
    yq <- models[[i]]
    yq <- yq[year %in% ndy & parameter == parm]
    with(yq, {
      polygon(c(year[tau == 0.5], rev(year[tau == 0.5])), c(fit[tau == 0.25], rev(fit[tau == 0.75])), border = NA, col = COL)
      lines(year[tau == 0.5], fit[tau == 0.5])
      lines(year[tau == 0.5], lower[tau == 0.5], lty = 3)
      lines(year[tau == 0.5], higher[tau == 0.5], lty = 3)
    })
    yq <- actuals[[i]]
    yq <- yq[year %in% ndy & parameter == parm]
    with(yq, {
      points(year, md)
      points(year, lower, pch = 4)
      points(year, upper, pch = 4)
    })
    lab <- if (i %in% c('England', 'Scotland', 'Wales')) '' else 'England: '
    lab <- if (i %in% c('Greater Glasgow & Clyde', 'Lothian & Tayside', 'Other regions of Scotland')) 'Scotland: ' else lab
    lab <- paste0(lab, i)
    text(mean(ndy), YLIM[2] * 0.95, lab)
    # using 'segments' here rather than 'rect' to ensure correct rendering in PDF
    segments(min(ndy), YLIM[1], x1 = max(ndy))
    segments(min(ndy), YLIM[2], x1 = max(ndy))
    segments(min(ndy), YLIM[1], y1 = YLIM[2])
    segments(max(ndy), YLIM[1], y1 = YLIM[2])
    if (match(i, regions) %in% ylabs) axis(2, gridy, pos = min(ndy), las = 2)
    if (match(i, regions) %in% xlabs) axis(1, gridx, pos = YLIM[1])
  }
  empty_plot()
  ys <- seq(1, 8, length.out = 4)
  segments(x0 = 0, y0 = ys[4] + c(-0.5, 0, 0.5), x1 = 2, lty = c(3, 1, 3))
  rect(0, ys[3] - 1, 2, ys[3] + 1, col = COL, border = NA)
  points(1, ys[2])
  points(1, ys[1], pch = 4)
  text(2.25, ys, c('Sample IQR', 'Sample median', 'Modelled IQR', 'Modelled median\n(95% CI)'), adj = 0)
  if (!is.na(YLAB)) {mtext(YLAB, side = 2, outer = T, line = 3, cex = 0.8)}
  if (!is.na(XLAB)) {mtext(XLAB, side = 1, outer = T, line = 3, cex = 0.8)}
}

#  :::::::::
#  read data
#  .........

act_uam <- read.csv(url("https://raw.githubusercontent.com/danlewer/uam_nesi/main/quantiles/uam_actual_quantiles_19jan2022.csv"), stringsAsFactors = F)
reg_uam <- read.csv(url("https://raw.githubusercontent.com/danlewer/uam_nesi/main/quantiles/uam_modelled_quantiles_19jan2022.csv"), stringsAsFactors = F)
act_nesi <- read.csv(url("https://raw.githubusercontent.com/danlewer/uam_nesi/main/quantiles/nesi_actual_quantiles_2oct2021.csv"), stringsAsFactors = F)
reg_nesi <- read.csv(url("https://raw.githubusercontent.com/danlewer/uam_nesi/main/quantiles/nesi_modelled_quantiles_2oct2021.csv"), stringsAsFactors = F)
setDT(act_uam); setDT(reg_uam); setDT(act_nesi); setDT(reg_nesi)
names(act_nesi)[1] <- 'year'
names(reg_nesi)[1] <- 'year'

act <- rbind(act_uam, act_nesi)
reg <- rbind(reg_uam, reg_nesi)

act$region[act$region == 'GCC'] <- 'Greater Glasgow & Clyde'
act$region[act$region == 'Lothian_and_Tayside'] <- 'Lothian & Tayside'
act$region[act$region == 'Other'] <- 'Other regions of Scotland'

reg$region[reg$region == 'GCC'] <- 'Greater Glasgow & Clyde'
reg$region[reg$region == 'Lothian_and_Tayside'] <- 'Lothian & Tayside'
reg$region[reg$region == 'Other'] <- 'Other regions of Scotland'

regions <- c(c('England', 'Scotland', 'Wales'), setdiff(unique(act$region), c('England', 'Scotland', 'Wales')))
regions <- setdiff(regions, 'Northern Ireland')

act <- split(act, f = act$region)
reg <- split(reg, f = reg$region)

#  :::::::::::::::::::::::::::::::::::::::::::::::::::::
#  panel plots by region (for Supplementary Information)
#  .....................................................

emf('age_by_region_panel.emf', height = 8, width = 8, family = 'Candara')
panel_plots(parm = 'age', YLIM = c(20, 50), ndy = 1990:2019)
dev.off()

emf('age_initiation_panel.emf', height = 8, width = 8, family = 'Candara')
panel_plots(parm = 'age_init', YLIM = c(15, 45), ndy = 1990:2019, YLAB = 'Age first injected', XLAB = 'Year first injected')
dev.off()

emf('duration_injecting_by_region_panel.emf', height = 8, width = 8, family = 'Candara')
panel_plots(parm = 'dur', YLIM = c(0, 25), ndy = 1990:2019, YLAB = 'Duration of injecting')
dev.off()

#  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#  effect of limiting to those initiated in past 3 years on age of initiation
#  ..........................................................................

# this is a sensitivity analysis that was reported during peer review but is not included in the published article

sampleInit <- dcast(act$England[parameter %in% c('age_init_whole_sample', 'age_init') & year >= 1990], year ~ parameter, value.var = 'md')
sampleInit[, age_init := age_init + 0.1]
sampleInit[, age_init_whole_sample - 0.1]
modelInit <- dcast(reg$England[parameter %in% c('age_init_whole_sample', 'age_init') & tau == 0.5], year ~ parameter, value.var = 'fit')

cols <- brewer.pal(3, 'Set1')

emf('age_of_initiation_whole_vs_past3year.emf', height = 5, width = 6, family = 'Candara')

plot(1, type = 'n', xlim = c(1990, 2019), ylim = c(18, 35), xlab = NA, ylab = NA, axes = F)
rect(1989, 18, 2020, 35)
with(modelInit, {
  lines(year, age_init, col = cols[1])
  lines(year, age_init_whole_sample, col = cols[2])
})
with(sampleInit, {
  points(year, age_init, pch = 4, col = cols[1])
  points(year, age_init_whole_sample, pch = 4, col = cols[2])
})
axis(1, seq(1990, 2015, 5), pos = 18)
axis(2, c(18, 20, 25, 30, 35), pos = 1989, las = 2)
title(xlab = 'Year of initiation')
title(ylab = 'Age of initiation')

dev.off()

#  :::::::::::::::::::::::::::::::::::
#  miscellaneous values for manuscript
#  ...................................

# median age in 1990, 2008, and 2019
dcast(reg$England[parameter == 'age' & year %in% c(1990, 2008, 2019)], year ~ tau, value.var = 'fit')
dcast(reg$Scotland[parameter == 'age' & year %in% c(1990, 2008, 2019)], year ~ tau, value.var = 'fit')
rbindlist(reg)[year == 2019 & parameter == 'age' & tau == 0.5][order(fit)]
dcast(rbindlist(reg)[year == 2019 & parameter == 'age'], region ~ tau, value.var = 'fit')[order(`0.5`)]

# median age of initiation in 1990, 2008, and 2019
dcast(reg$England[parameter == 'age_init' & year %in% c(1990, 2008, 2019)], year ~ tau, value.var = 'fit')
dcast(reg$Scotland[parameter == 'age_init' & year %in% c(1990, 2008, 2019)], year ~ tau, value.var = 'fit')
rbindlist(reg)[year == 2019 & parameter == 'age_init' & tau == 0.5][order(fit)]

# median duration of injecting in 1990, 2008, and 2019
dcast(reg$England[parameter == 'dur' & year %in% c(1990, 2008, 2019)], year ~ tau, value.var = 'fit')
dcast(reg$Scotland[parameter == 'dur' & year %in% c(1990, 2008, 2019)], year ~ tau, value.var = 'fit')
rbindlist(reg)[year == 2019 & parameter == 'dur' & tau == 0.5][order(fit)]


#  :::::::::::::::::::::::::::::::::::::
#  histograms of duration by survey year
#  .....................................

# -- England --

ndy <- 1993:2019
id <- read.csv(url("https://raw.githubusercontent.com/danlewer/uam_nesi/main/histogram/uam_histogram_dur_2oct2021.csv"), stringsAsFactors = F)
iy <- read.csv(url("https://raw.githubusercontent.com/danlewer/uam_nesi/main/histogram/uam_histogram_yearInit_2oct2021.csv"), stringsAsFactors = F)
id <- id[id$year %in% ndy,]
iy <- iy[iy$year %in% ndy,]
id <- split(id, f = id$year)
iy <- split(iy, f = iy$year)

cols <- colorRampPalette(brewer.pal(11, 'Spectral'))(length(ndy))

cairo_pdf('Figure3.pdf', height = 6, width = 12, family = 'Candara')

par(mfrow = c(1, 2), mar = c(0, 0, 0, 0), oma = c(4, 4, 0, 7), xpd = NA)

plot(1, type = 'n', xlim = c(0, 40), ylim = c(0, 0.1), axes = F, xlab = NA, ylab = NA)
rect(0, 0, 40, 0.1)
lapply(seq_along(ndy), function (x) {
  with(id[[x]], {
    points(injdurI, pc, col = cols[x], pch = 19, cex = 0.7)
    lines(injdurI, pc, col = cols[x])
  })
})
axis(1, 0:40, labels = F, tck = -0.01, pos = 0)
axis(1, 0:7 * 5, tck = -0.03, pos = 0)
axis(2, 0:10/100, paste0(0:10, '%'), las = 2, pos = 0)
title(xlab = 'Duration of injecting', line = 2.5)
title(ylab = 'Proportion of participants', line = 2.5)

plot(1, type = 'n', xlim = c(1970, 2020), ylim = c(0, 0.1), axes = F, xlab = NA, ylab = NA)
rect(1970, 0, 2020, 0.1)
lapply(seq_along(iy), function (x) {
  with(iy[[x]], {
    points(fyearI, pc, col = cols[x], pch = 19, cex = 0.7)
    lines(fyearI, pc, col = cols[x])
  })
})
axis(1, 1970 + 0:4 * 10, tck = -0.03, pos = 0)
axis(1, 1970:2019, labels = F, tck = -0.01, pos = 0)
title(xlab = 'Year of initiation', line = 2.5)

ys <- seq(0, 0.09, length.out = length(iy))
segments(2022, ys, 2026, col = cols)
points(rep(2024, length(ys)), ys, col = cols, pch = 19)
text(2027, ys[1993:2019 %% 5 == 0], (1993:2019)[1993:2019 %% 5 == 0], adj = 0)
text(2022, max(ys) * 1.05, 'Survey year', adj = 0)

dev.off()

# -- Scotland --

ndy <- 2008:2019
id <- read.csv(url("https://raw.githubusercontent.com/danlewer/uam_nesi/main/histogram/nesi_histogram_dur_2oct2021.csv"), stringsAsFactors = F, col.names = c('year', 'injdurI', 'pc'))
iy <- read.csv(url("https://raw.githubusercontent.com/danlewer/uam_nesi/main/histogram/nesi_histogram_yearInit_2oct2021.csv"), stringsAsFactors = F, col.names = c('year', 'fyearI', 'pc'))
id <- id[id$year %in% ndy,]
iy <- iy[iy$year %in% ndy,]
id <- split(id, f = id$year)
iy <- split(iy, f = iy$year)

cols <- colorRampPalette(brewer.pal(11, 'Spectral'))(length(ndy))

emf('nesi_duration_injecting_histograms_22jan2022.emf', height = 6, width = 12, family = 'Candara')

par(mfrow = c(1, 2), mar = c(0, 0, 0, 0), oma = c(4, 4, 0, 7), xpd = NA)

plot(1, type = 'n', xlim = c(0, 40), ylim = c(0, 0.1), axes = F, xlab = NA, ylab = NA)
rect(0, 0, 40, 0.1)
lapply(seq_along(ndy), function (x) {
  with(id[[x]], {
    points(injdurI, pc, col = cols[x], pch = 19, cex = 0.7)
    lines(injdurI, pc, col = cols[x])
  })
})
axis(1, 0:40, labels = F, tck = -0.01, pos = 0)
axis(1, 0:7 * 5, tck = -0.03, pos = 0)
axis(2, 0:10/100, paste0(0:10, '%'), las = 2, pos = 0)
title(xlab = 'Duration of injecting', line = 2.5)
title(ylab = 'Proportion of participants', line = 2.5)

plot(1, type = 'n', xlim = c(1970, 2020), ylim = c(0, 0.1), axes = F, xlab = NA, ylab = NA)
rect(1970, 0, 2020, 0.1)
lapply(seq_along(iy), function (x) {
  with(iy[[x]], {
    points(fyearI, pc, col = cols[x], pch = 19, cex = 0.7)
    lines(fyearI, pc, col = cols[x])
  })
})
axis(1, 1970 + 0:4 * 10, tck = -0.03, pos = 0)
axis(1, 1970:2019, labels = F, tck = -0.01, pos = 0)
title(xlab = 'Year of initiation', line = 2.5)

ys <- seq(0, 0.09, length.out = length(iy))
segments(2022, ys, 2026, col = cols)
points(rep(2024, length(ys)), ys, col = cols, pch = 19)
text(2027, ys, ndy, adj = 0)
text(2022, max(ys) * 1.07, 'Survey year', adj = 0)

dev.off()

#  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#  single plot with age of initiation, duration of injecting, and age of population
#  ................................................................................

pf <- function(actual = act, modelled = reg, parm = 'age', xlim = c(1993, 2019), ylim = c(20, 50), EngCol = 'red', ScotCol = 'blue', title = NA, alpha = 0.7) {
  plot(1, type = 'n', axes = F, xlab = NA, ylab = NA, xlim = xlim, ylim = ylim)
  rect(xlim[1], ylim[1], xlim[2], ylim[2], lwd = 0.5)
  yrs <- xlim[1]:xlim[2]
  axis(1, yrs, pos = ylim[1], labels = F, tck = -0.01, lwd = 0.5)
  axis(1, yrs[yrs %% 5 == 0], pos = ylim[1], tck = -0.03, lwd = 0.5)
  yl <- ylim[1]:ylim[2]
  axis(2, yl[yl %% 5 == 0], pos = xlim[1], las = 2, tck = -0.03, lwd = 0.5)
  axis(2, yl, labels = F, pos = xlim[1], tck = -0.01, lwd = 0.5)
  segments(xlim[1], yl[yl %% 5 == 0], x1 = xlim[2], lty = 3, lwd = 0.5)
  with(modelled$England[year %in% yrs & parameter == parm], {
    yrs2 <- unique(year)
    polygon(x = c(yrs2, rev(yrs2)), y = c(fit[tau == 0.25], rev(fit[tau == 0.75])), col = add.alpha(EngCol, alpha = alpha), border = NA)
    lines(yrs, fit[tau == 0.5], col = EngCol, lwd = 1)
  })
  with(modelled$Scotland[year %in% yrs & parameter == parm], {
    yrs2 <- unique(year)
    polygon(x = c(yrs2, rev(yrs2)), y = c(fit[tau == 0.25], rev(fit[tau == 0.75])), col = add.alpha(ScotCol, alpha = alpha), border = NA)
    lines(yrs2, fit[tau == 0.5], col = ScotCol, lwd = 1)
  })
  with(actual$England[year %in% yrs & parameter == parm], points(year, md, col = EngCol, pch = 4))
  with(actual$Scotland[year %in% yrs & parameter == parm], points(year, md, col = ScotCol, pch = 4))
  text(mean(xlim), (ylim[2] - ylim[1]) * 1.05 + ylim[1], title, font = 2)
}

cols <- brewer.pal(3, 'Set1')[2:3]

cairo_pdf('Figure2.pdf', height = 3.5, width = 8, family = 'Corbel')

layout(m = matrix(1:4, ncol = 4), widths = c(3, 3, 3, 3))
par(mar = c(3, 2, 3, 0), xpd = NA)
pf(EngCol = cols[1], ScotCol = cols[2], alpha = 0.3, title = 'Age of participants', ylim = c(15, 50))
pf(parm = 'age_init', alpha = 0.3, EngCol = cols[1], ScotCol = cols[2], title = 'Age of initiation', ylim = c(15, 50))
pf(parm = 'dur', EngCol = cols[1], ScotCol = cols[2], alpha = 0.3, title = 'Duration of injecting', ylim = c(0, 25))

par(mar = c(0, 0, 0, 0))
empty_plot()
ys <- seq(1, 9, length.out = 10)
rect(1, ys[7], 3, ys[9], col = add.alpha(cols[1], 0.3), border = NA)
segments(1, ys[8], 3, ys[8], col = cols[1], lwd = 1)
points(2, ys[6], col = cols[1], pch = 4)
xlab <- c('Actual median', 'Lower quartile', 'Median', 'Upper quartile')
text(4, ys[6:9], xlab, adj = 0)
text(1, ys[10], 'England, Wales\nand Northern Ireland', adj = 0, font = 2)
rect(1, ys[2], 3, ys[4], col = add.alpha(cols[2], 0.3), border = NA)
segments(1, ys[3], 3, ys[3], col = cols[2], lwd = 1)
points(2, ys[1], col = cols[2], pch = 4)
text(4, ys[1:4], xlab, adj = 0)
text(1, ys[5], 'Scotland', adj = 0, font = 2)

dev.off()

#  ::::::::::::::::::::::::::::::::::::::::::::::::::
#  sensitivity analysis by averaging regional results
#  ..................................................

# this is a sensitivity analysis that was reported during peer review but is not included in the published article

parms <- c('age', 'age_init', 'dur')

sens_regions <- reg_uam[!(region %in% c('England', 'Wales', 'Northern Ireland')) & tau == 0.5 & parameter != 'age_init_whole_sample', .(mean_reg = mean(fit)), c('year', 'parameter')][
  reg_uam[region == 'England' & tau == 0.5 & parameter != 'age_init_whole_sample', c('year', 'fit', 'parameter')], on = c('year', 'parameter')]
sens_regions[, parameter := factor(parameter, c('age', 'age_init', 'dur'), c('Age of participants', 'Age of initiation', 'Duration of injecting'))]
sens_regions <- split(sens_regions, f = sens_regions$parameter)

cols <- brewer.pal(3, 'Set1')

emf('sensitivity_averaging_regions.emf', height = 4, width = 9, family = 'Candara')

par(mfrow = c(1, 3))
lapply(sens_regions, function (x) {
  with(x, {
    plot(1, type = 'n', xlim = c(1990, 2019), ylim = c(0, 40), xlab = NA, ylab = NA)
    lines(year, fit, col = cols[1])
    lines(year, mean_reg, col = cols[2])
    title(main = x$parameter[1])
  })
})

dev.off()
