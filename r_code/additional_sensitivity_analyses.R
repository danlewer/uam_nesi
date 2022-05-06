library(devEMF)
library(data.table)
library(RColorBrewer)

#  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#  percent of participants who have previously participated
#  ........................................................

# this analysis was reported during peer review but is not included in the main article

pp <- read.csv(url('https://raw.githubusercontent.com/danlewer/uam_nesi/main/misc/previous_participation.csv'))
pp$pc <- pp$reporting_previous_participation / pp$participants

emf('percent_participation_previously.emf', height = 5, width = 5, family = 'Candara')

plot(1, type = 'n', xlim = c(1990, 2019), ylim = c(0, 0.6), axes = F, xlab = 'Survey year', ylab = 'Percent reporting previous participation')
rect(1990, 0, 2019, 0.6, col = 'grey97')
with(pp, {
  points(survey_year, pc, pch = 19)
  lines(survey_year, pc)
})
axis(1, seq(1990, 2015, 5), pos = 0)
axis(2, 0:6/10, 0:6 * 10, pos = 1990, las = 2)

dev.off()

#  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#  percent reporting recent injecting (for Supplementary Information)
#  ..................................................................

recent_injecting <- read.csv(url('https://raw.githubusercontent.com/danlewer/uam_nesi/main/misc/proportion_injected_past_year.csv'))
setDT(recent_injecting)

top <- t(apply(recent_injecting[, -'year'], 1, cumsum)) / rowSums(recent_injecting[, -'year'])
colnames(top) <- c('top_yes', 'top_no', 'top_missing')
recent_injecting <- cbind(recent_injecting, top)

recent_plot <- melt(recent_injecting, id.vars = 'year', measure.vars = c('top_yes', 'top_no', 'top_missing'), variable.name = 'recent', value.name = 'top')
recent_plot[, recent := gsub('top_', '', recent)]
recent_plot <- recent_plot[order(year)]
recent_plot[, bottom := shift(top, fill = 0), year]
cols <- brewer.pal(3, 'Pastel1')
recent_plot[, col := as.character(factor(recent, c('yes', 'no', 'missing'), cols))]
xlab <- seq(1990, 2015, 5)

emf('injected_past_year_19jan2022.emf', height = 3.5, width = 6, family = 'Candara')

par(mar = c(5, 4, 0, 7), xpd = NA)
plot(1, type = 'n', xlim = c(1990, 2020), ylim = c(0, 1), axes = F, xlab = NA, ylab = NA)
with(recent_plot, rect(year, bottom, year + 1, top, col = col))
axis(1, 1990:2020, labels = F, pos = 0, tck = -0.03)
axis(1, xlab, labels = F, pos = 0, tck = -0.06)
axis(1, xlab + 2.5, paste0(xlab, '\n', '-', xlab+4), tick = F, xlab, pos = -0.03)
axis(2, 0:5/5, paste0(0:5 * 20, '%'), las = 2, pos = 1990)
title(xlab = 'Survey year')
title(ylab = 'Percent of participants')
ys <- seq(0.3, 0.7, length.out = 4)
rect(2021, ys[-length(ys)], 2022, ys[-1], col = cols)
text(2022.5, ys[-length(ys)] + diff(ys)/2, c('Yes', 'No', 'Missing'), adj = 0)
text(2021, max(ys) + 0.12, 'Injected in\npast year', adj = 0)

dev.off()

#  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#  comparison of results with and without those reporting previous survey participation
#  ....................................................................................

# this analysis was reported during peer review but is not included in the main article

load(url('https://github.com/danlewer/uam_nesi/raw/main/misc/previousParticipationSensitivities_27jan2022.Rdata'))

cols <- brewer.pal(3, 'Set1')

plotPrevPartSens <- function (p, col = cols, ylim = c(0, 50), xlim = c(1990, 2019), main = 'Age', ...) {
  plot(1, type = 'n', xlim = xlim, ylim = ylim, axes = F, xlab = NA, ylab = NA, ...)
  # segments used here instead of 'rect' to ensure proper rendering in PDF
  segments(xlim[1], ylim[1], x1 = xlim[2])
  segments(xlim[1], ylim[2], x1 = xlim[2])
  segments(xlim[2], ylim[1], y1 = ylim[2])
  axis(1, seq(1990, 2015, 5), pos = ylim[1])
  with(p[sens == 'all'], {
    points(year, sample, col = cols[1])
    lines(year, fit, col = cols[1])
  })
  with(p[sens == 'new'], {
    points(year - 0.1, sample, col = cols[2])
    lines(year, fit, col = cols[2])
  })
  with(p[sens == 'prev'], {
    points(year + 0.1, sample, col = cols[3])
    lines(year, fit, col = cols[3])
  })
  text(1990, ylim[2] * 1.05, main, adj = 0)
}

emf('repeat_participation_sens_27jan2022.emf', height = 5, width = 10, family = 'Candara')

par(mfrow = c(1, 3), oma = c(5, 2, 5, 15), mar = c(0, 1, 0, 0), xpd = NA)
plotPrevPartSens(ageSpitSens, ylim = c(0, 40))
axis(2, seq(0, 40, 5), pos = 1990, las = 2)
title(xlab = 'Survey year')
plotPrevPartSens(durSpitSens, ylim = c(0, 40), main = 'Duration')
axis(2, seq(0, 40, 5), pos = 1990, las = 2, labels = NA)
title(xlab = 'Survey year')
plotPrevPartSens(iniSpitSens, ylim = c(0, 40), main = 'Age at initiation')
axis(2, seq(0, 40, 5), pos = 1990, las = 2, labels = NA)
title(xlab = 'Year of initiation')
ys <- seq(27, 35, length.out = 3)
segments(2021, ys, 2025, ys, col = cols)
text(2026, ys, c('All (main results)', 'New participants\nonly', 'Repeat participants\nonly'), adj = 0)

dev.off()
