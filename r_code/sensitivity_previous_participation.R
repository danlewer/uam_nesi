library(RColorBrewer)
library(data.table)

load(url('https://github.com/danlewer/uam_nesi/raw/main/misc/previousParticipationSensitivities_27jan2022.Rdata'))

cols <- brewer.pal(3, 'Set1')

plotPrevPartSens <- function (p, col = cols, ylim = c(0, 50), xlim = c(1990, 2019), main = 'Age', ...) {
  plot(1, type = 'n', xlim = xlim, ylim = ylim, axes = F, xlab = NA, ylab = NA, ...)
  rect(xlim[1], ylim[1], xlim[2], ylim[2])
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
