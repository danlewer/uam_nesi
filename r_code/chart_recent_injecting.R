library(data.table)
library(RColorBrewer)

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
