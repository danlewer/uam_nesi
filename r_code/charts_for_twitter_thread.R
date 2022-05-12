library(data.table)
library(RColorBrewer) # for plot colours

ndy <- 1993:2019
id <- read.csv(url("https://raw.githubusercontent.com/danlewer/uam_nesi/main/histogram/uam_histogram_dur_2oct2021.csv"), stringsAsFactors = F)
iy <- read.csv(url("https://raw.githubusercontent.com/danlewer/uam_nesi/main/histogram/uam_histogram_yearInit_2oct2021.csv"), stringsAsFactors = F)
id <- id[id$year %in% ndy,]
iy <- iy[iy$year %in% ndy,]
id <- split(id, f = id$year)
iy <- split(iy, f = iy$year)

cols <- colorRampPalette(brewer.pal(11, 'Spectral'))(length(ndy))

png('addiction_tweet3.png', height = 5, width = 5, units = 'in', res = 300)

par(xpd = NA)
plot(1, type = 'n', xlim = c(1970, 2020), ylim = c(0, 0.06), axes = F, xlab = NA, ylab = NA)
rect(1970, 0, 2020, 0.06)
with(iy$`2019`, {
  points(fyearI, pc, col = cols[length(cols)], pch = 19, cex = 0.7)
  lines(fyearI, pc, col = cols[length(cols)])
})
axis(1, 1970 + 0:4 * 10, tck = -0.03, pos = 0)
axis(1, 1970:2019, labels = F, tck = -0.01, pos = 0)
title(xlab = 'Year of initiation', line = 2.5)
title(ylab = 'Proportion of participants', line = 2.5)
axis(2, 0:6/100, paste0(0:6, '%'), las = 2, pos = 1970)
text(1970, 0.065, 'Year of initiation of injecting for UAM\nparticipants in 2019', adj = 0)

dev.off()

png('addiction_tweet4.png', height = 5, width = 5, units = 'in', res = 300)

par(xpd = NA)
plot(1, type = 'n', xlim = c(1980, 2020), ylim = c(0, 0.06), axes = F, xlab = NA, ylab = NA)
rect(1980, 0, 2020, 0.06)
pcs2019 <- iy$`2019`$pc[iy$`2019`$fyearI >= 1980]
rat <- exp(-(0:(length(pcs2019)-1))/15)
pcs2019rescale <- pcs2019 / rev(rat)
pcs2019rescale <- pcs2019rescale / sum(pcs2019rescale)
lines(1980:2019, pcs2019, col = cols[length(cols)])
points(1980:2019, pcs2019, col = cols[length(cols)], pch = 19, cex = 0.7)
lines(1980:2019, pcs2019rescale, col = cols[1], lty = 2)
points(1980:2019, pcs2019rescale, col = cols[1], pch = 19, cex = 0.7)
axis(1, 1980 + 0:3 * 10, tck = -0.03, pos = 0)
axis(1, 1980:2019, labels = F, tck = -0.01, pos = 0)
title(xlab = 'Year of initiation', line = 2.5)
title(ylab = 'Proportion of participants', line = 2.5)
axis(2, 0:6/100, paste0(0:6, '%'), las = 2, pos = 1980)
text(1980, 0.065, 'Year of initiation of injecting for UAM\nparticipants in 2019', adj = 0)
text(1990, 0.055, 'Adjusted for\nquitting', col = cols[1])
text(2010, 0.04, 'Crude', col = cols[length(cols)])

dev.off()


png('addiction_tweet5.png', height = 7, width = 8, units = 'in', res = 300)

par(mar = c(5, 5, 3, 7), xpd = NA)
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
axis(2, 0:10/100, paste0(0:10, '%'), las = 2, pos = 1970)

ys <- seq(0, 0.09, length.out = length(iy))
segments(2022, ys, 2026, col = cols)
points(rep(2024, length(ys)), ys, col = cols, pch = 19)
text(2027, ys[1993:2019 %% 5 == 0], (1993:2019)[1993:2019 %% 5 == 0], adj = 0)
text(2022, max(ys) * 1.05, 'Survey year', adj = 0)
title(ylab = 'Proportion of participants', line = 2.5)
text(1970, 0.105, 'Year of initiation of injecting for UAM\nparticipants', adj = 0)

dev.off()
