library(devEMF)

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
