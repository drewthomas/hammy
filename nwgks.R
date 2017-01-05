library(mgcv)

d <- read.table("output/nwgks_multicos.dat")
d <- d[order(d$V1),]

draw_nwgks_fit_with_error_curves <- function(d, mean_shift=FALSE)
{
	y <- d$V3
	if (mean_shift) {
		y <- y - mean(d$V2)
	}
	lines(d$V1, y, col="red")
	lines(d$V1, y - sqrt(d$V4), col="red", lty="dotted")
	lines(d$V1, y + sqrt(d$V4), col="red", lty="dotted")
	lines(d$V1, y - sd(d$V3 - d$V2), col="blue", lty="dotted")
	lines(d$V1, y + sd(d$V3 - d$V2), col="blue", lty="dotted")
}

par(las=1, mar=c(4.9, 4, 0.1, 0.1), mfrow=c(2,1))

if ((length(d$V1) / length(unique(d$V1))) > 60) {
	plot(jitter(d$V1, 0.7), d$V2, pch=21, cex=0.5,
	     col="#000000aa", bg="#0000001a")
} else {
	plot(d$V1, d$V2, cex=0.5)
}
grid()
draw_nwgks_fit_with_error_curves(d)

plot(gam(V2 ~ s(V1), gaussian, d), seWithMean=TRUE)
#plot(gam(V2 ~ s(V1), poisson(link="identity"), d), seWithMean=TRUE)
points(d$V1, d$V2 - mean(d$V2), pch=21, cex=0.5,
       col="#0000004f", bg="#0000001f")
draw_nwgks_fit_with_error_curves(d, TRUE)
