# Take `hammy` output and compare it against the output of R's standard GAM
# package `mgcv`. This script compares the partial residuals and
# individual smooth functions estimated by `hammy` (blue & red) to those
# estimated by `mgcv`'s `gam` function (grey & black).

library(mgcv)

h <- read.table("out-hammy.dat")

k <- ncol(h) / 3

draw_hammy_term_with_error_curves <- function(h, select)
{
	h <- h[order(h[, select]),]
	idx <- k + 3 + (2 * (select - 1))
	y <- h[, idx]
	lines(h[, select], y, col="red")
	lines(h[, select], y - (2 * sqrt(h[, 1 + idx])), col="red", lty="dotted")
	lines(h[, select], y + (2 * sqrt(h[, 1 + idx])), col="red", lty="dotted")
	grid(col="#0000001f")
}

gam_mod <- gam(as.formula(paste(colnames(h)[k],
                                "~",
                                paste(paste0("s(", colnames(h)[1:(k-1)], ")"),
                                      collapse=" + "))),
               gaussian, h[, 1:k])

#par(las=1, mar=c(4.9, 4, 0.1, 0.1), mfrow=c(2,3))
par(las=1, mar=c(4.9, 4, 0.1, 0.1), mfrow=c(1,1))

for (i in 1:(k-1)) {
	plot(gam_mod, select=i, col="#000000df")
	partial_residuals <- residuals(gam_mod, type="working") +
	                     predict(gam_mod, type="terms")[,i]
	points(h[,i], partial_residuals, pch=21, col="#0000004f", bg="#00000004")
#	points(jitter(h[,i], 0.7), h[,k] - mean(h[,k]), pch=4, cex=0.5,
#	       col="#00bfff7a")
	points(jitter(h[,i], 0.7),
	       h[,k] - mean(h[,k]) - apply(h[, seq(k+3, ncol(h), 2)], 1, sum)
	             + h[, k + 3 + (2 * (i - 1))],
	       pch=4, cex=0.7, col="#00bfff7a")
	draw_hammy_term_with_error_curves(h, i)
	readline("Press Enter to see next plot.")
}
