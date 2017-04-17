loo_co2 <- read.table("output/nwgksloo_CO2.dat")
loo_cos <- read.table("output/nwgksloo_cosine.dat")
loo_dif <- read.table("output/nwgksloo_diffract.dat")
loo_js <- read.table("output/nwgksloo_JS.dat")
loo_mul <- read.table("output/nwgksloo_multicos.dat")
loo_rpd <- read.table("output/nwgksloo_RPDI.dat")
loo_rvn <- read.table("output/nwgksloo_RVN63.dat")
loo_zig <- read.table("output/nwgksloo_zigzag.dat")

par(las=1, mar=c(5,4,1,1), mfrow=c(4,2))

plo <- function(loo, lab)
{
	plot(loo, xlab=expression(lambda), ylab="LOO RMS", log="xy", type="l")
	grid()
	y <- loo[is.finite(loo[,2]), 2]

	# Label the plot directly.
	text(quantile(loo[,1], 0.13), min(y) + 0.22 * diff(range(y)),
	     lab, cex=1.2)

	# Draw a vertical dotted line representing `nwgks`'s initial
	# `lambda` guess.
	abline(v=prod(range(loo[is.finite(loo[,2]), 1])^c(0.75, 0.25)),
	       lty="dotted")
}

plo(loo_cos, "synthetic\ncosine")
plo(loo_zig, "synthetic\nzigzag")
plo(loo_mul, "synthetic\nmulti-cosine")
plo(loo_co2, "CO2")
plo(loo_dif, "diffraction")
plo(loo_js, "Jill\nStein")
plo(loo_rpd, "UK RPDI")
plo(loo_rvn, "RVN '63\nanthro.")

par(mfrow=c(1,1))
