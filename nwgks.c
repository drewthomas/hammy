/* Nadaraya-Watson Gaussian-kernel smoothing */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "data.h"
#include "nwgks.h"

/* Relative tightness of the quick bounds on lambda. Larger values relax the
   bounds, risking longer runtime but making it (even) less likely that the
   RMS-minimizing lambda lies outside the quick bounds. */
#define REL_BOUND_LOOSENESS 70.0

/* Relative bound-based tolerance to use to decide to stop iteratively
   converging on the LOO RMS-minimizing lambda. */
#define NWGKS_REL_TOL_BO 2e-9

/* Relative size of the iterative RMS minimizer's first step away from the
   initial lambda guess. */
#define REL_INITIAL_STEP 4e-3

/* Parabolic interpolation requires three samples of the RMS LOO function
   that are reasonably close to the local minimum. When samples become
   available for which the upper abscissa divided by the lower abscissa is
   at most this number, use parabolic interpolation instead of trial-and-
   error stepping. */
#define PI_TRIGGER_RATIO 2.0

/* Compute the mean of `x`, an array of `double`s of length `n`. */
double mean_of_array_of_doubles(const double* x, unsigned int n)
{
	unsigned int i;
	double m;

	m = x[0];
	for (i = 1; i < n; i++) {
		m = ((i / (double) (1 + i)) * m) + (x[i] / (1 + i));
	}

	return m;
}

/* Using the data `b`, compute a Nadaraya-Watson Gaussian-kernel estimate
   of the ordinate at the abscissa `at` using the smoothing parameter
   `lambda`, higher values of which produce more aggressive smoothing. */
double nwgks_at(const Bivariate* b, double at, double lambda)
{
	unsigned int i;
	double sum_w = 0.0;
	double w;
	double y = 0.0;

	if (isinf(lambda)) {
		/* In the limit of infinite smoothing, everything gets a weight
		   of exp(0) = 1, and hence the estimate is an unweighted mean. */
		return mean_of_array_of_doubles(b->y, b->n);
	}

	for (i = 0; i < b->n; i++) {
		w = (b->x[i] - at) / lambda;
		w = exp(-w * w / 2.0);
		sum_w += w;
		y += w * b->y[i];
	}

	return y / sum_w;
}

/* Using the data `b`, compute a NWGK estimate of the ordinate at the
   abscissa `at` using the smoothing parameter `lambda`, storing an estimate
   of the ordinate's variance in `estimator_variance` (computed from the
   error-variance estimate `error_variance`). */
double nwgks_at_with_var(const Bivariate* b, double at, double lambda, double* estimator_variance, double error_variance)
{
	unsigned int i;
	double sum_w = 0.0;
	double w;
	double y = 0.0;

	*estimator_variance = 0.0;

	for (i = 0; i < b->n; i++) {
		w = exp(-(b->x[i] - at) * (b->x[i] - at) / (2.0 * lambda * lambda));
		sum_w += w;
		y += w * b->y[i];
		*estimator_variance += w * w * error_variance;
	}

	*estimator_variance /= (sum_w * sum_w);
	return y / sum_w;
}

/* Using the data `b`, compute a NWGK estimate of the ordinate at the
   abscissa `at`, leaving out the datum with index `leave_out_idx`.
   NB: this function may return NaN if `lambda` is so low that all of the
   retained data points have weights represented as 0 in FP arithmetic.
   The caller is responsible for checking this. */
double nwgks_at_leaving_out(const Bivariate* b, double at, double lambda, unsigned int leave_out_idx)
{
	unsigned int i;
	double sum_w = 0.0;
	double w;
	double y = 0.0;

	for (i = 0; i < b->n; i++) {
		if (i == leave_out_idx) {
			continue;
		}
		w = exp(-(b->x[i] - at) * (b->x[i] - at) / (2.0 * lambda * lambda));
		sum_w += w;
		y += w * b->y[i];
	}

	return y / sum_w;
}

/* Compute the leave-one-out root-mean-square residual of an NWG kernel
   smooth of `b` with smoothing parameter `lambda`. Specifically, run through
   each datum, estimating the datum's ordinate (y-value) using the other
   data points, and taking the RMS error of that estimate over the data. */
double nwgks_loo_rms(const Bivariate* b, double lambda)
{
	unsigned int i;
	double mean_of_y;
	double residual;
	double ss = 0.0;

	if (isinf(lambda)) {
		/* Infinite smoothing implies that everything gets a weight of 1
		   and the NWG kernel estimate is y's sample mean. */
		mean_of_y = mean_of_array_of_doubles(b->y, b->n);
		for (i = 0; i < b->n; i++) {
			residual = mean_of_y - b->y[i];
			ss += residual * residual;
		}
	} else {
		for (i = 0; i < b->n; i++) {
			residual = nwgks_at_leaving_out(b, b->x[i], lambda, i) - b->y[i];
			ss += residual * residual;
		}
	}
	if (isnan(ss)) {
		return INFINITY;
	}

	return sqrt(ss / b->n);
}

/* Comparison function for use with `qsort`.
   Returns -1 if the double at `first` is less than that at `second`;
   returns +1 if the double at `first` is greater;
   returns 0 if the doubles are equal. */
int compare_two_doubles(const void* first, const void* second)
{
	if (*((double*) first) < *((double*) second)) {
		return -1;
	} else if (*((double*) first) > *((double*) second)) {
		return 1;
	} else {
		return 0;
	}
}

/* Assuming the optimal lambda is finite, compute lower & upper bounds
   on the plausible optimal lambda for the data `b`, storing those
   bounds in `bound`. */
bool bound_nwgks_lambda(const Bivariate* b, double* bound)
{
	unsigned int i;
	double cur_abs_diff;
	double min_abs_diff = INFINITY;
	double* sorted_x;
	double x_max;
	double x_min;

	/* Make a sorted list of the data's abscissae (x-values). */
	if ((sorted_x = malloc(b->n * sizeof(double))) == NULL) {
		return false;
	}
	memcpy(sorted_x, b->x, b->n * sizeof(double));
	qsort(sorted_x, b->n, sizeof(double), compare_two_doubles);

	x_min = sorted_x[0];
	x_max = sorted_x[b->n - 1];

	/* Find the minimum gap between adjacent values in the sorted list. */
	for (i = 1; i < b->n; i++) {
		cur_abs_diff = sorted_x[i] - sorted_x[i - 1];
		if ((cur_abs_diff > 0.0) && (cur_abs_diff < min_abs_diff)) {
			min_abs_diff = cur_abs_diff;
		}
	}

	bound[0] = min_abs_diff / REL_BOUND_LOOSENESS;  /* lower */
	bound[1] = 0.5 * REL_BOUND_LOOSENESS * (x_max - x_min);  /* upper */

	free(sorted_x);

	return true;
}

/* Deduce the parabola through three distinct (x, y) points and return
   the x-coordinate (abscissa) of the parabola's turning point. */
double parab_interp(const double *x, const double *y)
{
	double A = (y[1] - y[0]) / (x[1] - x[0]);
	double B;

	B = ((y[2] - y[0]) / (x[2] - x[0])) - ((y[1] - y[0]) / (x[1] - x[0]));
	B /= (x[2] - x[1]);

	return ((B * (x[0] + x[1])) - A) / (2.0 * B);
}

/* Solve for a parabola's turning point as `parab_interp` does, but with
   `x` and `y` taken on logarithmic scales. This seems to be a better form
   of parabolic interpolation for minimizing the LOO RMS. */
double parab_interp_log_scale(double* x, double* y)
{
	double log_x[3] = { log(x[0]), log(x[1]), log(x[2]) };
	double log_y[3] = { log(y[0]), log(y[1]), log(y[2]) };

	return exp(parab_interp(log_x, log_y));
}

/* Fit a Nadaraya-Watson Gaussian-kernel smooth to the bivariate dataset `b`,
   iteratively estimating the parameter `lambda` and recording the final
   fit's leave-one-out RMS error in `loo_rms`, the iterations taken to
   converge in `iterations_for_caller`, and the residuals' estimated
   variance in `est_resid_var`.

   `iterations_for_caller` may be `NULL`. If so, this function does not
   record the number of iterations at that address, but merely counts
   iterations internally.

   If `lambda` is a normal nonzero floating-point number, this function uses
   `lambda` as a starting value; note that in this case the caller must also
   initialize `loo_rms` to the corresponding LOO RMS (or an overestimate
   thereof; an overestimate is OK because it will simply cause this function
   to experiment with other lambda values). If `lambda` is zero, infinite, or
   subnormal, this function picks its own starting value. */

bool nwgks_fit_by_loo(const Bivariate* b, double* lambda, double* loo_rms, unsigned int* iterations_for_caller, double* est_resid_var)
{
	double b_rms[2];
	double bound[2];
	double caller_lambda = 0.0;
	double caller_loo_rms = INFINITY;
	unsigned int iterations = 0;
	double lb_bound[2];  /* bounds on `bound[0]` */
	double multiplier;
	double prev_rms;
	double prev_lambda;
	char step_type = '\0';  /* 'p' for para. interp., 'o' otherwise */

	/* Variables for parabolic interpolation: an abscissa difference,
	   lambda-RMS pairs, a flag used to decide how & when to generate
	   the next trial lambda value by parabolic interpolation. */
	double absc_diff;
	double nearby_lambda[3] = { -INFINITY, 0.0, INFINITY };
	double nearby_rms[3] = { INFINITY, INFINITY, INFINITY };
	bool parabolic_interpolation_worthwhile;

	/* These are used only for computing the residual variance (which is
	   done only if `est_resid_var` is non-NULL. */
	unsigned int i;
	double resid;
	double mean_squ_resid;

	if (!bound_nwgks_lambda(b, bound)) {
		/* The attempt to bound lambda for this dataset failed. */
		if (iterations_for_caller != NULL) {
			*iterations_for_caller = 0;
		}
		return false;
	}

	if (fpclassify(*lambda) == FP_NORMAL) {
		caller_lambda = *lambda;
		caller_loo_rms = *loo_rms;
	}

	b_rms[0] = nwgks_loo_rms(b, bound[0]);
	b_rms[1] = nwgks_loo_rms(b, bound[1]);
	#ifdef NWGKS_DEBUG_OUTPUT_TO_STDERR
	fprintf(stderr, "NFBL: bound %g has LOO RMS %g\n", bound[0], b_rms[0]);
	fprintf(stderr, "NFBL: bound %.15g has LOO RMS %.15g\n",
	        bound[1], b_rms[1]);
	#endif

	prev_lambda = bound[0];
	prev_rms = b_rms[0];
	*lambda = bound[1];
	*loo_rms = b_rms[1];

	/* `bound[0]`, the quick lower bound on lambda, is typically so low that
	   its LOO RMS error is infinity for numerical reasons. If so, that lower
	   bound is too low; iteratively raise it until its RMS is finite.

	   This procedure does risk cutting out a small interval of permissible
	   lambda values just below the revised lower bound, but using a fine
	   tolerance forces that interval to be tiny, and for every test input
	   I've never seen signs of the optimal lambda being in that interval.
	   And even in the very unlikely worst-case scenario where the optimal
	   lambda is in the excluded interval, the result is just that the
	   final lambda will be marginally suboptimal. */

	if (isinf(b_rms[0])) {
		lb_bound[0] = bound[0];
		lb_bound[1] = bound[1];
		while (((lb_bound[1] - lb_bound[0]) > (NWGKS_REL_TOL_BO * (bound[1] - bound[0]))) || isinf(*loo_rms)) {
			iterations++;
			*lambda = sqrt(lb_bound[0] * lb_bound[1]);
			*loo_rms = nwgks_loo_rms(b, *lambda);
			#ifdef NWGKS_DEBUG_OUTPUT_TO_STDERR
			fprintf(stderr,
			        "NFBL: lb_bound, *lambda, *loo_rms: (%g %g) %g %g\n",
			        lb_bound[0], lb_bound[1], *lambda, *loo_rms);
			#endif
			if (isinf(*loo_rms)) {
				lb_bound[0] = *lambda;
			} else {
				lb_bound[1] = *lambda;
			}
			if (*loo_rms < prev_rms) {
				/* Try to further exploit the RMS calculations being done in
				   this loop; if the RMS error for `*lambda` improves on the
				   error of the lambda guess made previously, substitute
				   `*lambda` for that earlier guess. */
				prev_lambda = *lambda;
				prev_rms = *loo_rms;
			}
		}
		#ifdef NWGKS_DEBUG_OUTPUT_TO_STDERR
		fprintf(stderr,
		        "NFBL: lb_bound: (%g %g)\n",
		        lb_bound[0], lb_bound[1]);
		#endif
		if (!isinf(*loo_rms)) {
			bound[0] = *lambda;
			b_rms[0] = *loo_rms;
		} else {
			bound[0] = lb_bound[1];
			b_rms[0] = nwgks_loo_rms(b, lb_bound[1]);
		}
		#ifdef NWGKS_DEBUG_OUTPUT_TO_STDERR
		fprintf(stderr, "NFBL: new lower bound %g has LOO RMS %g\n",
		        bound[0], b_rms[0]);
		#endif
	}

	nearby_lambda[0] = bound[0];
	nearby_rms[0] = b_rms[0];
	nearby_lambda[2] = bound[1];
	nearby_rms[2] = b_rms[1];

	/* `bound` should now have decent lower & upper bounds on the LOO-RMS-
	   minimizing value of lambda. Iteratively converge on that minimizing
	   value of lambda.

	   Firstly, try out a weighted geometric mean of the bounds, weighted
	   towards the lower bound (there's no theory backing this, it just
	   seems a good choice for a few datasets).

	   Then try to iteratively converge on the optimal lambda, step by step.
	   When a step successfully reduces the RMS, increase the step
	   size; if a step raises the RMS instead, shrink the step size and
	   double back. */

	*lambda = pow(bound[0], 0.75) * pow(bound[1], 0.25);
	*loo_rms = nwgks_loo_rms(b, pow(bound[0], 0.75) * pow(bound[1], 0.25));
	#ifdef NWGKS_DEBUG_OUTPUT_TO_STDERR
	fprintf(stderr,
	        "NFBL: 75/25 wei. geo. avg. `lambda` = %g, LOO RMS = %g\n",
	        *lambda, *loo_rms);
	#endif
	if (*loo_rms > prev_rms) {
		/* The weighted geometric average didn't improve on `prev_lambda`,
		   so roll `lambda` (and its accompanying `loo_rms`) back to
		   `prev_lambda` (and its accompanying LOO RMS). */
		*lambda = prev_lambda;
		*loo_rms = prev_rms;
	}

	#ifdef NWGKS_DEBUG_OUTPUT_TO_STDERR
	fprintf(stderr, "NFBL: `prev_lambda` = %g (%g), `lambda` = %g (%g)\n",
	        prev_lambda, prev_rms, *lambda, *loo_rms);
	#endif
	if ((*lambda <= bound[0]) || (*lambda >= bound[1])) {
		/* The preceding statements in this function couldn't come up with
		   a good starting point for lambda. Use the geometric mean of
		   lambda's bounds. */
		*lambda = sqrt(bound[0] * bound[1]);
		*loo_rms = nwgks_loo_rms(b, *lambda);
	}
	if (caller_lambda != 0.0) {
		*lambda = caller_lambda;
		*loo_rms = caller_loo_rms;
	}
	#ifdef NWGKS_DEBUG_OUTPUT_TO_STDERR
	fprintf(stderr, "NFBL: `prev_lambda` = %g (%g), `lambda` = %g (%g)\n",
	        prev_lambda, prev_rms, *lambda, *loo_rms);
	#endif
/*
	step = -REL_INITIAL_STEP * (bound[1] - bound[0]);
*/
	multiplier = 1.0 - REL_INITIAL_STEP;
	do {
		prev_lambda = *lambda;
		prev_rms = *loo_rms;
		iterations++;
		parabolic_interpolation_worthwhile =
			isnormal(nearby_lambda[0])
			&& isnormal(nearby_lambda[2])
			&& ((nearby_lambda[2] / nearby_lambda[0]) < PI_TRIGGER_RATIO);
		if (parabolic_interpolation_worthwhile) {
			*lambda = parab_interp_log_scale(nearby_lambda, nearby_rms);
			step_type = 'p';
			#ifdef NWGKS_DEBUG_OUTPUT_TO_STDERR
			fprintf(stderr, "NFBL: parabolic interpolation gave %.15g\n",
			        *lambda);
			#endif
			/* Once the optimizer's made it into a roughly quadratic
			   bowl of the LOO function near a local minimum, the
			   parabolic interpolator often undershoots the minimum.
			   Try to adjust for that on the assumption that each
			   parabolic-interpolation step only steps halfway to the
			   minimum. */
//			*lambda += (*lambda - nearby_lambda[1]);
			#ifdef NWGKS_DEBUG_OUTPUT_TO_STDERR
			fprintf(stderr, "NFBL: accelerated interpol'd value is %.15g\n",
			        *lambda);
			#endif
			absc_diff = nearby_lambda[2] - nearby_lambda[0];
			parabolic_interpolation_worthwhile =
				(*lambda != prev_lambda)
				&& (*lambda >= (nearby_lambda[0] - (0.5 * absc_diff)))
				&& (*lambda <= (nearby_lambda[2] + (0.5 * absc_diff)));
			if (!parabolic_interpolation_worthwhile) {
				/* The parabolic interpolation, meant to generate a new,
				   better `lambda`, failed, either by producing the same
				   `lambda` as before or by jumping some way outside the
				   region through which it's interpolating. Fall back on
				   trying iterative descent. */
/*
				*lambda = prev_lambda + step;
*/
				*lambda = multiplier * prev_lambda;
				step_type = 'o';
				#ifdef NWGKS_DEBUG_OUTPUT_TO_STDERR
				fprintf(stderr,
				        "NFBL: rejecting parab. interp., using lambda = %g\n",
				        *lambda);
				#endif
			}
		} else {
/*
			*lambda = prev_lambda + step;
*/
			*lambda = multiplier * prev_lambda;
			step_type = 'o';
		}
		while ((*lambda < bound[0]) || (*lambda > bound[1])) {
/*
			step /= 2.0;
			*lambda = prev_lambda + step;
*/
			multiplier = (1.0 + multiplier) / 2.0;
			#ifdef NWGKS_DEBUG_OUTPUT_TO_STDERR
			fprintf(stderr, "NFBL: {1} `multiplier` = %g\n", multiplier);
			#endif
			*lambda = multiplier * prev_lambda;
		}
		*loo_rms = nwgks_loo_rms(b, *lambda);
		#ifdef NWGKS_DEBUG_OUTPUT_TO_STDERR
		fprintf(stderr, "NFBL: `lambda` from %g to %.9g (LR %.15g)\n",
		        prev_lambda, *lambda, *loo_rms);
		#endif
/*
		if (*loo_rms < prev_rms) {
			step *= 2.0;
		} else {
			step /= -2.0;
		}
*/
		if (step_type == 'o') {
			if (*loo_rms < prev_rms) {
//				multiplier *= multiplier;
//				multiplier = pow(multiplier, 2.3);
				multiplier += (multiplier - 1.0);
				#ifdef NWGKS_DEBUG_OUTPUT_TO_STDERR
				fprintf(stderr, "NFBL: {2} `multiplier` = %g\n", multiplier);
				#endif
			} else {
//				multiplier = 1.0 / sqrt(multiplier);
//				multiplier = pow(multiplier, -0.7);
				multiplier = 1.0 + ((1.0 - multiplier) / 2.0);
				#ifdef NWGKS_DEBUG_OUTPUT_TO_STDERR
				fprintf(stderr, "NFBL: {3} `multiplier` = %g\n", multiplier);
				#endif
			}
		}
		if (*loo_rms < nearby_rms[1]) {
			if (*lambda > nearby_lambda[1]) {
				nearby_lambda[0] = nearby_lambda[1];
				nearby_rms[0] = nearby_rms[1];
			} else {
				nearby_lambda[2] = nearby_lambda[1];
				nearby_rms[2] = nearby_rms[1];
			}
			nearby_lambda[1] = *lambda;
			nearby_rms[1] = *loo_rms;
		} else if ((*loo_rms < nearby_rms[2]) && (*lambda > nearby_lambda[1])) {
			nearby_lambda[2] = *lambda;
			nearby_rms[2] = *loo_rms;
		} else if ((*loo_rms < nearby_rms[0]) && (*lambda < nearby_lambda[1])) {
			nearby_lambda[0] = *lambda;
			nearby_rms[0] = *loo_rms;
		}
		#ifdef NWGKS_DEBUG_OUTPUT_TO_STDERR
		fprintf(stderr, "NFBL: parabola points are:\n"
		        "\t1: (%.16g, %.15g)\n"
		        "\t2: (%.16g, %.15g)\n"
		        "\t3: (%.16g, %.15g)\n",
		        nearby_lambda[0], nearby_rms[0],
		        nearby_lambda[1], nearby_rms[1],
		        nearby_lambda[2], nearby_rms[2]);
		#endif
	} while (fabs(*lambda - prev_lambda) > (NWGKS_REL_TOL_BO * (bound[1] - bound[0])));

	if (b_rms[1] <= *loo_rms) {
		/* The upper bound on lambda is at least as good (in terms of LOO
		   RMS error) as the now-converged-upon estimate `*lambda`.
		   That suggests that the true lambda is effectively infinite;
		   inform the caller. */
		*lambda = INFINITY;
		*loo_rms = b_rms[1];
	}

	/* If the caller's requested an estimate of the residuals' variance,
	   compute an estimate by applying Bessel's correction to the mean
	   square of the residuals. */
	if (est_resid_var != NULL) {
		resid = b->y[0] - nwgks_at(b, b->x[0], *lambda);
		mean_squ_resid = resid * resid;
		for (i = 1; i < b->n; i++) {
			resid = b->y[i] - nwgks_at(b, b->x[i], *lambda);
			/* This seems like a way to compute the mean square that
			   is (i) online and (ii) less liable to hit an overflow
			   if many residuals are huge and/or n is arbitrarily
			   large. The second term in the sum will eventually vanish
			   in the infinite-n limit, but infinity's far away. */
			mean_squ_resid = ((i / (double) (1 + i)) * mean_squ_resid)
			                 + ((resid * resid) / (1 + i));
		}
		*est_resid_var = b->n * (mean_squ_resid / (b->n - 1));
	}

	if (iterations_for_caller != NULL) {
		*iterations_for_caller = iterations;
	}

	return true;
}

#ifdef MAKE_STANDALONE_NWGKS

bool invoked_as_nwgksloo(const char* argv_zero)
{
	return (!strcoll(argv_zero, "nwgksloo"))
	       || !strcoll(argv_zero + strlen(argv_zero) - 9, "/nwgksloo");
}

bool probe_loo_rms_function(const Bivariate* b)
{
	unsigned int samples;

	double bound[2];
	unsigned int i;
	double lambda;
	double multiplier;

	if (!bound_nwgks_lambda(b, bound)) {
		/* The attempt to bound lambda for this dataset failed. */
		return false;
	}

	samples = 150;
	if (b->n < 100) {
		samples *= 5;
	}

	/* Compute and write the LOO RMS at logarithmically equally spaced
	   values of lambda over the range of plausible values. */
	lambda = bound[0];
	multiplier = pow(bound[1] / bound[0], 1.0 / (double) samples);
	for (i = 0; i < samples; i++) {
		printf("%g\t%g\n", lambda, nwgks_loo_rms(b, lambda));
		fflush(stdout);  /* make it easier to see intermediate output */
		lambda *= multiplier;
	}

	return true;
}

/* In case you want a standalone kernel-smoothing program...just feed
   lines of data into `stdin`, one whitespace-separated pair of numbers on
   each line.

   Bonus: if invoked as `nwgksloo`, the program computes (and writes to
   standard output) the LOO RMS for the input data, at a sample of possible
   values for the smoothing parameter. */

int main(int argc, char* argv[])
{
	Bivariate b;
	int i;
	double x;
	double y;

	double est_err_var;
	double est_estim_var;
	unsigned int iterations;
	double loo_rms = INFINITY;
	double lambda = 0.0;  /* a value of 0.0 indicates `lambda`'s unknown */

	if (!init_Bivariate(&b)) {
		fputs("Can't initialize bivariate data storage, exiting.\n", stderr);
		return 1;
	}

	do {
		i = scanf("%lg %lg", &x, &y);
		if (i == 2) {
			if (!append_Bivariate(&b, x, y)) {
				fprintf(stderr, "Ran out of memory reading datum %u.\n",
				        1 + b.n);
				break;
			}
		} else if (i != EOF) {
			do {
				i = getchar();
			} while ((i != '\n') && (i != EOF));
		}
	} while ((!feof(stdin)) && !ferror(stdin));

	if (b.n < 3) {
		fprintf(stderr, "Read in only %u row", b.n);
		if (b.n > 1) {
			putc('s', stderr);
		}
		fputs(" of data, but > 2 are necessary.\n", stderr);
		free_Bivariate(&b);
		return 2;
	}

	/* If this program was invoked as `nwgksloo`, compute the approximate
	   shape of the LOO RMS function. The check of `argc` is unnecessary
	   but suppresses a pointless `gcc` warning. */
	if (argc && invoked_as_nwgksloo(argv[0])) {
		/* Just mechanically calculate how the LOO RMS varies
		   varies with lambda, then finish. */
		i = (int) !probe_loo_rms_function(&b);
		free_Bivariate(&b);
		return i;
	}

	nwgks_fit_by_loo(&b, &lambda, &loo_rms, &iterations, &est_err_var);
	printf("# %u iterations: lambda %g, LOO RMS %g, and est. err. var. %g\n",
	       iterations, lambda, loo_rms, est_err_var);

	for (i = 0; ((unsigned int) i) < b.n; i++) {
		printf("%g %.8g %.8g ",
		       b.x[i], b.y[i],
		       nwgks_at_with_var(&b, b.x[i], lambda,
		                         &est_estim_var, est_err_var));
		printf("%g\n", est_estim_var);
	}

	free_Bivariate(&b);

	return 0;
}

#endif
