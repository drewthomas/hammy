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

/* Using the data `b`, compute a Nadaraya-Watson Gaussian-kernel estimate
   of the ordinate at the abscissa `at` using the smoothing parameter
   `lambda`, higher values of which produce more aggressive smoothing. */
double nwgks_at(const Bivariate* b, double at, double lambda)
{
	unsigned int i;
	double sum_w = 0.0;
	double w;
	double y = 0.0;

	for (i = 0; i < b->n; i++) {
		w = exp(-(b->x[i] - at) * (b->x[i] - at) / (2.0 * lambda * lambda));
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
	double residual;
	double ss = 0.0;

	for (i = 0; i < b->n; i++) {
		residual = nwgks_at_leaving_out(b, b->x[i], lambda, i) - b->y[i];
		ss += residual * residual;
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

/* Fit a Nadaraya-Watson Gaussian-kernel smooth to the bivariate dataset `b`,
   iteratively estimating the parameter `lambda` and recording the final fit's
   leave-one-out RMS error in `loo_rms`, the iterations taken to converge
   in `iterations_for_caller`, and the residuals' estimated variance in
   `est_resid_var`.

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
	unsigned int iterations = 0;
	double lb_bound[2];  /* bounds on `bound[0]` */
	double prev_rms;
	double prev_lambda;
	double step;

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
	}

	b_rms[0] = nwgks_loo_rms(b, bound[0]);
	b_rms[1] = nwgks_loo_rms(b, bound[1]);
	#ifdef NWGKS_DEBUG_OUTPUT_TO_STDERR
	fprintf(stderr, "NFBL: bound %g has LOO RMS %g\n", bound[0], b_rms[0]);
	fprintf(stderr, "NFBL: bound %.15g has LOO RMS %.16g\n",
	        bound[1], b_rms[1]);
	#endif
	prev_lambda = bound[1];
	prev_rms = b_rms[1];
	if ((caller_lambda != 0.0) && (*loo_rms < prev_rms)) {
		prev_lambda = caller_lambda;
		prev_rms = *loo_rms;
	}

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
		if (caller_lambda == 0.0) {
			/* The caller didn't supply an initial estimate of `lambda`.
			   Use the geometric mean of the lower bound's bounds. */
			*lambda = sqrt(lb_bound[0] * lb_bound[1]);
		}
		while ((lb_bound[1] - lb_bound[0]) > (NWGKS_REL_TOL_BO * (bound[1] - bound[0]))) {
			iterations++;
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
				/* Try to exploit the RMS calculations being done in this
				   loop further; if the RMS error for `*lambda` improves on
				   the error of the lambda guess made previously, substitute
				   `*lambda` for that earlier guess. */
				prev_lambda = *lambda;
				prev_rms = *loo_rms;
			}
			*lambda = sqrt(lb_bound[0] * lb_bound[1]);
		}
		#ifdef NWGKS_DEBUG_OUTPUT_TO_STDERR
		fprintf(stderr,
		        "NFBL: lb_bound, *lambda, *loo_rms: (%g %g) %g %g\n",
		        lb_bound[0], lb_bound[1], *lambda, *loo_rms);
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

	/* `bound` should now have usable lower & upper bounds on the LOO-RMS-
	   minimizing value of lambda. Iteratively converge on that value of
	   lambda. When a step successfully reduces the RMS, double the step
	   size; if a step raises the RMS instead, halve the step size and
	   double back. */
	if (prev_lambda == bound[1]) {
		/* The preceding statements in this function couldn't come up with
		   a better starting point for lambda than lambda's upper bound. Try
		   a more realistic one. */
		prev_lambda = sqrt(bound[0] * bound[1]);
		prev_rms = nwgks_loo_rms(b, prev_lambda);
	}
	#ifdef NWGKS_DEBUG_OUTPUT_TO_STDERR
	fprintf(stderr, "NFBL: `prev_lambda` = %g, `lambda` = %g\n",
	        prev_lambda, *lambda);
	#endif
	step = -REL_INITIAL_STEP * (bound[1] - bound[0]);
	do {
		iterations++;
		*lambda = prev_lambda + step;
		while ((*lambda < bound[0]) || (*lambda > bound[1])) {
			step /= 2.0;
			*lambda = prev_lambda + step;
		}
		*loo_rms = nwgks_loo_rms(b, *lambda);
		#ifdef NWGKS_DEBUG_OUTPUT_TO_STDERR
		fprintf(stderr, "NFBL: from %g (%g) to %.9g (%.16g)\n",
		        prev_lambda, prev_rms, *lambda, *loo_rms);
		#endif
		if (*loo_rms < prev_rms) {
			prev_rms = *loo_rms;
			prev_lambda = *lambda;
			step *= 2.0;
		} else {
			step /= -2.0;
		}
	} while (fabs(step) > (NWGKS_REL_TOL_BO * (bound[1] - bound[0])));

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
			   if many residuals are huge and/or n is arbitrarily.
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

/* In case you want a standalone kernel-smoothing program...just feed
   lines of data into `stdin`, one whitespace-separated pair of numbers on
   each line. */
int main(void)
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
