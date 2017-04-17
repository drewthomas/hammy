/* homoscedastic additive modelling */

#include <stdio.h>
#include <stdlib.h>

#include "data.h"
#include "nwgks.h"

/* Maximum number of columns of input to accept. */
#define MAX_COLS 99

/* When iteratively fitting a homoscedastic additive model, stop the fitting
   when the relative change in leave-one-out RMS error is less than
   `REL_LOO_RMS_CHANGE_TOL` (on the assumption that the fitting has
   essentially converged at that point). */
#define REL_LOO_RMS_CHANGE_TOL 1e-5

/* Maximum number of fitting iterations to carry out after the initial fit. */
#define MAX_REFITS 999

/* Define this macro to have `hammy` write debugging output to standard
   error as it runs. */
#define HAMMY_DEBUG_OUTPUT_TO_STDERR

void fit_for_initial_param_estimates(Bivariate* b, const Multivariate* w, unsigned int max_predictor_col_idx, double* lambda, double* loo_rms, unsigned int* fit_iterations, double* est_err_var, double mean_y)
{
	unsigned int i;
	unsigned int j;

	/* Mean-centre the dependent (y) variable. */
	for (i = 0; i < w->n; i++) {
		w->x[w->k - 1][i] -= mean_y;
	}

	#ifdef HAMMY_DEBUG_OUTPUT_TO_STDERR
	fputs(" fit     lambda it's     LOO RMS   prev LR     err. var\n",
	      stderr);
	#endif

	/* For each predictor variable in turn, run a one-dimensional smooth
	   of the y variable against that predictor variable. Then subtract the
	   smoothed estimates of y from the y values in the original
	   multivariate input `w`. The idea here is to subtract off each
	   predictor's contribution in turn. */

	for (i = 0; i < max_predictor_col_idx; i++) {

		extract_Bivariate_from_Multivariate(&(b[i]), w, i,
		                                    max_predictor_col_idx);

		lambda[i] = 0.0;  /* indicate that `lambda[i]` is unknown */
	    nwgks_fit_by_loo(&(b[i]), &(lambda[i]), &(loo_rms[i]),
		                 fit_iterations, &(est_err_var[i]));

		#ifdef HAMMY_DEBUG_OUTPUT_TO_STDERR
		fprintf(stderr, " 1.%u: %9g %4u %11g           %12g\n",
		        1 + i, lambda[i], *fit_iterations, loo_rms[i],
		        est_err_var[i]);
	    #endif

		for (j = 0; j < w->n; j++) {
			w->x[max_predictor_col_idx][j] -= nwgks_at(&(b[i]), b[i].x[j],
			                                           lambda[i]);
		}

	}
}

bool fit(const Multivariate* d, Multivariate* w, double* lambda)
{
	Bivariate b[d->k - 1];
	unsigned int col_idx;
	unsigned int i;
	unsigned int j;
	double mean_y = d->x[d->k - 1][0];
	double y;

	double best_loo_rms = 1e300;
	double est_estim_var[MAX_COLS - 1];
	double est_err_var[MAX_COLS - 1];
	unsigned int fit_iterations;
	double loo_rms[MAX_COLS - 1];
	double prev_loo_rms = best_loo_rms;
	unsigned int re_fits = 0;
	double total_est_estim_var;

	if (!copy_and_init_Multivariate(d, w)) {
		fputs("Can't make working copy of input data.\n", stderr);
		return false;
	}

	for (i = 0; i < (d->k - 1); i++) {
		if (!init_Bivariate(&(b[i]))) {
			fputs("Can't initialize temporary bivariate data storage.\n",
			      stderr);
			return false;
		}
	}

	for (i = 1; i < d->n; i++) {
		mean_y = ((i / (double) (1 + i)) * mean_y)
		         + (d->x[d->k - 1][i] / (1 + i));
	}

	/* Estimate the smooth function (lambda's value, effectively) for each
	   predictor in turn, as a starting point to improve. */
	fit_for_initial_param_estimates(b, w, d->k - 1, lambda, loo_rms,
	                                &fit_iterations, est_err_var, mean_y);

	/* Repeatedly fit one-dimensional smooths for each predictor
	   variable, using the results on one iteration as the starting
	   point of the next. */

	do {

		prev_loo_rms = loo_rms[d->k - 2];
		if (prev_loo_rms < best_loo_rms) {
			best_loo_rms = prev_loo_rms;
		}
		re_fits++;

		for (i = 0; i < (d->k - 1); i++) {

			extract_Bivariate_from_Multivariate(&(b[i]), d, i, d->k - 1);
			for (j = 0; j < b[i].n; j++) {
				b[i].y[j] -= mean_y;
			}
			for (col_idx = 0; col_idx < (d->k - 1); col_idx++) {
				if (col_idx == i) {
					continue;
				}
				for (j = 0; j < b[i].n; j++) {
					b[i].y[j] -= nwgks_at(&(b[col_idx]), b[col_idx].x[j],
					                      lambda[col_idx]);
				}
			}
			nwgks_fit_by_loo(&(b[i]), &(lambda[i]), &(loo_rms[i]),
			                 &fit_iterations, &(est_err_var[i]));
			#ifdef HAMMY_DEBUG_OUTPUT_TO_STDERR
		    fprintf(stderr,
		            "%2u.%u: %9g %4u %11g %9.4g %12g\n",
		            1 + re_fits, 1 + i, lambda[i],
			        fit_iterations, loo_rms[i], best_loo_rms, est_err_var[i]);
			#endif
		}

	} while ((re_fits < MAX_REFITS) && ((loo_rms[d->k - 2] / prev_loo_rms) < (1.0 - REL_LOO_RMS_CHANGE_TOL)));

	/* Contrary to hopes & expectations, the fitting iterations may not
	   have improved the fit in the end! Allow for that possibility. */
	if ((re_fits == 1) && (loo_rms[d->k - 2] > best_loo_rms)) {
		copy_Multivariate(d, w);
		fit_for_initial_param_estimates(b, w, d->k - 1, lambda, loo_rms,
		                                &fit_iterations, est_err_var, mean_y);
		fputs("Reverting to very first fit's param. estimates\n", stderr);
	}

	/* Write the results for each datum to standard output: the value of
	   each predictor, the final estimate of y from the combined,
	   multivariate smooth, the variance of that smoothed estimate of y (NB:
	   that is not the same as the variance in the observed y, which includes
	   the error variance!), and the estimate & variance of each predictor's
	   additive contribution to the smoothed y value. */
	for (i = 0; i < d->n; i++) {
		for (j = 0; j < d->k; j++) {
			printf("%g ", d->x[j][i]);
		}
		y = mean_y;
		total_est_estim_var = 0.0;
		for (j = 0; j < d->k - 1; j++) {
			y += nwgks_at_with_var(&(b[j]),
			                       b[j].x[i],
			                       lambda[j],
			                       &(est_estim_var[j]),
			                       est_err_var[j]);
			total_est_estim_var += est_estim_var[j];
		}
		printf("%g %g ", y, total_est_estim_var);
		for (j = 0; j < (d->k - 1); j++) {
			printf("%g ", nwgks_at_with_var(&(b[j]), b[j].x[i], lambda[j],
			                                &(est_estim_var[j]),
			                                est_err_var[j]));
			printf("%g", est_estim_var[j]);
			if (j == (d->k - 2)) {
				putchar('\n');
			} else {
				putchar(' ');
			}
		}
	}

	for (i = 0; i < (d->k - 1); i++) {
		free_Bivariate(&(b[i]));
	}

	return true;
}

int main(int argc, char* argv[])
{
	unsigned int cols_so_far;
	Multivariate d;
	int i;
	unsigned int input_cols;
	double row[MAX_COLS];
	Multivariate working;

	double lambda[MAX_COLS - 1];

	if (argc < 2) {
		fprintf(stderr, "Usage: %s NUM_Y_COLS\n", argv[0]);
		return 1;
	}

	input_cols = strtoul(argv[1], NULL, 10);
	if (input_cols < 2) {
		fputs("There must be at least 2 columns of input.\n", stderr);
		return 2;
	} else if (input_cols > MAX_COLS) {
		fprintf(stderr,
		        "Sorry, %u columns of input is too many; "
		        "the maximum is %u.\n",
		        input_cols, MAX_COLS);
		return 3;
	}

	if (!init_Multivariate(&d, input_cols)) {
		fputs("Can't initialize multivariate data storage, exiting.\n",
		      stderr);
		return 4;
	}

	do {
		i = 0;  /* stop the compiler whinging about `i` maybe being undef'd */
		for (cols_so_far = 0; cols_so_far < input_cols; cols_so_far++) {
			i = scanf("%lg", &(row[cols_so_far]));
			if (i != 1) {
				break;
			}
		}
		if (cols_so_far == input_cols) {
			/* A complete row of data was just read, try to append it. */
			if (!append_Multivariate(&d, row)) {
				fprintf(stderr, "Ran out of memory reading datum %u.\n",
				        1 + d.n);
				break;
			}
		} else if (i != EOF) {
			do {
				i = getchar();
			} while ((i != '\n') && (i != EOF));
		}
	} while ((!feof(stdin)) && !ferror(stdin));

	if (d.n < 4) {
		fprintf(stderr, "Read in only %u row", d.n);
		if (d.n > 1) {
			putc('s', stderr);
		}
		fputs(" of data, but > 3 are necessary.\n", stderr);
		return 5;
	}

	if (!fit(&d, &working, lambda)) {
		fputs("Failed to allocate memory when attempting fit.\n", stderr);
		free_Multivariate(&d);
		return 6;
	}

	fputs("lambda (", stderr);
	for (i = 0; i < (int) (d.k - 2); i++) {
		fprintf(stderr, "%g, ", lambda[i]);
	}
	fprintf(stderr, "%g)\n", lambda[d.k - 2]);

	free_Multivariate(&working);
	free_Multivariate(&d);

	return 0;
}
