/* Define this macro to have `nwgks_fit_by_loo` write debugging output to
   `stderr` as it runs. */
/*
#define NWGKS_DEBUG_OUTPUT_TO_STDERR
*/

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

double nwgks_at(const Bivariate* b, double at, double lambda);
double nwgks_at_with_var(const Bivariate* b, double at, double lambda, double* estimator_variance, double error_variance);
double nwgks_at_leaving_out(const Bivariate* b, double at, double lambda, unsigned int leave_out_idx);

double nwgks_loo_rms(const Bivariate* b, double lambda);

int compare_two_doubles(const void* first, const void* second);

bool bound_nwgks_lambda(const Bivariate* b, double* bound);

bool nwgks_fit_by_loo(const Bivariate* b, double* lambda, double* loo_rms, unsigned int* iterations_for_caller, double* est_resid_var);
