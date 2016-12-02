/* Define this macro to have `nwgks_fit_by_loo` write debugging output to
   `stderr` as it runs. */
/*
#define NWGKS_DEBUG_OUTPUT_TO_STDERR
*/

double nwgks_at(const Bivariate* b, double at, double lambda);
double nwgks_at_with_var(const Bivariate* b, double at, double lambda, double* estimator_variance, double error_variance);
double nwgks_at_leaving_out(const Bivariate* b, double at, double lambda, unsigned int leave_out_idx);

double nwgks_loo_rms(const Bivariate* b, double lambda);

int compare_two_doubles(const void* first, const void* second);

bool bound_nwgks_lambda(const Bivariate* b, double* bound);

bool nwgks_fit_by_loo(const Bivariate* b, double* lambda, double* loo_rms, unsigned int* iterations_for_caller, double* est_resid_var);
