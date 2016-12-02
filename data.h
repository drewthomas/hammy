#include <stdbool.h>

typedef struct Bivariate {
	unsigned int a;  /* number of data points for which space is alloc'd */
	unsigned int n;  /* actual number of data points */
	double* x;
	double* y;
} Bivariate;

typedef struct Multivariate {
	unsigned int a;  /* number of data points for which space is alloc'd */
	unsigned int n;  /* actual number of data points */
	unsigned int k;  /* number of columns */
	double** x;
} Multivariate;

bool init_Bivariate(Bivariate* d);
bool expand_Bivariate_memory(Bivariate* d, unsigned int new_a);
bool append_Bivariate(Bivariate* d, double new_x, double new_y);
void free_Bivariate(Bivariate* d);

bool init_Multivariate(Multivariate* d, unsigned int k);
bool append_Multivariate(Multivariate* d, const double* new_x);
void copy_Multivariate(const Multivariate* d, Multivariate* new_d);
bool copy_and_init_Multivariate(const Multivariate* d, Multivariate* new_d);
void print_sample_Multivariate(FILE* fp, const Multivariate* d);
bool extract_Bivariate_from_Multivariate(Bivariate* b, const Multivariate* m, unsigned int x_col_idx, unsigned y_col_idx);
void free_Multivariate(Multivariate* d);
