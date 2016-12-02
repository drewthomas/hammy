#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "data.h"

/* When initializing a `Bivariate` or `Multivariate` data matrix, allocate
   memory for this many rows. */
#define INITIAL_ROW_ALLOCATION 1000

bool init_Bivariate(Bivariate* d)
{
	d->a = INITIAL_ROW_ALLOCATION;
	d->n = 0;
	if ((d->x = malloc(2 * d->a * sizeof(double))) == NULL) {
		return false;
	}
	d->y = d->x + d->a;
	return true;
}

bool expand_Bivariate_memory(Bivariate* d, unsigned int new_a)
{
	double* old_x_ptr = d->x;

	/* Check that it's actually necessary to expand the memory `d`'s
	   already allocated. */
	if (new_a <= d->a) {
		return true;
	}

	/* Record the new number of rows for which memory's needed, and
	   reallocate memory accordingly. */
	d->a = new_a;
	if ((d->x = realloc(d->x, 2 * d->a * sizeof(double))) == NULL) {
		/* Can't allocate more memory for `d`. Roll it back to how it
		   was before this function was called. */
		d->a = d->n;
		d->x = old_x_ptr;
		return false;
	}

	/* `realloc` is mischievous; it might've moved the region of memory
	   that's being used! Allow for that possibility. */
	if (d->x != old_x_ptr) {
		memmove(d->x, old_x_ptr, d->n * sizeof(double));
	}

	d->y = d->x + d->a;
	memmove(d->y, old_x_ptr + d->n, d->n * sizeof(double));

	return true;
}

bool append_Bivariate(Bivariate* d, double new_x, double new_y)
{
	/* Is `d` full? Allocate more memory for it if so. */
	if (d->n == d->a) {
		if (!expand_Bivariate_memory(d, 2 * (1 + d->a))) {
			return false;  /* failed to allocate more memory */
		}
	}

	/* Put the newest (x, y) pair in the next available slot in `d`. */
	d->x[d->n] = new_x;
	d->y[d->n] = new_y;
	d->n++;

	return true;
}

void free_Bivariate(Bivariate* d)
{
	free(d->x);
	d->a = 0;
	d->n = 0;
	d->x = NULL;
	d->y = NULL;
}

bool init_Multivariate(Multivariate* d, unsigned int k)
{
	unsigned int col;

	d->a = INITIAL_ROW_ALLOCATION;
	d->n = 0;
	d->k = k;

	if ((d->x = malloc(d->k * sizeof(double*))) == NULL) {
		return false;
	}
	if ((d->x[0] = malloc(d->k * d->a * sizeof(double))) == NULL) {
		free(d->x);
		return false;
	}

	for (col = 1; col < k; col++) {
		d->x[col] = d->x[0] + (col * d->a);
	}

	return true;
}

bool append_Multivariate(Multivariate* d, const double* new_x)
{
	unsigned int col;
	double* old_x_ptr = d->x[0];

	/* Is `d` full? Allocate more memory for it if so. */
	if (d->n == d->a) {

		/* Try (approximately) doubling the number of slots for data rows. */
		d->a = (2 * (1 + d->a));
		if ((d->x[0] = realloc(d->x[0],
		                       d->k * d->a * sizeof(double))) == NULL) {
			/* Can't allocate more memory for `d`. Roll it back to how it
			   was before this function was called. */
			d->a = d->n;
			d->x[0] = old_x_ptr;
			return false;
		}

		/* `realloc` is mischievous; it might've moved the region of memory
		   that's being used! Allow for that possibility. */
		if (d->x[0] != old_x_ptr) {
			memmove(d->x[0], old_x_ptr, d->k * d->n * sizeof(double));
		}

		/* Excepting the first column, the pointer to each column's data
		   has to be adjusted to point to the new memory slot, and then the
		   actual column data themselves have to be moved. Do this in
		   reverse order of column index, otherwise an earlier column's
		   data move can clobber the data for a later column before the
		   later column's data move occurs.*/
		for (col = d->k - 1; col > 0; col--) {
			d->x[col] = d->x[0] + (col * d->a);
			memmove(d->x[col],
			        old_x_ptr + (col * d->n),
			        d->n * sizeof(double));
		}
	}

	/* Copy the new row of data into the next available empty row in `d`. */
	for (col = 0; col < d->k; col++) {
		d->x[col][d->n] = new_x[col];
	}
	d->n++;

	return true;
}

void copy_Multivariate(const Multivariate* d, Multivariate* new_d)
{
	unsigned int col;

	if (new_d->a != d->a) {
		/* Caller wrongly called this function even though `new_d` doesn't
		   have the right amount of memory allocated to match `d`. Try to
		   proceed anyway -- not this function's problem, guv! */
	}

	new_d->n = d->n;
	new_d->k = d->k;

	memmove(new_d->x[0], d->x[0], d->k * d->a * sizeof(double));

	for (col = 1; col < d->k; col++) {
		new_d->x[col] = new_d->x[0] + (d->x[col] - d->x[0]);
	}
}

bool copy_and_init_Multivariate(const Multivariate* d, Multivariate* new_d)
{
	new_d->a = d->a;
	new_d->n = d->n;
	new_d->k = d->k;

	if ((new_d->x = malloc(d->k * sizeof(double*))) == NULL) {
		return false;
	}
	if ((new_d->x[0] = malloc(d->k * d->a * sizeof(double))) == NULL) {
		free(new_d->x);
		return false;
	}

	copy_Multivariate(d, new_d);

	return true;
}

void print_sample_Multivariate(FILE* fp, const Multivariate* d)
{
	unsigned int col;
	unsigned int row_idx;

	for (row_idx = 0; row_idx < d->n; row_idx += (1 + (d->n / 10))) {
		fprintf(fp, "%4u: ", 1 + row_idx);
		for (col = 0; col < (d->k - 1); col++) {
			fprintf(fp, "% 9g ", d->x[col][row_idx]);
		}
		fprintf(fp, "% 9g\n", d->x[d->k - 1][row_idx]);
	}
}

bool extract_Bivariate_from_Multivariate(Bivariate* b, const Multivariate* m, unsigned int x_col_idx, unsigned y_col_idx)
{
	if (!expand_Bivariate_memory(b, m->n)) {
		return false;
	}

	b->n = m->n;

	memcpy(b->x, m->x[x_col_idx], m->n * sizeof(double));
	memcpy(b->y, m->x[y_col_idx], m->n * sizeof(double));

	return true;
}

void free_Multivariate(Multivariate* d)
{
	free(d->x[0]);
	free(d->x);
	d->x = NULL;

	d->a = 0;
	d->n = 0;
	d->k = 0;
}
