// This file is distributed under the GNU AGPL.  See file LICENSE for details.
// Copyright (C) 2013 Enric Meinhardt-Llopis <http://dev.ipol.im/~coco>

#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "iio.h"


// auxiliary function: exit the program with an error message
static void fail(const char *fmt, ...)

{
	va_list argp;
	fprintf(stderr, "\nFAIL: ");
	va_start(argp, fmt);
	vfprintf(stderr, fmt, argp);
	va_end(argp);
	fprintf(stderr, "\n\n");
	fflush(NULL);
	exit(-1);
}

// auxiliary function (like malloc, but always returns a valid pointer)
static void *xmalloc(size_t size)
{
	void *p = malloc(size);
	if (!p)
		exit(fprintf(stderr,
			"ERROR: out of mem requesting %zu butes\n", size));
	return p;
}

// y[k] = (1/n) * sum_i x[i][k]
static void float_avg(float *y, float *xx, int d, int n)
{
	float (*x)[d] = (void*)xx;
	for (int k = 0; k < d; k++)
	{
		y[k] = 0;
		for (int i = 0; i < n; i++)
			y[k] += x[i][k]/n;
	}
}

// euclidean distance between the vectors x and y, regularized around 0
static float fdiste(float *x, float *y, int n, float e)
{
	return n ? hypot(*x - *y, fdiste(x + 1, y + 1, n - 1, e)) : e;
}

// Parameters for the Weiszfeld algorithm
//
// These values are fairly conservative for computing the median
// of up to a few thousand vectors on a cube of side [0,255]
#define WEISZ_NITER 6
#define WEISZ_EPSILON 1e-5

// y[k] = euclidean median of the vectors x[i][k]
static void float_weisz(float *y, float *x, int d, int n)
{
	float_avg(y, x, d, n);
	int niter = WEISZ_NITER;
	for (int k = 0; k < niter; k++) {
		float a[d], b = 0;
		for (int l = 0; l < d; l++)
			a[l] = 0;
		for (int i = 0; i < n; i++) {
			float dxy = fdiste(x + i*d, y, d, WEISZ_EPSILON);
			for (int l = 0; l < d; l++)
				a[l] += x[i*d + l]/dxy;
			b += 1/dxy;
		}
		for (int l = 0; l < d; l++)
			y[l] = a[l]/b;
	}
}

// whether all the components of a vector are finite (e.g. not NAN or INF)
static bool isgood(float *x, int n)
{
	for (int i = 0; i < n; i++)
		if (!isfinite(x[i]))
			return false;
	return true;
}

// auxiliary function to parse command line arguments
//
// @c pointer to original argc
// @v pointer to original argv
// @o option name (after hyphen)
// @d default value
static char *pick_option(int *c, char ***v, char *o, char *d)
{
	int argc = *c;
	char **argv = *v;
	int id = d ? 1 : 0;
	for (int i = 0; i < argc - id; i++)
		if (argv[i][0] == '-' && 0 == strcmp(argv[i]+1, o))
		{
			char *r = argv[i+id]+1-id;
			*c -= id+1;
			for (int j = i; j < argc - id; j++)
				(*v)[j] = (*v)[j+id+1];
			return r;
		}
	return d;
}

int main(int c, char *v[])
{
	// parse command line arguments
	char *filename_out = pick_option(&c, &v, "o", "-");
	if (c < 4) {
		fprintf(stderr,
		"usage:\n\t%s {sum|min|max|avg|weisz} [v1 ...] [-o out]\n", *v);
		//          0  1                          2  3
		return EXIT_FAILURE;
	}
	int n = c - 2;
	char *operation_name = v[1];
	void (*f)(float*,float*,int,int) = NULL;
	if (0 == strcmp(operation_name, "avg"))   f = float_avg;
	if (0 == strcmp(operation_name, "weisz"))   f = float_weisz;
	if (!f) fail("unrecognized operation \"%s\"", operation_name);

	// read input images
	float *x[n];
	int w[n], h[n], pd[n];
	for (int i = 0; i < n; i++)
		x[i] = iio_read_image_float_vec(v[i+2], w + i, h + i, pd + i);
	for (int i = 0; i < n; i++) {
		if (w[i] != *w || h[i] != *h || pd[i] != *pd)
			fail("%dth image sizes mismatch\n", i);
	}

	// alloc space for output image
	float (*y) = xmalloc(*w * *h * *pd * sizeof*y);

	// perform the computation
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int i = 0; i < *w * *h; i++) {
		float tmp[n][*pd];
		int ngood = 0;
		for (int j = 0; j < n; j++)
			if (isgood(x[j]+i**pd, *pd)) {
				for (int k = 0; k < *pd; k++)
					tmp[ngood][k] = x[j][i**pd+k];
				ngood += 1;
			}
		f(y + i**pd, tmp[0], *pd, ngood);
	}

	// save result, cleanup and exit
	iio_save_image_float_vec(filename_out, y, *w, *h, *pd);
	free(y);
	for (int i = 0; i < n; i++)
		free(x[i]);
	return EXIT_SUCCESS;
}
