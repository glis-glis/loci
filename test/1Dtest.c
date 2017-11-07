/* Copyright (C) 2017 Andreas Füglistaler <andreas.fueglistaler@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "../include/loci.h"
#include "testhelpers.h"
#include "runtests.h"

/* Helpers */

/* Map function `f` on `xs` with length `len` and write to `out`. */
static void mapf(double (*f)(double x), int len,
            const double xs[], double out[])    {
    for(int i = 0; i < len; i++) {
        out[i]  = (*f)(xs[i]);
    }
}

/* Convert range `rg` to array `out`. */
static void r2a(const Range *rg, double out[]) {
    for(int i = 0; i < rg->len; i++) {
        out[i]  = rg->x0 + rg->dx*i;
    }
}

/* Test functions and derivative. */
static double f1(double x)      {return 1;}
static double f1_x(double x)    {return 0;}

static double f2(double x)      {return x;}
static double f2_x(double x)    {return 1;}

static double f3(double x)      {return 2*x + 1;}
static double f3_x(double x)    {return 2;}

static double f4(double x)      {return log(x*x + 1);}
static double f4_x(double x)    {return 2*x/(x*x + 1);}

/* Test function arrays. */
static double (*funs[4])(double x)    = {f1, f2, f3, f4};
static double (*funs_x[4])(double x)  = {f1_x, f2_x, f3_x, f4_x};

/* Test basic functionality:
 * create interpolation
 * interpolate a point
 * interpolate derivative
 * map interpolation
 * map derivation. */
void basic1D()  {
    OUT(0, "1D Basic Test");

    /* Interpolation range */
    Range rg    = {.x0 = 1, .dx = 0.3, .len = 20};
    /* Test range */
    Range rgt   = {rg.x0 + rg.dx/3, rg.dx, rg.len-1};

    /* Test value */
    double xt   = rg.x0 + rg.dx*rg.len/2 + rg.dx/7;

    int l       = rg.len;   
    int lt      = rgt.len;  

    /* allocate memorry for input and outpout values */
    double *xs0     = malloc(l*sizeof(double));
    double *fs      = malloc(l*sizeof(double));
    double *fs_x    = malloc(l*sizeof(double));

    double *xst     = malloc(lt*sizeof(double));
    double *fst     = malloc(lt*sizeof(double));
    double *fst_x   = malloc(lt*sizeof(double));

    double *res     = malloc(lt*sizeof(double));

    /* Transorm range to array */
    r2a(&rg, xs0);
    r2a(&rgt, xst);

    for(int i = 0; i < 4; i++)  {
        OUT(1, "Function %i", i+1);

        /* Exact values */
        double ft      = funs[i](xt);
        double ft_x    = funs_x[i](xt);

        /* Map exact values on grid points */
        mapf(funs[i], l, xs0, fs);
        mapf(funs_x[i], l, xs0, fs_x);

        /* Create interpolation */
        Interpolation ip    = interpolation(&rg, fs, fs_x);

        /* Compare interpolated value to exact value */
        double err  = fabs(interpolate(&ip, xt) - ft);
        assert((err < tols[1] + tols[1]*fabs(ft)) &&
               "interpolation(xt)"); 

        /* Compare interpolated derivative to exact value */
        err  = fabs(diff_x(&ip, xt) - ft_x);
        assert((err < tols[2] + tols[2]*fabs(ft)) &&
               "diff_x(xt)"); 

        /**/

        /* Norm of arrays */
        double nfs      = norm(lt, fs);
        double nfs_x    = norm(lt, fs_x);

        /* Compare interpolated values to exact values on grid-points */
        map(&ip, lt, xs0, res);
        err  = errnorm(lt, fs, res);
        assert((err < tols[0] + tols[0]*nfs) &&
               "map(interpolation, xs0)"); 

        /* Compare interpolated derivatives to exact values on grid-points */
        map_x(&ip, lt, xs0, res);
        err  = errnorm(lt, fs_x, res);
        assert((err < tols[0] + tols[0]*nfs_x) &&
               "map(diff_x, xs0)"); 

        /**/

        /* Map exact values on test points */
        mapf(funs[i], lt, xst, fst);
        mapf(funs_x[i], lt, xst, fst_x);
        double nfst     = norm(lt, fst);
        double nfst_x   = norm(lt, fst_x);

        /* Compare interpolated values to exact values on test-points */
        map(&ip, lt, xst, res);
        err  = errnorm(lt, fst, res);
        assert((err < tols[1] + tols[1]*nfst) &&
               "map(interpolation, xst)"); 

        /* Compare interpolated derivatives to exact values on test-points */
        map_x(&ip, lt, xst, res);
        err  = errnorm(lt, fs_x, res);
        err  = errnorm(lt, fst_x, res);
        assert((err < tols[2] + tols[2]*nfst_x) &&
               "map(diff_x, xst)"); 

        loci_free(&ip);
    }



    free(xs0);
    free(fs);
    free(fs_x);
    free(xst);
    free(fst);
    free(fst_x);

    free(res);

    OUT(0, "Pass\n");
}

/* Test higher order derivative on exponantial function. */
void exp1D()  {
    OUT(0, "1D Exp Test");

    /* See basic1D() */
    Range rg    = {.x0 = 1, .dx = 0.3, .len = 20};
    Range rgt   = {rg.x0 + rg.dx/3, rg.dx, rg.len - 1};

    double xt   = rg.x0 + rg.dx*rg.len/2 + rg.dx/7;
    double ft   = exp(xt);

    int l       = rg.len;
    int lt      = rgt.len;

    double *fs      = malloc(l*sizeof(double));
    double *xs0     = malloc(l*sizeof(double));
    double *xst     = malloc(lt*sizeof(double));
    double *fst     = malloc(lt*sizeof(double));
    double *res     = malloc(lt*sizeof(double));

    r2a(&rg, xs0);
    r2a(&rgt, xst);

    mapf(exp, l, xs0, fs);

    mapf(exp, lt, xst, fst);
    double nfs  = norm(lt, fst);

    /* Create interpolation */
    Interpolation ip    = interpolation(&rg, fs, fs);

    for(int i = 0; i < 3; i++)  {
        OUT(1, "d^%if/dx^%i", i, i);

        double ti   = tols[i+1];

        /* Compare interpolated derivative to exact value */
        double err  = fabs(diff(&ip, i, xt) - ft);
        assert((err < ti + ti*fabs(ft)) &&
               "diff(i, xt)"); 

        /* Compare interpolated derivatives to exact values on test-points */
        map_diff(&ip, i, lt, xst, res);
        err  = errnorm(lt, fst, res);

        assert((err < ti + ti*nfs) &&
           "map(diff, i, xs0)"); 
    }

    free(fs);
    free(xs0);
    free(xst);
    free(fst);
    free(res);

    loci_free(&ip);

    OUT(0, "Pass\n");
}

/* Test out of bounds behaviour. */
void boundary1D()  {
    OUT(0, "1D Boundary Test");

    /* See basic1D() */
    double big  = 1e6;
    double fr   = 89./97;

    Range rg    = {.x0 = 1, .dx = 0.3, .len = 20};
    int l       = rg.len;
    double xe   = rg.x0 + l*rg.dx;

    double *xs0     = malloc(l*sizeof(double));
    double *fs      = malloc(l*sizeof(double));
    double *fs_x    = malloc(l*sizeof(double));

    r2a(&rg, xs0);
    mapf(f4, l, xs0, fs);
    mapf(f4_x, l, xs0, fs_x);

    /* Points out of bounds */
    double xts[6]   = {xe    + 1./big,
                       rg.x0 - 1./big,
                       xe    + fr*rg.dx,
                       rg.x0 - fr*rg.dx,
                       xe    + big,
                       rg.x0 - big
    };

    /* Create interpolation */
    Interpolation ip    = interpolation(&rg, fs, fs_x);

    for(int i = 0; i < 6; i++)  {
        double xi   = xts[i];

        /* Test whether out of bounds == nan */
        assert(isnan(interpolate(&ip, xi)) &&
               "interpolate(xi) out of bounds"); 

        assert(isnan(diff_x(&ip, xi)) &&
               "diff_x(xi) out of bounds"); 

        assert(isnan(diff(&ip, 2, xi)) &&
               "diff(2, xi) out of bounds"); 

    }

    free(xs0);
    free(fs);
    free(fs_x);

    loci_free(&ip);

    OUT(0, "Pass\n");
}
