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

/* Map function `f` on grid defined by ranges `rx`, `ry` and write to `out`. */
static void mapgrid(double (*f)(double x, double y), const Range *rx, 
               const Range *ry, double out[]) {
    for(int i = 0; i < rx->len; i++) {
        for(int j = 0; j < ry->len; j++) {
            int ij      = i*ry->len + j;
            double xi   = rx->x0 + rx->dx*i;
            double yj   = ry->x0 + ry->dx*j;
            out[ij]     = f(xi, yj);
        }
    }
}

/* Convert ranges `rx`, `ry` to arrays `xs`, `ys`. */
static void rs2as(const Range *rx, const Range *ry, double xs[], double ys[]) {
    for(int i = 0; i < rx->len; i++) {
        for(int j = 0; j < ry->len; j++) {
            int ij  = i*ry->len + j;
            xs[ij]  = rx->x0 + rx->dx*i;
            ys[ij]  = ry->x0 + ry->dx*j;
        }
    }
}

static double e2(double x, double y) {return exp(x+y);}

/* Test functions and derivative. */
static const double A  = 2.;
static const double B  = 0.5;

static double f1(double x, double y)    {return x + y;}
static double f1_x(double x, double y)  {return 1;}
static double f1_y(double x, double y)  {return 1;}
static double f1_xy(double x, double y) {return 0;}

static double f2(double x, double y)    {return x*y;}
static double f2_x(double x, double y)  {return y;}
static double f2_y(double x, double y)  {return x;}
static double f2_xy(double x, double y) {return 1;}

static double f3(double x, double y)    {return x*x*y*y;}
static double f3_x(double x, double y)  {return 2*x*y*y;}
static double f3_y(double x, double y)  {return 2*x*x*y;}
static double f3_xy(double x, double y) {return 4*x*y;}

static double f4(double x, double y)    {return log(A*x*x + B*y*y + 1);}
static double f4_x(double x, double y)  {return 2*A*x/(A*x*x + B*y*y + 1);}
static double f4_y(double x, double y)  {return 2*B*y/(A*x*x + B*y*y + 1);}
static double f4_xy(double x, double y) {
    return -4*A*B*x*y/SQ(A*x*x + B*y*y + 1);
}

/* Test function arrays. */
static double (*funs[4])(double x, double y)    = {f1, f2, f3, f4};
static double (*funs_x[4])(double x, double y)  = {f1_x, f2_x, f3_x, f4_x};
static double (*funs_y[4])(double x, double y)  = {f1_y, f2_y, f3_y, f4_y};
static double (*funs_xy[4])(double x, double y) = {f1_xy, f2_xy, f3_xy, f4_xy};

/* Test basic functionality:
 * create interpolation
 * interpolate a point
 * interpolate derivative
 * map interpolation
 * map derivation. */
void basic2D()  {
    OUT(0, "2D Basic Test");

    /* Interpolation grid (rx*ry) */
    Range rx    = {.x0 = 1, .dx = 0.3, .len = 20};
    Range ry    = {.x0 = -5, .dx = 0.2, .len = 30};

    /* Grid points without last one */
    Range r0x   = {rx.x0, rx.dx, rx.len-1};
    Range r0y   = {ry.x0, ry.dx, ry.len-1};

    /* Test points */
    Range rtx   = {rx.x0 + rx.dx/3, rx.dx, rx.len-1};
    Range rty   = {ry.x0 + ry.dx/3, ry.dx, ry.len-1};

    /* Test values */
    double xt   = rx.x0 + rx.dx*rx.len/2 + rx.dx/7;
    double yt   = ry.x0 + ry.dx*ry.len/2 + ry.dx/7;

    int lx      = rx.len;   
    int ly      = ry.len;   

    int ltx     = rtx.len;  
    int lty     = rty.len;  
    int lt      = ltx*lty;

    /* allocate memorry for input and outpout values */
    double *fs     = malloc(lx*ly*sizeof(double));
    double *fs_x   = malloc(lx*ly*sizeof(double));
    double *fs_y   = malloc(lx*ly*sizeof(double));
    double *fs_xy  = malloc(lx*ly*sizeof(double));

    double *xs0     = malloc(lt*sizeof(double));
    double *ys0     = malloc(lt*sizeof(double));
    double *xst     = malloc(lt*sizeof(double));
    double *yst     = malloc(lt*sizeof(double));
    double *res     = malloc(lt*sizeof(double));

    double *fs0    = malloc(lt*sizeof(double));
    double *fs0_x  = malloc(lt*sizeof(double));
    double *fs0_y  = malloc(lt*sizeof(double));

    /* Transorm range to array */
    rs2as(&r0x, &r0y, xs0, ys0);
    rs2as(&rtx, &rty, xst, yst);

    for(int i = 0; i < 4; i++)  {
        OUT(1, "Function %i", i+1);

        /* Exact values */
        double ft       = funs[i](xt, yt);
        double ft_x     = funs_x[i](xt, yt);
        double ft_y     = funs_y[i](xt, yt);

        /* Map exact values on grid points */
        mapgrid(funs[i], &rx, &ry, fs);
        mapgrid(funs_x[i], &rx, &ry, fs_x);
        mapgrid(funs_y[i], &rx, &ry, fs_y);
        mapgrid(funs_xy[i], &rx, &ry, fs_xy);

        /* Create interpolation */
        Interpolation ip    = interpolation(&rx, &ry, fs, fs_x, fs_y,
                                                fs_xy);

        /* Compare interpolated value to exact value */
        double err  = fabs(interpolate(&ip, xt, yt) - ft);
        assert((err < tols[1] + tols[1]*fabs(ft)) &&
               "interpolation(xt, yt)"); 

        /* Compare interpolated derivative to exact value */
        err  = fabs(diff_x(&ip, xt, yt) - ft_x);
        assert((err < tols[2] + tols[2]*fabs(ft_x)) &&
               "diff_x(xt, yt)"); 

        err  = fabs(diff_y(&ip, xt, yt) - ft_y);
        assert((err < tols[2] + tols[2]*fabs(ft_y)) &&
               "diff_y(xt, yt)"); 

        /**/


        /* Map exact values on grid points (without last one,
         * as it is out of bounds) */
        mapgrid(funs[i], &r0x, &r0y, fs0);
        mapgrid(funs_x[i], &r0x, &r0y, fs0_x);
        mapgrid(funs_y[i], &r0x, &r0y, fs0_y);

        double nfs  = norm(lt, fs0);
        double nfs_x= norm(lt, fs0_x);
        double nfs_y= norm(lt, fs0_y);

        /* Compare interpolated values to exact values on grid-points */
        map(&ip, lt, xs0, ys0, res);
        err  = errnorm(lt, fs0, res);
        assert((err < tols[0] + tols[0]*nfs) &&
               "map(interpolation, xs0, ys0)"); 

        /* Compare interpolated derivatives to exact values on grid-points */
        map_x(&ip, lt, xs0, ys0, res);
        err  = errnorm(lt, fs0_x, res);
        assert((err < tols[0] + tols[0]*nfs_x) &&
               "map(diff_x, xs0, ys0)"); 

        map_y(&ip, lt, xs0, ys0, res);
        err  = errnorm(lt, fs0_y, res);
        assert((err < tols[0] + tols[0]*nfs_y) &&
               "map(diff_y, xs0, ys0)"); 

        /**/

        /* Map exact values on test points */
        mapgrid(funs[i], &rtx, &rty, fs0);
        mapgrid(funs_x[i], &rtx, &rty, fs0_x);
        mapgrid(funs_y[i], &rtx, &rty, fs0_y);

        nfs     = norm(lt, fs0);
        nfs_x   = norm(lt, fs0_x);
        nfs_y   = norm(lt, fs0_y);

        /* Compare interpolated values to exact values on test-points */
        map(&ip, lt, xst, yst, res);
        err  = errnorm(lt, fs0, res);
        assert((err < tols[1] + tols[1]*nfs) &&
               "map(interpolation, xst, yst)"); 

        /* Compare interpolated derivatives to exact values on test-points */
        map_x(&ip, lt, xst, yst, res);
        err  = errnorm(lt, fs0_x, res);
        assert((err < tols[2] + tols[2]*nfs_x) &&
               "map(diff_x, xst, yst)"); 

        map_y(&ip, lt, xst, yst, res);
        err  = errnorm(lt, fs0_y, res);
        assert((err < tols[2] + tols[2]*nfs_y) &&
               "map(diff_y, xst, yst)"); 

        loci_free(&ip);
    }



    free(fs);
    free(fs_x);
    free(fs_y);
    free(fs_xy);

    free(xs0);
    free(ys0);
    free(xst);
    free(yst);
    free(res);

    free(fs0);
    free(fs0_x);
    free(fs0_y);

    OUT(0, "Pass\n");
}

/* Test higher order derivative on exponantial function. */
void exp2D()  {
    OUT(0, "2D Exp Test");

    /* See basic2D() */
    Range rx    = {.x0 = 1, .dx = 0.3, .len = 20};
    Range ry    = {.x0 = -5, .dx = 0.2, .len = 30};

    Range rtx   = {rx.x0 + rx.dx/3, rx.dx, rx.len-1};
    Range rty   = {ry.x0 + ry.dx/3, ry.dx, ry.len-1};

    double xt   = rx.x0 + rx.dx*rx.len/2 + rx.dx/7;
    double yt   = ry.x0 + ry.dx*ry.len/2 + ry.dx/7;

    int lx      = rx.len;   
    int ly      = ry.len;   

    int ltx     = rtx.len;  
    int lty     = rty.len;  
    int lt      = ltx*lty;

    double *fs     = malloc(lx*ly*sizeof(double));
    double *fst    = malloc(lt*sizeof(double));
    double *xst    = malloc(lt*sizeof(double));
    double *yst    = malloc(lt*sizeof(double));
    double *res    = malloc(lt*sizeof(double));

    mapgrid(e2, &rx, &ry, fs);

    mapgrid(e2, &rtx, &rty, fst);
    double nfs  = norm(lt, fst);

    /* Create interpolation */
    Interpolation ip    = interpolation(&rx, &ry, fs, fs, fs, fs);

    rs2as(&rtx, &rty, xst, yst);

    double ft   = e2(xt, yt);

    for(int i = 0; i < 4; i++)  {
        for(int j = 0; j < 4; j++)  {
            OUT(1, "d^%if/dx^%i/dy^%i", i+j, i, j);
            double ti   = tols[MAX(i, j) + 1];

            /* Compare interpolated derivative to exact value */
            double err  = fabs(diff(&ip, i, j, xt, yt) - ft);
            assert((err < ti + ti*fabs(ft)) &&
                   "diff(i, j, xt, yt)"); 

            /* Compare interpolated derivatives to exact values on 
             * test-points */
            map_diff(&ip, i, j, lt, xst, yst, res);
            err  = errnorm(lt, fst, res);

            assert((err < ti + ti*nfs) &&
               "map(diff, i, j, xst, yst)"); 
        }
    }

    free(fs);
    free(fst);
    free(xst);
    free(yst);
    free(res);

    loci_free(&ip);

    OUT(0, "Pass\n");
}

/* Test out of bounds behaviour. */
void boundary2D()  {
    OUT(0, "2D Boundary Test");

    double big  = 1e6;
    double fr   = 89./97;

    /* See basic2D() */
    Range rx    = {.x0 = 1, .dx = 0.3, .len = 20};
    Range ry    = {.x0 = -5, .dx = 0.2, .len = 30};

    int lx      = rx.len;   
    int ly      = ry.len;   
    int l       = lx*ly;
    double xe   = rx.x0 + lx*rx.dx;
    double ye   = ry.x0 + ly*ry.dx;

    double *fs      = malloc(l*sizeof(double));
    double *fs_x    = malloc(l*sizeof(double));
    double *fs_y    = malloc(l*sizeof(double));
    double *fs_xy   = malloc(l*sizeof(double));

    mapgrid(f4, &rx, &ry, fs);
    mapgrid(f4_x, &rx, &ry, fs_x);
    mapgrid(f4_y, &rx, &ry, fs_y);
    mapgrid(f4_xy, &rx, &ry, fs_xy);

    /* Points out of bounds */
    double xts[7]   = {rx.x0 + fr*rx.dx,    /* in bounds */
                       xe    - 1./big,
                       rx.x0 - 1./big,
                       xe    + fr*rx.dx,
                       rx.x0 - fr*rx.dx,
                       xe    + big,
                       rx.x0 - big
    };

    double yts[7]   = {ry.x0 + fr*ry.dx,    /* in bounds */
                       ye    + 1./big,
                       ry.x0 - 1./big,
                       ye    + fr*ry.dx,
                       ry.x0 - fr*ry.dx,
                       ye    + big,
                       ry.x0 - big
    };

    /* Create interpolation */
    Interpolation ip    = interpolation(&rx, &ry, fs, fs_x, fs_y, fs_xy);

    for(int i = 0; i < 7; i++)  {
        for(int j = 0; j < 7; j++)  {
            if(i+j > 0)    {
                double xi   = xts[i];
                double yj   = yts[j];

                /* Test whether out of bounds == nan */
                assert(isnan(interpolate(&ip, xi, yj)) &&
                       "interpolate(xi, yj) out of bounds"); 

                assert(isnan(diff_x(&ip, xi, yj)) &&
                       "diff_x(xi, yj) out of bounds"); 

                assert(isnan(diff_y(&ip, xi, yj)) &&
                       "diff_z(xi, yj) out of bounds"); 

                assert(isnan(diff(&ip, 2, 2, xi, yj)) &&
                       "diff(2, 2, xi, yj) out of bounds"); 
            }
        }

    }

    free(fs);
    free(fs_x);
    free(fs_y);
    free(fs_xy);

    loci_free(&ip);

    OUT(0, "Pass\n");
}
