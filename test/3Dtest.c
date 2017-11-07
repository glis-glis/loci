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

/* Map function `f` on grid defined by ranges `rx`, `ry`, `rz` and write 
 * to `out`. */
static void mapgrid(double (*f)(double x, double y, double z), const Range *rx, 
                    const Range *ry, const Range *rz, double out[]) {
    for(int i = 0; i < rx->len;i++) {
     for(int j = 0; j < ry->len; j++) {
      for(int k = 0; k < rz->len; k++) {
        int ijk     = i*ry->len*rz->len + j*rz->len + k;
        double xi   = rx->x0 + rx->dx*i;
        double yj   = ry->x0 + ry->dx*j;
        double zk   = rz->x0 + rz->dx*k;
        out[ijk]    = f(xi, yj, zk);
      }
     }
    }
}

/* Convert ranges `rx`, `ry`, `rz` to arrays `xs`, `ys`, `zs`. */
static void rs2as(const Range *rx, const Range *ry, const Range *rz, 
           double xs[], double ys[], double zs[]) {
    for(int i = 0; i < rx->len; i++) {
     for(int j = 0; j < ry->len; j++) {
      for(int k = 0; k < rz->len; k++) {
        int ijk = i*ry->len*rz->len + j*rz->len + k;
        xs[ijk] = rx->x0 + rx->dx*i;
        ys[ijk] = ry->x0 + ry->dx*j;
        zs[ijk] = rz->x0 + rz->dx*k;
      }
     }
    }
}

static double e3(double x, double y, double z) {return exp(x+y+z);}

/* Test functions and derivative. */
static const double A  = 2.;
static const double B  = 0.5;
static const double C  = 0.75;

static double f1(double x, double y, double z)    {return x + y + z;}
static double f1_x(double x, double y, double z)    {return 1;}
static double f1_y(double x, double y, double z)    {return 1;}
static double f1_z(double x, double y, double z)    {return 1;}
static double f1_xy(double x, double y, double z)    {return 0;}
static double f1_xz(double x, double y, double z)    {return 0;}
static double f1_yz(double x, double y, double z)    {return 0;}
static double f1_xyz(double x, double y, double z)    {return 0;}

static double f2(double x, double y, double z)    {return x*y*z;}
static double f2_x(double x, double y, double z)    {return y*z;}
static double f2_y(double x, double y, double z)    {return x*z;}
static double f2_z(double x, double y, double z)    {return x*y;}
static double f2_xy(double x, double y, double z)    {return z;}
static double f2_xz(double x, double y, double z)    {return y;}
static double f2_yz(double x, double y, double z)    {return x;}
static double f2_xyz(double x, double y, double z)    {return 1;}

static double f3(double x, double y, double z)   {return x*x*y*y*z*z;}
static double f3_x(double x, double y, double z)   {return 2*x*y*y*z*z;}
static double f3_y(double x, double y, double z)   {return 2*x*x*y*z*z;}
static double f3_z(double x, double y, double z)   {return 2*x*x*y*y*z;}
static double f3_xy(double x, double y, double z)   {return 4*x*y*z*z;}
static double f3_xz(double x, double y, double z)   {return 4*x*y*y*z;}
static double f3_yz(double x, double y, double z)   {return 4*x*x*y*z;}
static double f3_xyz(double x, double y, double z)   {return 8*x*y*z;}

static double f4(double x, double y, double z)   {
    return log(A*x*x + B*y*y + C*z*z + 1);
}
static double f4_x(double x, double y, double z)   {
    return 2*A*x/(A*x*x + B*y*y + C*z*z + 1);
}
static double f4_y(double x, double y, double z)   {
    return 2*B*y/(A*x*x + B*y*y + C*z*z + 1);
}
static double f4_z(double x, double y, double z)   {
    return 2*C*z/(A*x*x + B*y*y + C*z*z + 1);
}
static double f4_xy(double x, double y, double z)   {
    return -4*A*B*x*y/SQ(A*x*x + B*y*y + C*z*z + 1);
}
static double f4_xz(double x, double y, double z)   {
    return -4*A*C*x*z/SQ(A*x*x + B*y*y + C*z*z + 1);
}
static double f4_yz(double x, double y, double z)   {
    return -4*B*C*y*z/SQ(A*x*x + B*y*y + C*z*z + 1);
}
static double f4_xyz(double x, double y, double z)   {
    return 16*A*B*C*x*y*z/CU(A*x*x + B*y*y + C*z*z + 1);
}

/* Test function arrays. */
static double (*funs[4])(double x, double y, double z)     = 
    {f1, f2, f3, f4};
static double (*funs_x[4])(double x, double y, double z)   = 
    {f1_x, f2_x, f3_x, f4_x};
static double (*funs_y[4])(double x, double y, double z)   = 
    {f1_y, f2_y, f3_y, f4_y};
static double (*funs_z[4])(double x, double y, double z)   = 
    {f1_z, f2_z, f3_z, f4_z};
static double (*funs_xy[4])(double x, double y, double z)  = 
    {f1_xy, f2_xy, f3_xy, f4_xy};
static double (*funs_xz[4])(double x, double y, double z)  = 
    {f1_xz, f2_xz, f3_xz, f4_xz};
static double (*funs_yz[4])(double x, double y, double z)  = 
    {f1_yz, f2_yz, f3_yz, f4_yz};
static double (*funs_xyz[4])(double x, double y, double z) = 
    {f1_xyz, f2_xyz, f3_xyz, f4_xyz};

/* Test basic functionality:
 * create interpolation
 * interpolate a point
 * interpolate derivative
 * map interpolation
 * map derivation. */
void basic3D()  {
    OUT(0, "3D Basic Test");

    /* Interpolation grid (rx*ry*rz) */
    Range rx    = {.x0 = 1, .dx = 0.3, .len = 20};
    Range ry    = {.x0 = -5, .dx = 0.2, .len = 30};
    Range rz    = {.x0 = 2, .dx = 0.1, .len = 10};

    /* Grid points without last one */
    Range r0x   = {rx.x0, rx.dx, rx.len-1};
    Range r0y   = {ry.x0, ry.dx, ry.len-1};
    Range r0z   = {rz.x0, rz.dx, rz.len-1};

    /* Test points */
    Range rtx   = {rx.x0 + rx.dx/3, rx.dx, rx.len-1};
    Range rty   = {ry.x0 + ry.dx/3, ry.dx, ry.len-1};
    Range rtz   = {rz.x0 + rz.dx/3, rz.dx, rz.len-1};

    /* Test values */
    double xt   = rx.x0 + rx.dx*rx.len/2 + rx.dx/7;
    double yt   = ry.x0 + ry.dx*ry.len/2 + ry.dx/7;
    double zt   = rz.x0 + rz.dx*rz.len/2 + rz.dx/7;

    int lx      = rx.len;   
    int ly      = ry.len;   
    int lz      = rz.len;
    int l       = lx*ly*lz;

    int ltx     = rtx.len;  
    int lty     = rty.len;  
    int ltz     = rtz.len;  
    int lt      = ltx*lty*ltz;

    /* allocate memorry for input and outpout values */
    double *fs     = malloc(l*sizeof(double));
    double *fs_x   = malloc(l*sizeof(double));
    double *fs_y   = malloc(l*sizeof(double));
    double *fs_z   = malloc(l*sizeof(double));
    double *fs_xy  = malloc(l*sizeof(double));
    double *fs_xz  = malloc(l*sizeof(double));
    double *fs_yz  = malloc(l*sizeof(double));
    double *fs_xyz = malloc(l*sizeof(double));

    double *xs0     = malloc(lt*sizeof(double));
    double *ys0     = malloc(lt*sizeof(double));
    double *zs0     = malloc(lt*sizeof(double));
    double *xst     = malloc(lt*sizeof(double));
    double *yst     = malloc(lt*sizeof(double));
    double *zst     = malloc(lt*sizeof(double));
    double *res     = malloc(lt*sizeof(double));

    double *fs0    = malloc(lt*sizeof(double));
    double *fs0_x  = malloc(lt*sizeof(double));
    double *fs0_y  = malloc(lt*sizeof(double));
    double *fs0_z  = malloc(lt*sizeof(double));

    /* Transform range to array */
    rs2as(&r0x, &r0y, &r0z, xs0, ys0, zs0);
    rs2as(&rtx, &rty, &rtz, xst, yst, zst);

    for(int i = 0; i < 4; i++)  {
        OUT(1, "Function %i", i+1);

        /* Exact values */
        double ft       = funs[i](xt, yt, zt);
        double ft_x     = funs_x[i](xt, yt, zt);
        double ft_y     = funs_y[i](xt, yt, zt);
        double ft_z     = funs_z[i](xt, yt, zt);

        /* Map exact values on grid points */
        mapgrid(funs[i], &rx, &ry, &rz, fs);
        mapgrid(funs_x[i], &rx, &ry, &rz, fs_x);
        mapgrid(funs_y[i], &rx, &ry, &rz, fs_y);
        mapgrid(funs_z[i], &rx, &ry, &rz, fs_z);
        mapgrid(funs_xy[i], &rx, &ry, &rz, fs_xy);
        mapgrid(funs_xz[i], &rx, &ry, &rz, fs_xz);
        mapgrid(funs_yz[i], &rx, &ry, &rz, fs_yz);
        mapgrid(funs_xyz[i], &rx, &ry, &rz, fs_xyz);

        /* Create interpolation */
        Interpolation ip    = interpolation(&rx, &ry, &rz, fs, 
                                                fs_x, fs_y, fs_z, 
                                                fs_xy, fs_xz, fs_yz, fs_xyz);

        /* Compare interpolated value to exact value */
        double err  = fabs(interpolate(&ip, xt, yt, zt) - ft);
        assert((err < tols[1] + tols[1]*fabs(ft)) &&
               "interpolation(xt, yt, zt)"); 

        /* Compare interpolated derivative to exact value */
        err  = fabs(diff_x(&ip, xt, yt, zt) - ft_x);
        assert((err < tols[2] + tols[2]*fabs(ft_x)) &&
               "diff_x(xt, yt, zt)"); 

        err  = fabs(diff_y(&ip, xt, yt, zt) - ft_y);
        assert((err < tols[2] + tols[2]*fabs(ft_y)) &&
               "diff_y(xt, yt, zt)"); 

        err  = fabs(diff_z(&ip, xt, yt, zt) - ft_z);
        assert((err < tols[2] + tols[2]*fabs(ft_z)) &&
               "diff_z(xt, yt, zt)"); 

        /**/

        /* Map exact values on grid points (without last one,
         * as it is out of bounds) */
        mapgrid(funs[i], &r0x, &r0y, &r0z, fs0);
        mapgrid(funs_x[i], &r0x, &r0y, &r0z, fs0_x);
        mapgrid(funs_y[i], &r0x, &r0y, &r0z, fs0_y);
        mapgrid(funs_z[i], &r0x, &r0y, &r0z, fs0_z);

        double nfs  = norm(lt, fs0);
        double nfs_x= norm(lt, fs0_x);
        double nfs_y= norm(lt, fs0_y);
        double nfs_z= norm(lt, fs0_z);

        /* Compare interpolated values to exact values on grid-points */
        map(&ip, lt, xs0, ys0, zs0, res);
        err  = errnorm(lt, fs0, res);
        assert((err < tols[0] + tols[0]*nfs) &&
               "map(interpolation, xs0, ys0, zs0)"); 

        /* Compare interpolated derivatives to exact values on grid-points */
        map_x(&ip, lt, xs0, ys0, zs0, res);
        err  = errnorm(lt, fs0_x, res);
        assert((err < tols[0] + tols[0]*nfs_x) &&
               "map(diff_x, xs0, ys0, zs0)"); 

        map_y(&ip, lt, xs0, ys0, zs0, res);
        err  = errnorm(lt, fs0_y, res);
        assert((err < tols[0] + tols[0]*nfs_y) &&
               "map(diff_y, xs0, ys0, zs0)"); 

        map_z(&ip, lt, xs0, ys0, zs0, res);
        err  = errnorm(lt, fs0_z, res);
        assert((err < tols[0] + tols[0]*nfs_z) &&
               "map(diff_z, xs0, ys0, zs0)"); 

        /**/

        /* Map exact values on test points */
        mapgrid(funs[i], &rtx, &rty, &rtz, fs0);
        mapgrid(funs_x[i], &rtx, &rty, &rtz, fs0_x);
        mapgrid(funs_y[i], &rtx, &rty, &rtz, fs0_y);
        mapgrid(funs_z[i], &rtx, &rty, &rtz, fs0_z);

        nfs     = norm(lt, fs0);
        nfs_x   = norm(lt, fs0_x);
        nfs_y   = norm(lt, fs0_y);
        nfs_z   = norm(lt, fs0_z);

        /* Compare interpolated values to exact values on test-points */
        map(&ip, lt, xst, yst, zst, res);
        err  = errnorm(lt, fs0, res);
        assert((err < tols[1] + tols[1]*nfs) &&
               "map(interpolation, xst, yst, zst)"); 

        /* Compare interpolated derivatives to exact values on test-points */
        map_x(&ip, lt, xst, yst, zst, res);
        err  = errnorm(lt, fs0_x, res);
        assert((err < tols[1] + tols[1]*nfs_x) &&
               "map(diff_x, xst, yst, zst)"); 

        map_y(&ip, lt, xst, yst, zst, res);
        err  = errnorm(lt, fs0_y, res);
        assert((err < tols[1] + tols[1]*nfs_y) &&
               "map(diff_y, xst, yst, zst)"); 

        map_z(&ip, lt, xst, yst, zst, res);
        err  = errnorm(lt, fs0_z, res);
        assert((err < tols[1] + tols[1]*nfs_z) &&
               "map(diff_z, xst, yst, zst)"); 

        loci_free(&ip);
    }

    free(fs);
    free(fs_x);
    free(fs_y);
    free(fs_z);
    free(fs_xy);
    free(fs_xz);
    free(fs_yz);
    free(fs_xyz);

    free(xs0);
    free(ys0);
    free(zs0);
    free(xst);
    free(yst);
    free(zst);
    free(res);

    free(fs0);
    free(fs0_x);
    free(fs0_y);
    free(fs0_z);

    OUT(0, "Pass\n");
}

/* Test higher order derivative on exponantial function. */
void exp3D()  {
    OUT(0, "3D Exp Test");

    /* See basic3D() */
    Range rx    = {.x0 = 1, .dx = 0.3, .len = 30};
    Range ry    = {.x0 = -5, .dx = 0.2, .len = 20};
    Range rz    = {.x0 = 2, .dx = 0.1, .len = 10};

    Range rtx   = {rx.x0 + rx.dx/3, rx.dx, rx.len-1};
    Range rty   = {ry.x0 + ry.dx/3, ry.dx, ry.len-1};
    Range rtz   = {rz.x0 + rz.dx/3, rz.dx, rz.len-1};

    double xt   = rx.x0 + rx.dx*rx.len/2 + rx.dx/7;
    double yt   = ry.x0 + ry.dx*ry.len/2 + ry.dx/7;
    double zt   = rz.x0 + rz.dx*rz.len/2 + rz.dx/7;

    int lx      = rx.len;   
    int ly      = ry.len;   
    int lz      = rz.len;
    int l       = lx*ly*lz;

    int ltx     = rtx.len;  
    int lty     = rty.len;  
    int ltz     = rtz.len;  
    int lt      = ltx*lty*ltz;

    double *fs     = malloc(l*sizeof(double));
    double *fst    = malloc(lt*sizeof(double));
    double *xst    = malloc(lt*sizeof(double));
    double *yst    = malloc(lt*sizeof(double));
    double *zst    = malloc(lt*sizeof(double));
    double *res    = malloc(lt*sizeof(double));

    mapgrid(e3, &rx, &ry, &rz, fs);

    mapgrid(e3, &rtx, &rty, &rtz, fst);
    double nfs  = norm(lt, fst);

    /* Create interpolation */
    Interpolation ip    = interpolation(&rx, &ry, &rz, fs, fs, fs, fs,
                                            fs, fs, fs, fs);

    rs2as(&rtx, &rty, &rtz, xst, yst, zst);

    double ft   = e3(xt, yt, zt);


    for(int i = 0; i < 3; i++)  {
     for(int j = 0; j < 3; j++)  {
      for(int k = 0; k < 3; k++)  {
        OUT(1, "d^%if/dx^%i/dy^%i/dz^%i", i+j+k, i, j, k);
        double ti   = tols[MAX(i, MAX(j, k)) + 2];

        /* Compare interpolated derivative to exact value */
        double err  = fabs(diff(&ip, i, j, k, xt, yt, zt) - ft);
        assert((err < ti + ti*fabs(ft)) &&
               "diff(i, j, k, xt, yt, zt)"); 

        /* Compare interpolated derivatives to exact values on test-points */
        map_diff(&ip, i, j, k, lt, xst, yst, zst, res);
        err  = errnorm(lt, fst, res);

        assert((err < ti + ti*nfs) &&
           "map(diff, i, j, k, xst, yst, zst)"); 
      }
     }
    }

    free(fs);
    free(fst);
    free(xst);
    free(yst);
    free(zst);
    free(res);

    loci_free(&ip);

    OUT(0, "Pass\n");
}


/* Test out of bounds behaviour. */
void boundary3D()  {
    OUT(0, "3D Boundary Test");

    double big  = 1e6;
    double fr   = 89./97;

    /* See basic3D() */
    Range rx    = {.x0 = 1, .dx = 0.3, .len = 30};
    Range ry    = {.x0 = -5, .dx = 0.2, .len = 20};
    Range rz    = {.x0 = 2, .dx = 0.1, .len = 10};

    int lx      = rx.len;   
    int ly      = ry.len;   
    int lz      = rz.len;
    int l       = lx*ly*lz;

    double xe   = rx.x0 + lx*rx.dx;
    double ye   = ry.x0 + ly*ry.dx;
    double ze   = rz.x0 + lz*rz.dx;

    double *fs      = malloc(l*sizeof(double));
    double *fs_x    = malloc(l*sizeof(double));
    double *fs_y    = malloc(l*sizeof(double));
    double *fs_z    = malloc(l*sizeof(double));
    double *fs_xy   = malloc(l*sizeof(double));
    double *fs_xz   = malloc(l*sizeof(double));
    double *fs_yz   = malloc(l*sizeof(double));
    double *fs_xyz  = malloc(l*sizeof(double));

    mapgrid(f4, &rx, &ry, &rz, fs);
    mapgrid(f4_x, &rx, &ry, &rz, fs_x);
    mapgrid(f4_y, &rx, &ry, &rz, fs_y);
    mapgrid(f4_z, &rx, &ry, &rz, fs_z);
    mapgrid(f4_xy, &rx, &ry, &rz, fs_xy);
    mapgrid(f4_xz, &rx, &ry, &rz, fs_xz);
    mapgrid(f4_yz, &rx, &ry, &rz, fs_yz);
    mapgrid(f4_xyz, &rx, &ry, &rz, fs_xyz);

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

    double zts[7]   = {rz.x0 + fr*rz.dx,    /* in bounds */
                       ze    + 1./big,
                       rz.x0 - 1./big,
                       ze    + fr*rz.dx,
                       rz.x0 - fr*rz.dx,
                       ze    + big,
                       rz.x0 - big
    };

    /* Create interpolation */
    Interpolation ip    = interpolation(&rx, &ry, &rz, fs, 
                                           fs_x, fs_y, fs_z,
                                           fs_xy, fs_xz, fs_yz, fs_xyz);

    for(int i = 0; i < 7; i++)  {
     for(int j = 0; j < 7; j++)  {
      for(int k = 0; k < 7; k++)  {
        if(i+j+k > 0)    {
            double xi   = xts[i];
            double yj   = yts[j];
            double zk   = zts[k];

            /* Test whether out of bounds == nan */
            assert(isnan(interpolate(&ip, xi, yj, zk)) &&
                   "interpolate(xi, yj, zk) out of bounds"); 

            assert(isnan(diff_x(&ip, xi, yj, zk)) &&
                   "diff_x(xi, yj, zk) out of bounds"); 

            assert(isnan(diff_y(&ip, xi, yj, zk)) &&
                   "diff_y(xi, yj, zk) out of bounds"); 

            assert(isnan(diff_z(&ip, xi, yj, zk)) &&
                   "diff_z(xi, yj, zk) out of bounds"); 

            assert(isnan(diff(&ip, 2, 2, 2, xi, yj, zk)) &&
                   "diff(2, 2, 2, xi, yj, zk) out of bounds"); 
        }
      }
     }
    }

    free(fs);
    free(fs_x);
    free(fs_y);
    free(fs_z);
    free(fs_xy);
    free(fs_xz);
    free(fs_yz);
    free(fs_xyz);

    loci_free(&ip);

    OUT(0, "Pass\n");
}
