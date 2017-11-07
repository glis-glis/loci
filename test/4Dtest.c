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

/* Map function `f` on grid defined by ranges `rx`, `ry`, `rz`, `rw` and write 
 * to `out`. */
static void mapgrid(double (*f)(double x, double y, double z, double w), 
                    const Range *rx, const Range *ry, const Range *rz, 
                    const Range *rw, double out[]) {
    for(int i = 0; i < rx->len;i++) {
     for(int j = 0; j < ry->len; j++) {
      for(int k = 0; k < rz->len; k++) {
       for(int l = 0; l < rw->len; l++) {
        int ijkl    = i*ry->len*rz->len*rw->len + j*rz->len*rw->len + 
                      k*rw->len + l;
        double xi   = rx->x0 + rx->dx*i;
        double yj   = ry->x0 + ry->dx*j;
        double zk   = rz->x0 + rz->dx*k;
        double wl   = rw->x0 + rw->dx*l;
        out[ijkl]   = f(xi, yj, zk, wl);
       }
      }
     }
    }
}

/* Convert ranges `rx`, `ry`, `rz`, `rw` to arrays `xs`, `ys`, `zs`, `ws`. */
static void rs2as(const Range *rx, const Range *ry, 
                  const Range *rz, const Range *rw, 
                  double xs[], double ys[], double zs[], double ws[]) {
    for(int i = 0; i < rx->len; i++) {
     for(int j = 0; j < ry->len; j++) {
      for(int k = 0; k < rz->len; k++) {
       for(int l = 0; l < rw->len; l++) {
        int ijkl    = i*ry->len*rz->len*rw->len + j*rz->len*rw->len + 
                      k*rw->len + l;
        xs[ijkl]    = rx->x0 + rx->dx*i;
        ys[ijkl]    = ry->x0 + ry->dx*j;
        zs[ijkl]    = rz->x0 + rz->dx*k;
        ws[ijkl]    = rw->x0 + rw->dx*l;
       }
      }
     }
    }
}

static double e4(double x, double y, double z, double w) {return exp(x+y+z+w);}

/* Test functions */
static const double A  = 2;
static const double B  = 0.5;
static const double C  = 0.75;
static const double D  = 0.4;

static double f1(double x, double y, double z, double w) {return x + y + z + w;}
static double f1_x(double x, double y, double z, double w)    {return 1;}
static double f1_y(double x, double y, double z, double w)    {return 1;}
static double f1_z(double x, double y, double z, double w)    {return 1;}
static double f1_w(double x, double y, double z, double w)    {return 1;}
static double f1_xy(double x, double y, double z, double w)    {return 0;}
static double f1_xz(double x, double y, double z, double w)    {return 0;}
static double f1_xw(double x, double y, double z, double w)    {return 0;}
static double f1_yz(double x, double y, double z, double w)    {return 0;}
static double f1_yw(double x, double y, double z, double w)    {return 0;}
static double f1_zw(double x, double y, double z, double w)    {return 0;}
static double f1_xyz(double x, double y, double z, double w)    {return 0;}
static double f1_xyw(double x, double y, double z, double w)    {return 0;}
static double f1_xzw(double x, double y, double z, double w)    {return 0;}
static double f1_yzw(double x, double y, double z, double w)    {return 0;}
static double f1_xyzw(double x, double y, double z, double w)    {return 0;}

static double f2(double x, double y, double z, double w)    {return x*y*z*w;}
static double f2_x(double x, double y, double z, double w)    {return y*z*w;}
static double f2_y(double x, double y, double z, double w)    {return x*z*w;}
static double f2_z(double x, double y, double z, double w)    {return x*y*w;}
static double f2_w(double x, double y, double z, double w)    {return x*y*z;}
static double f2_xy(double x, double y, double z, double w)    {return z*w;}
static double f2_xz(double x, double y, double z, double w)    {return y*w;}
static double f2_xw(double x, double y, double z, double w)    {return y*z;}
static double f2_yz(double x, double y, double z, double w)    {return x*w;}
static double f2_yw(double x, double y, double z, double w)    {return x*z;}
static double f2_zw(double x, double y, double z, double w)    {return x*y;}
static double f2_xyz(double x, double y, double z, double w)    {return w;}
static double f2_xyw(double x, double y, double z, double w)    {return z;}
static double f2_xzw(double x, double y, double z, double w)    {return y;}
static double f2_yzw(double x, double y, double z, double w)    {return x;}
static double f2_xyzw(double x, double y, double z, double w)    {return 1;}

static double f3(double x, double y, double z, double w)   {
    return SQ(x)*SQ(y)*SQ(z)*SQ(w);
}
static double f3_x(double x, double y, double z, double w)   {
    return 2*x*SQ(y)*SQ(z)*SQ(w);
}
static double f3_y(double x, double y, double z, double w)   {
    return 2*y*SQ(x)*SQ(z)*SQ(w);
}
static double f3_z(double x, double y, double z, double w)   {
    return 2*z*SQ(x)*SQ(y)*SQ(w);
}
static double f3_w(double x, double y, double z, double w)   {
    return 2*w*SQ(x)*SQ(y)*SQ(z);
}
static double f3_xy(double x, double y, double z, double w)   {
    return 4*x*y*SQ(z)*SQ(w);
}
static double f3_xz(double x, double y, double z, double w)   {
    return 4*x*z*SQ(y)*SQ(w);
}
static double f3_xw(double x, double y, double z, double w)   {
    return 4*x*w*SQ(y)*SQ(z);
}
static double f3_yz(double x, double y, double z, double w)   {
    return 4*y*z*SQ(x)*SQ(w);
}
static double f3_yw(double x, double y, double z, double w)   {
    return 4*y*w*SQ(x)*SQ(z);
}
static double f3_zw(double x, double y, double z, double w)   {
    return 4*z*w*SQ(x)*SQ(y);
}
static double f3_xyz(double x, double y, double z, double w)   {
    return 8*x*y*z*SQ(w);
}
static double f3_xyw(double x, double y, double z, double w)   {
    return 8*x*y*w*SQ(z);
}
static double f3_xzw(double x, double y, double z, double w)   {
    return 8*x*z*w*SQ(y);
}
static double f3_yzw(double x, double y, double z, double w)   {
    return 8*y*z*w*SQ(x);
}
static double f3_xyzw(double x, double y, double z, double w)   {
    return 16*x*y*z*w;
}

static double f4(double x, double y, double z, double w)   {
    return log(A*x*x + B*y*y + C*z*z + D*w*w + 1);
}
static double f4_x(double x, double y, double z, double w)   {
    return 2*A*x/(A*x*x + B*y*y + C*z*z + D*w*w + 1);
}
static double f4_y(double x, double y, double z, double w)   {
    return 2*B*y/(A*x*x + B*y*y + C*z*z + D*w*w + 1);
}
static double f4_z(double x, double y, double z, double w)   {
    return 2*C*z/(A*x*x + B*y*y + C*z*z + D*w*w + 1);
}
static double f4_w(double x, double y, double z, double w)   {
    return 2*D*w/(A*x*x + B*y*y + C*z*z + D*w*w + 1);
}
static double f4_xy(double x, double y, double z, double w)   {
    return -4*A*B*x*y/SQ(A*x*x + B*y*y + C*z*z + D*w*w + 1);
}
static double f4_xz(double x, double y, double z, double w)   {
    return -4*A*C*x*z/SQ(A*x*x + B*y*y + C*z*z + D*w*w + 1);
}
static double f4_xw(double x, double y, double z, double w)   {
    return -4*A*D*x*w/SQ(A*x*x + B*y*y + C*z*z + D*w*w + 1);
}
static double f4_yz(double x, double y, double z, double w)   {
    return -4*B*C*y*z/SQ(A*x*x + B*y*y + C*z*z + D*w*w + 1);
}
static double f4_yw(double x, double y, double z, double w)   {
    return -4*B*D*y*w/SQ(A*x*x + B*y*y + C*z*z + D*w*w + 1);
}
static double f4_zw(double x, double y, double z, double w)   {
    return -4*C*D*z*w/SQ(A*x*x + B*y*y + C*z*z + D*w*w + 1);
}
static double f4_xyz(double x, double y, double z, double w)   {
    return 16*A*B*C*x*y*z/CU(A*x*x + B*y*y + C*z*z + D*w*w + 1);
}
static double f4_xyw(double x, double y, double z, double w)   {
    return 16*A*B*D*x*y*w/CU(A*x*x + B*y*y + C*z*z + D*w*w + 1);
}
static double f4_xzw(double x, double y, double z, double w)   {
    return 16*A*C*D*x*z*w/CU(A*x*x + B*y*y + C*z*z + D*w*w + 1);
}
static double f4_yzw(double x, double y, double z, double w)   {
    return 16*A*C*D*x*z*w/CU(A*x*x + B*y*y + C*z*z + D*w*w + 1);
}
static double f4_xyzw(double x, double y, double z, double w)   {
    return -96*A*B*C*D*x*y*z*w/QU(A*x*x + B*y*y + C*z*z + D*w*w + 1);
}

/* Test function arrays. */
static double (*funs[4])(double x, double y, double z, double w)     = 
    {f1, f2, f3, f4};

static double (*funs_x[4])(double x, double y, double z, double w)   = 
    {f1_x, f2_x, f3_x, f4_x};
static double (*funs_y[4])(double x, double y, double z, double w)   = 
    {f1_y, f2_y, f3_y, f4_y};
static double (*funs_z[4])(double x, double y, double z, double w)   = 
    {f1_z, f2_z, f3_z, f4_z};
static double (*funs_w[4])(double x, double y, double z, double w)   = 
    {f1_w, f2_w, f3_w, f4_w};

static double (*funs_xy[4])(double x, double y, double z, double w)  = 
    {f1_xy, f2_xy, f3_xy, f4_xy};
static double (*funs_xz[4])(double x, double y, double z, double w)  = 
    {f1_xz, f2_xz, f3_xz, f4_xz};
static double (*funs_xw[4])(double x, double y, double z, double w)  = 
    {f1_xw, f2_xw, f3_xw, f4_xw};
static double (*funs_yz[4])(double x, double y, double z, double w)  = 
    {f1_yz, f2_yz, f3_yz, f4_yz};
static double (*funs_yw[4])(double x, double y, double z, double w)  = 
    {f1_yw, f2_yw, f3_yw, f4_yw};
static double (*funs_zw[4])(double x, double y, double z, double w)  = 
    {f1_zw, f2_zw, f3_zw, f4_zw};

static double (*funs_xyz[4])(double x, double y, double z, double w) = 
    {f1_xyz, f2_xyz, f3_xyz, f4_xyz};
static double (*funs_xyw[4])(double x, double y, double z, double w) = 
    {f1_xyw, f2_xyw, f3_xyw, f4_xyw};
static double (*funs_xzw[4])(double x, double y, double z, double w) = 
    {f1_xzw, f2_xzw, f3_xzw, f4_xzw};
static double (*funs_yzw[4])(double x, double y, double z, double w) = 
    {f1_yzw, f2_yzw, f3_yzw, f4_yzw};

static double (*funs_xyzw[4])(double x, double y, double z, double w) = 
    {f1_xyzw, f2_xyzw, f3_xyzw, f4_xyzw};

/* Test basic functionality:
 * create interpolation
 * interpolate a point
 * interpolate derivative
 * map interpolation
 * map derivation. */
void basic4D()  {
    OUT(0, "4D Basic Test");

    /* Interpolation grid (rx*ry*rz*rw) */
    Range rx    = {.x0 = 1, .dx = 0.3, .len = 20};
    Range ry    = {.x0 = -5, .dx = 0.2, .len = 30};
    Range rz    = {.x0 = 2, .dx = 0.1, .len = 10};
    Range rw    = {.x0 = -3, .dx = 0.4, .len = 15};

    /* Grid points without last one */
    Range r0x   = {rx.x0, rx.dx, rx.len-1};
    Range r0y   = {ry.x0, ry.dx, ry.len-1};
    Range r0z   = {rz.x0, rz.dx, rz.len-1};
    Range r0w   = {rw.x0, rw.dx, rw.len-1};

    /* Test points */
    Range rtx   = {rx.x0 + rx.dx/3, rx.dx, rx.len-1};
    Range rty   = {ry.x0 + ry.dx/3, ry.dx, ry.len-1};
    Range rtz   = {rz.x0 + rz.dx/3, rz.dx, rz.len-1};
    Range rtw   = {rw.x0 + rw.dx/3, rw.dx, rw.len-1};

    /* Test values */
    double xt   = rx.x0 + rx.dx*rx.len/2 + rx.dx/7;
    double yt   = ry.x0 + ry.dx*ry.len/2 + ry.dx/7;
    double zt   = rz.x0 + rz.dx*rz.len/2 + rz.dx/7;
    double wt   = rw.x0 + rw.dx*rw.len/2 + rw.dx/7;

    int lx      = rx.len;   
    int ly      = ry.len;   
    int lz      = rz.len;
    int lw      = rw.len;
    int l       = lx*ly*lz*lw;

    int ltx     = rtx.len;  
    int lty     = rty.len;  
    int ltz     = rtz.len;  
    int ltw     = rtw.len;  
    int lt      = ltx*lty*ltz*ltw;

    /* allocate memorry for input and outpout values */
    double *fs     = malloc(l*sizeof(double));
    double *fs_x   = malloc(l*sizeof(double));
    double *fs_y   = malloc(l*sizeof(double));
    double *fs_z   = malloc(l*sizeof(double));
    double *fs_w   = malloc(l*sizeof(double));
    double *fs_xy  = malloc(l*sizeof(double));
    double *fs_xz  = malloc(l*sizeof(double));
    double *fs_xw  = malloc(l*sizeof(double));
    double *fs_yz  = malloc(l*sizeof(double));
    double *fs_yw  = malloc(l*sizeof(double));
    double *fs_zw  = malloc(l*sizeof(double));
    double *fs_xyz = malloc(l*sizeof(double));
    double *fs_xyw = malloc(l*sizeof(double));
    double *fs_xzw = malloc(l*sizeof(double));
    double *fs_yzw = malloc(l*sizeof(double));
    double *fs_xyzw= malloc(l*sizeof(double));

    double *xs0     = malloc(lt*sizeof(double));
    double *ys0     = malloc(lt*sizeof(double));
    double *zs0     = malloc(lt*sizeof(double));
    double *ws0     = malloc(lt*sizeof(double));
    double *xst     = malloc(lt*sizeof(double));
    double *yst     = malloc(lt*sizeof(double));
    double *zst     = malloc(lt*sizeof(double));
    double *wst     = malloc(lt*sizeof(double));
    double *res     = malloc(lt*sizeof(double));

    double *fs0    = malloc(lt*sizeof(double));
    double *fs0_x  = malloc(lt*sizeof(double));
    double *fs0_y  = malloc(lt*sizeof(double));
    double *fs0_z  = malloc(lt*sizeof(double));
    double *fs0_w  = malloc(lt*sizeof(double));

    /* Transform range to array */
    rs2as(&r0x, &r0y, &r0z, &r0w, xs0, ys0, zs0, ws0);
    rs2as(&rtx, &rty, &rtz, &rtw, xst, yst, zst, wst);

    for(int i = 0; i < 4; i++)  {
        OUT(1, "Function %i", i+1);

        /* Exact values */
        double ft       = funs[i](xt, yt, zt, wt);
        double ft_x     = funs_x[i](xt, yt, zt, wt);
        double ft_y     = funs_y[i](xt, yt, zt, wt);
        double ft_z     = funs_z[i](xt, yt, zt, wt);
        double ft_w     = funs_w[i](xt, yt, zt, wt);

        /* Map exact values on grid points */
        mapgrid(funs[i], &rx, &ry, &rz, &rw, fs);
        mapgrid(funs_x[i], &rx, &ry, &rz, &rw, fs_x);
        mapgrid(funs_y[i], &rx, &ry, &rz, &rw, fs_y);
        mapgrid(funs_z[i], &rx, &ry, &rz, &rw, fs_z);
        mapgrid(funs_w[i], &rx, &ry, &rz, &rw, fs_w);
        mapgrid(funs_xy[i], &rx, &ry, &rz, &rw, fs_xy);
        mapgrid(funs_xz[i], &rx, &ry, &rz, &rw, fs_xz);
        mapgrid(funs_xw[i], &rx, &ry, &rz, &rw, fs_xw);
        mapgrid(funs_yz[i], &rx, &ry, &rz, &rw, fs_yz);
        mapgrid(funs_yw[i], &rx, &ry, &rz, &rw, fs_yw);
        mapgrid(funs_zw[i], &rx, &ry, &rz, &rw, fs_zw);
        mapgrid(funs_xyz[i], &rx, &ry, &rz, &rw, fs_xyz);
        mapgrid(funs_xyw[i], &rx, &ry, &rz, &rw, fs_xyw);
        mapgrid(funs_xzw[i], &rx, &ry, &rz, &rw, fs_xzw);
        mapgrid(funs_yzw[i], &rx, &ry, &rz, &rw, fs_yzw);
        mapgrid(funs_xyzw[i], &rx, &ry, &rz, &rw, fs_xyzw);


        /* Create interpolation */
        Interpolation ip    = interpolation(&rx, &ry, &rz, &rw, fs, 
                                  fs_x, fs_y, fs_z, fs_w, fs_xy, fs_xz, fs_xw, 
                                  fs_yz, fs_yw, fs_zw, fs_xyz, fs_xyw, fs_xzw, 
                                  fs_yzw, fs_xyzw);

        /* Compare interpolated value to exact value */
        double err  = fabs(interpolate(&ip, xt, yt, zt, wt) - ft);
        assert((err < tols[1] + tols[1]*fabs(ft)) &&
               "interpolation(xt, yt, zt, wt)"); 

        /* Compare interpolated derivative to exact value */
        err  = fabs(diff_x(&ip, xt, yt, zt, wt) - ft_x);
        assert((err < tols[2] + tols[2]*fabs(ft_x)) &&
               "diff_x(xt, yt, zt, wt)"); 

        err  = fabs(diff_y(&ip, xt, yt, zt, wt) - ft_y);
        assert((err < tols[2] + tols[2]*fabs(ft_y)) &&
               "diff_y(xt, yt, zt, wt)"); 

        err  = fabs(diff_z(&ip, xt, yt, zt, wt) - ft_z);
        assert((err < tols[2] + tols[2]*fabs(ft_z)) &&
               "diff_z(xt, yt, zt, wt)"); 

        err  = fabs(diff_w(&ip, xt, yt, zt, wt) - ft_w);
        assert((err < tols[2] + tols[2]*fabs(ft_w)) &&
               "diff_w(xt, yt, zt, wt)"); 

        /**/

        /* Map exact values on grid points (without last one,
         * as it is out of bounds) */
        mapgrid(funs[i], &r0x, &r0y, &r0z, &r0w, fs0);
        mapgrid(funs_x[i], &r0x, &r0y, &r0z, &r0w, fs0_x);
        mapgrid(funs_y[i], &r0x, &r0y, &r0z, &r0w, fs0_y);
        mapgrid(funs_z[i], &r0x, &r0y, &r0z, &r0w, fs0_z);
        mapgrid(funs_w[i], &r0x, &r0y, &r0z, &r0w, fs0_w);

        double nfs  = norm(lt, fs0);
        double nfs_x= norm(lt, fs0_x);
        double nfs_y= norm(lt, fs0_y);
        double nfs_z= norm(lt, fs0_z);
        double nfs_w= norm(lt, fs0_w);

        /* Compare interpolated values to exact values on grid-points */
        map(&ip, lt, xs0, ys0, zs0, ws0, res);
        err  = errnorm(lt, fs0, res);
        assert((err < tols[0] + tols[0]*nfs) &&
               "map(interpolation, xs0, ys0, zs0, ws0)"); 

        /* Compare interpolated derivatives to exact values on grid-points */
        map_x(&ip, lt, xs0, ys0, zs0, ws0, res);
        err  = errnorm(lt, fs0_x, res);
        assert((err < tols[0] + tols[0]*nfs_x) &&
               "map(diff_x, xs0, ys0, zs0)"); 

        map_y(&ip, lt, xs0, ys0, zs0, ws0, res);
        err  = errnorm(lt, fs0_y, res);
        assert((err < tols[0] + tols[0]*nfs_y) &&
               "map(diff_y, xs0, ys0, zs0)"); 

        map_z(&ip, lt, xs0, ys0, zs0, ws0, res);
        err  = errnorm(lt, fs0_z, res);
        assert((err < tols[0] + tols[0]*nfs_z) &&
               "map(diff_z, xs0, ys0, zs0)"); 

        map_w(&ip, lt, xs0, ys0, zs0, ws0, res);
        err  = errnorm(lt, fs0_w, res);
        assert((err < tols[0] + tols[0]*nfs_w) &&
               "map(diff_w, xs0, ys0, zs0)"); 

        /**/

        /* Map exact values on test points */
        mapgrid(funs[i], &rtx, &rty, &rtz, &rtw, fs0);
        mapgrid(funs_x[i], &rtx, &rty, &rtz, &rtw, fs0_x);
        mapgrid(funs_y[i], &rtx, &rty, &rtz, &rtw, fs0_y);
        mapgrid(funs_z[i], &rtx, &rty, &rtz, &rtw, fs0_z);
        mapgrid(funs_w[i], &rtx, &rty, &rtz, &rtw, fs0_w);

        nfs     = norm(lt, fs0);
        nfs_x   = norm(lt, fs0_x);
        nfs_y   = norm(lt, fs0_y);
        nfs_z   = norm(lt, fs0_z);
        nfs_w   = norm(lt, fs0_w);

        /* Compare interpolated values to exact values on test-points */
        map(&ip, lt, xst, yst, zst, wst, res);
        err  = errnorm(lt, fs0, res);
        assert((err < tols[1] + tols[1]*nfs) &&
               "map(interpolation, xst, yst, zst)"); 

        /* Compare interpolated derivatives to exact values on test-points */
        map_x(&ip, lt, xst, yst, zst, wst, res);
        err  = errnorm(lt, fs0_x, res);
        assert((err < tols[1] + tols[1]*nfs_x) &&
               "map(diff_x, xst, yst, zst)"); 

        map_y(&ip, lt, xst, yst, zst, wst, res);
        err  = errnorm(lt, fs0_y, res);
        assert((err < tols[1] + tols[1]*nfs_y) &&
               "map(diff_y, xst, yst, zst)"); 

        map_z(&ip, lt, xst, yst, zst, wst, res);
        err  = errnorm(lt, fs0_z, res);
        assert((err < tols[1] + tols[1]*nfs_z) &&
               "map(diff_z, xst, yst, zst)"); 

        map_w(&ip, lt, xst, yst, zst, wst, res);
        err  = errnorm(lt, fs0_w, res);
        assert((err < tols[1] + tols[1]*nfs_w) &&
               "map(diff_w, xst, yst, zst)"); 

        loci_free(&ip);
    }

    free(fs);
    free(fs_x);
    free(fs_y);
    free(fs_z);
    free(fs_w);
    free(fs_xy);
    free(fs_xz);
    free(fs_xw);
    free(fs_yz);
    free(fs_yw);
    free(fs_zw);
    free(fs_xyz);
    free(fs_xyw);
    free(fs_xzw);
    free(fs_yzw);
    free(fs_xyzw);

    free(xs0);
    free(ys0);
    free(zs0);
    free(ws0);
    free(xst);
    free(yst);
    free(zst);
    free(wst);

    free(res);

    free(fs0);
    free(fs0_x);
    free(fs0_y);
    free(fs0_z);
    free(fs0_w);

    OUT(0, "Pass\n");
}


/* Test higher order derivative on exponantial function. */
void exp4D()  {
    OUT(0, "4D Exp Test");

    /* See basic4D() */
    Range rx    = {.x0 = 1, .dx = 0.3, .len = 20};
    Range ry    = {.x0 = -5, .dx = 0.2, .len = 30};
    Range rz    = {.x0 = 2, .dx = 0.1, .len = 10};
    Range rw    = {.x0 = -3, .dx = 0.4, .len = 15};

    Range rtx   = {rx.x0 + rx.dx/3, rx.dx, rx.len-1};
    Range rty   = {ry.x0 + ry.dx/3, ry.dx, ry.len-1};
    Range rtz   = {rz.x0 + rz.dx/3, rz.dx, rz.len-1};
    Range rtw   = {rw.x0 + rw.dx/3, rw.dx, rw.len-1};

    double xt   = rx.x0 + rx.dx*rx.len/2 + rx.dx/7;
    double yt   = ry.x0 + ry.dx*ry.len/2 + ry.dx/7;
    double zt   = rz.x0 + rz.dx*rz.len/2 + rz.dx/7;
    double wt   = rw.x0 + rw.dx*rw.len/2 + rw.dx/7;

    int lx      = rx.len;   
    int ly      = ry.len;   
    int lz      = rz.len;
    int lw      = rw.len;
    int le      = lx*ly*lz*lw;

    int ltx     = rtx.len;  
    int lty     = rty.len;  
    int ltz     = rtz.len;  
    int ltw     = rtw.len;  
    int lt      = ltx*lty*ltz*ltw;

    double *fs     = malloc(le*sizeof(double));
    double *fst    = malloc(lt*sizeof(double));
    double *xst    = malloc(lt*sizeof(double));
    double *yst    = malloc(lt*sizeof(double));
    double *zst    = malloc(lt*sizeof(double));
    double *wst    = malloc(lt*sizeof(double));
    double *res    = malloc(lt*sizeof(double));

    mapgrid(e4, &rx, &ry, &rz, &rw, fs);

    mapgrid(e4, &rtx, &rty, &rtz, &rtw, fst);
    double nfs  = norm(lt, fst);

    /* Create interpolation */
    Interpolation ip    = interpolation(&rx, &ry, &rz, &rw, fs, fs, fs, fs,
                                            fs, fs, fs, fs, fs, fs, fs, fs,
                                            fs, fs, fs, fs);

    rs2as(&rtx, &rty, &rtz, &rtw, xst, yst, zst, wst);

    double ft   = e4(xt, yt, zt, wt);


    for(int i = 0; i < 3; i++)  {
     for(int j = 0; j < 3; j++)  {
      for(int k = 0; k < 3; k++)  {
       for(int l = 0; l < 3; l++)  {
        OUT(1, "d^%if/dx^%i/dy^%i/dz^%i/dw^%i", i+j+k+l, i, j, k, l);
        double ti   = tols[MAX(i, MAX(j, MAX(k, l))) + 2];

        /* Compare interpolated derivative to exact value */
        double err  = fabs(diff(&ip, i, j, k, l, xt, yt, zt, wt) - ft);
        assert((err < ti + ti*fabs(ft)) &&
               "diff(i, j, k, l, xt, yt, zt, wt)"); 

        /* Compare interpolated derivatives to exact values on test-points */
        map_diff(&ip, i, j, k, l, lt, xst, yst, zst, wst, res);
        err  = errnorm(lt, fst, res);

        assert((err < ti + ti*nfs) &&
           "map(diff, i, j, k, l, xst, yst, zst, wst)"); 
       }
      }
     }
    }

    free(fs);
    free(fst);
    free(xst);
    free(yst);
    free(zst);
    free(wst);
    free(res);

    loci_free(&ip);

    OUT(0, "Pass\n");
}

/* Test out of bounds behaviour. */
void boundary4D()  {
    OUT(0, "4D Boundary Test");

    double big  = 1e6;
    double fr   = 89./97;

    /* See basic4D() */
    Range rx    = {.x0 = 1, .dx = 0.3, .len = 20};
    Range ry    = {.x0 = -5, .dx = 0.2, .len = 30};
    Range rz    = {.x0 = 2, .dx = 0.1, .len = 10};
    Range rw    = {.x0 = -3, .dx = 0.4, .len = 15};
    int lx      = rx.len;   
    int ly      = ry.len;   
    int lz      = rz.len;
    int lw      = rw.len;
    int le      = lx*ly*lz*lw;

    double xe   = rx.x0 + lx*rx.dx;
    double ye   = ry.x0 + ly*ry.dx;
    double ze   = rz.x0 + lz*rz.dx;
    double we   = rw.x0 + lw*rw.dx;

    double *fs      = malloc(le*sizeof(double));
    double *fs_x    = malloc(le*sizeof(double));
    double *fs_y    = malloc(le*sizeof(double));
    double *fs_z    = malloc(le*sizeof(double));
    double *fs_w    = malloc(le*sizeof(double));
    double *fs_xy   = malloc(le*sizeof(double));
    double *fs_xz   = malloc(le*sizeof(double));
    double *fs_xw   = malloc(le*sizeof(double));
    double *fs_yz   = malloc(le*sizeof(double));
    double *fs_yw   = malloc(le*sizeof(double));
    double *fs_zw   = malloc(le*sizeof(double));
    double *fs_xyz  = malloc(le*sizeof(double));
    double *fs_xyw  = malloc(le*sizeof(double));
    double *fs_xzw  = malloc(le*sizeof(double));
    double *fs_yzw  = malloc(le*sizeof(double));
    double *fs_xyzw = malloc(le*sizeof(double));

    mapgrid(f4, &rx, &ry, &rz, &rw, fs);
    mapgrid(f4_x, &rx, &ry, &rz, &rw, fs_x);
    mapgrid(f4_y, &rx, &ry, &rz, &rw, fs_y);
    mapgrid(f4_z, &rx, &ry, &rz, &rw, fs_z);
    mapgrid(f4_w, &rx, &ry, &rz, &rw, fs_w);
    mapgrid(f4_xy, &rx, &ry, &rz, &rw, fs_xy);
    mapgrid(f4_xz, &rx, &ry, &rz, &rw, fs_xz);
    mapgrid(f4_xw, &rx, &ry, &rz, &rw, fs_xw);
    mapgrid(f4_yz, &rx, &ry, &rz, &rw, fs_yz);
    mapgrid(f4_yw, &rx, &ry, &rz, &rw, fs_yw);
    mapgrid(f4_zw, &rx, &ry, &rz, &rw, fs_zw);
    mapgrid(f4_xyz, &rx, &ry, &rz, &rw, fs_xyz);
    mapgrid(f4_xyw, &rx, &ry, &rz, &rw, fs_xyw);
    mapgrid(f4_xzw, &rx, &ry, &rz, &rw, fs_xzw);
    mapgrid(f4_yzw, &rx, &ry, &rz, &rw, fs_yzw);
    mapgrid(f4_xyzw, &rx, &ry, &rz, &rw, fs_xyzw);

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

    double wts[7]   = {rw.x0 + fr*rw.dx,    /* in bounds */
                       we    + 1./big,
                       rw.x0 - 1./big,
                       we    + fr*rw.dx,
                       rw.x0 - fr*rw.dx,
                       we    + big,
                       rw.x0 - big
    };

    /* Create interpolation */
   Interpolation ip    = interpolation(&rx, &ry, &rz, &rw, fs, 
                                  fs_x, fs_y, fs_z, fs_w, fs_xy, fs_xz, fs_xw, 
                                  fs_yz, fs_yw, fs_zw, fs_xyz, fs_xyw, fs_xzw, 
                                  fs_yzw, fs_xyzw);


    for(int i = 0; i < 7; i++)  {
     for(int j = 0; j < 7; j++)  {
      for(int k = 0; k < 7; k++)  {
       for(int l = 0; l < 7; l++)  {
        if(i+j+k+l > 0)    {
            double xi   = xts[i];
            double yj   = yts[j];
            double zk   = zts[k];
            double wl   = wts[l];

            /* Test whether out of bounds == nan */
            assert(isnan(interpolate(&ip, xi, yj, zk, wl)) &&
                   "interpolate(xi, yj, zk, wk) out of bounds"); 

            assert(isnan(diff_x(&ip, xi, xi, zk, wl)) &&
                   "diff_x(xi, yj, zk, wk) out of bounds"); 

            assert(isnan(diff_y(&ip, xi, yj, zk, wl)) &&
                   "diff_y(xi, yj, zk, wk) out of bounds"); 

            assert(isnan(diff_z(&ip, xi, yj, zk, wl)) &&
                   "diff_z(xi, yj, zk, wk) out of bounds"); 

            assert(isnan(diff_w(&ip, xi, yj, zk, wl)) &&
                   "diff_z(xi, yj, zk, wk) out of bounds"); 

            assert(isnan(diff(&ip, 2, 2, 2, 2, xi, yj, zk, wl)) &&
                   "diff(2, 2, 2, 2, xi, yj, zk, wk) out of bounds"); 
        }
       }
      }
     }
    }

    free(fs);
    free(fs_x);
    free(fs_y);
    free(fs_z);
    free(fs_w);
    free(fs_xy);
    free(fs_xz);
    free(fs_xw);
    free(fs_yw);
    free(fs_yz);
    free(fs_zw);
    free(fs_xyz);
    free(fs_xyw);
    free(fs_xzw);
    free(fs_yzw);
    free(fs_xyzw);

    loci_free(&ip);

    OUT(0, "Pass\n");
}
