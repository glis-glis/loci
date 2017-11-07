/* loci
 * Local cubic interpolations in up to 4 dimensions
 *
 * Copyright (C) 2017 Andreas Füglistaler <andreas.fueglistaler@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <stdlib.h>

#include "helpers.h"
#include "../include/loci.h"
#include "../include/2D.h"

static const int NKS    = 16;

/* Declare private functions */
static void cellcoeffs(int ly, const real fs[], const real fs_x[],
                      const real fs_y[], const real fs_xy[], real ks[NKS]);
static real poly(const real ks[NKS], int edx, int edy, real x, real y);

/* "Constructor" */
Interpolation interpolation2D(const Range *rx, const Range *ry, 
                              const real fs[], const real fs_x[], 
                              const real fs_y[], const real fs_xy[]) {

    int lx      = rx->len;
    int ly      = ry->len;
    real *ds_x  = malloc(lx*ly*sizeof(real));
    real *ds_y  = malloc(lx*ly*sizeof(real));
    real *ds_xy = malloc(lx*ly*sizeof(real));

    for(int i = 0; i < lx; i++)  {
        for(int j = 0; j < ly; j++)  {
            int ij    = i*ly + j;
            ds_x[ij]  = rx->dx*fs_x[ij];
            ds_y[ij]  = ry->dx*fs_y[ij];
            ds_xy[ij] = rx->dx*ry->dx*fs_xy[ij];
        }
    }

    Interpolation ip;
    ip.x0s  = malloc(2*sizeof(real)); 
    ip.dxs  = malloc(2*sizeof(real)); 
    ip.lens = malloc(2*sizeof(int)); 

    ip.x0s[0]   = rx->x0;
    ip.x0s[1]   = ry->x0;
    ip.dxs[0]   = rx->dx;
    ip.dxs[1]   = ry->dx;
    ip.lens[0]  = rx->len - 1;
    ip.lens[1]  = ry->len - 1;

    ip.ks   = malloc(NKS*ip.lens[0]*ip.lens[1]*sizeof(real));

    for(int i = 0; i < ip.lens[0]; i++)  {
        for(int j = 0; j < ip.lens[1]; j++)  {
            int ijk = NKS*(i*ip.lens[1] + j);
            int ijf = i*ly + j;
            cellcoeffs(ly, &fs[ijf], &ds_x[ijf], 
                       &ds_y[ijf], &ds_xy[ijf], &ip.ks[ijk]);
        }
    }

    free(ds_x);
    free(ds_y);
    free(ds_xy);

    return ip;
}

real interpolate2D(const Interpolation* ip, real x, real y)    {
    int i, j;
    real xp, yp;

    xp_i(&xp, &i, x, ip->x0s[0], ip->dxs[0]);
    xp_i(&yp, &j, y, ip->x0s[1], ip->dxs[1]);


    if(i>= 0 && i < ip->lens[0] && j>= 0 && j < ip->lens[1])   {
        int ij  = NKS*(i*ip->lens[1] + j);
        return poly(&ip->ks[ij], 0, 0, xp, yp);
    }
    else    {
        return NAN;
    }
}

real diff2D_x(const Interpolation* ip, real x, real y)    {
    int i, j;
    real xp, yp;

    xp_i(&xp, &i, x, ip->x0s[0], ip->dxs[0]);
    xp_i(&yp, &j, y, ip->x0s[1], ip->dxs[1]);

    if(i>= 0 && i < ip->lens[0] && j>= 0 && j < ip->lens[1])   {
        int ij  = NKS*(i*ip->lens[1] + j);
        return poly(&ip->ks[ij], 1, 0, xp, yp)/ip->dxs[0];
    }
    else    {
        return NAN;
    }
}

real diff2D_y(const Interpolation* ip, real x, real y)    {
    int i, j;
    real xp, yp;

    xp_i(&xp, &i, x, ip->x0s[0], ip->dxs[0]);
    xp_i(&yp, &j, y, ip->x0s[1], ip->dxs[1]);

    if(i>= 0 && i < ip->lens[0] && j>= 0 && j < ip->lens[1])   {
        int ij  = NKS*(i*ip->lens[1] + j);
        return poly(&ip->ks[ij], 0, 1, xp, yp)/ip->dxs[1];
    }
    else    {
        return NAN;
    }
}

real diff2D(const Interpolation* ip, int edx, int edy, real x, real y)    {
    int i, j;
    real xp, yp;

    xp_i(&xp, &i, x, ip->x0s[0], ip->dxs[0]);
    xp_i(&yp, &j, y, ip->x0s[1], ip->dxs[1]);

    if(i>= 0 && i < ip->lens[0] && j>= 0 && j < ip->lens[1])   {
        int ij  = NKS*(i*ip->lens[1] + j);
        return poly(&ip->ks[ij], edx, edy, xp, yp)/ipow(ip->dxs[0], edx)/
                   ipow(ip->dxs[1], edy);
    }
    else    {
        return NAN;
    }
}

void map2D(const Interpolation* ip, int len, 
           const real xs[], const real ys[], real out[])  {
    for(int i = 0; i < len; i++) {
        out[i]  = interpolate2D(ip, xs[i], ys[i]);
    }
}

void map2D_x(const Interpolation* ip, int len, 
           const real xs[], const real ys[], real out[])  {
    for(int i = 0; i < len; i++) {
        out[i]  = diff2D_x(ip, xs[i], ys[i]);
    }
}

void map2D_y(const Interpolation* ip, int len, 
           const real xs[], const real ys[], real out[])  {
    for(int i = 0; i < len; i++) {
        out[i]  = diff2D_y(ip, xs[i], ys[i]);
    }
}

void map2D_diff(const Interpolation *ip, int edx, int edy, int len,
                const real xs[], const real ys[], real out[]) {
    for(int i = 0; i < len; i++) {
        out[i]  = diff2D(ip, edx, edy, xs[i], ys[i]);
    }
}

/* Private functions */

/* Create polynome.
 * Any modern compiler should unroll the loop. */
real poly(const real ks[NKS], int edx, int edy, real x, real y) {
    real r    = 0;
    for(int i = 0; i < 4; i++)  {
        real xx = xi(x, edx, i);
        for(int j = 0; j < 4; j++)  {
            int ij  = i*4 + j;
            r   += ks[ij]*xx*xi(y, edy, j); 
        }
    }
    return r;
}

void cellcoeffs(int ly, const real fs[], const real fs_x[],
                const real fs_y[], const real fs_xy[], real ks[NKS]) {

    int i00 = 0;
    int i01 = 1;
    int i10 = ly;
    int i11 = ly + 1;

    /* Maple generated code, do not touch */

    real t1 = 2 * fs_y[i00];
    real t2 = 3 * fs[i00];
    real t5 = 2 * fs[i00];
    real t8 = 2 * fs_xy[i00];
    real t9 = 3 * fs_x[i00];
    real t10 = 3 * fs_x[i01];
    real t12 = 2 * fs_x[i00];
    real t13 = 2 * fs_x[i01];
    real t17 = 3 * fs_y[i00];
    real t18 = 3 * fs_y[i10];
    real t26 = 3 * fs_x[i10];
    real t27 = 3 * fs_x[i11];
    real t29 = 2 * fs_xy[i01];
    real t30 = 2 * fs_xy[i10];
    real t32 = 3 * fs_y[i01];
    real t34 = 3 * fs_y[i11];
    real t35 = 9 * fs[i00] - 9 * fs[i01] - 9 * fs[i10] + 9 * fs[i11] + 6 * fs_x[i00] - 6 * fs_x[i01] + t26 - t27 + 4 * fs_xy[i00] + t29 + t30 + fs_xy[i11] + 6 * fs_y[i00] + t32 - 6 * fs_y[i10] - t34;
    real t38 = 6 * fs[i00];
    real t39 = 6 * fs[i01];
    real t40 = 6 * fs[i10];
    real t41 = 6 * fs[i11];
    real t42 = 2 * fs_x[i11];
    real t43 = 2 * fs_x[i10];
    real t44 = -t8 - 4 * fs_x[i00] + 4 * fs_x[i01] - t38 + t39 + t40 - t17 + t18 - t41 - fs_xy[i10] + t42 - t43 + t34 - t29 - t32 - fs_xy[i11];
    real t47 = 2 * fs_y[i10];
    real t51 = 2 * fs_y[i11];
    real t52 = 2 * fs_y[i01];
    real t53 = -4 * fs_y[i00] - t8 + 4 * fs_y[i10] - t38 + t39 - t9 + t10 + t40 - t41 - t30 + t27 - t26 + t51 - fs_xy[i01] - t52 - fs_xy[i11];
    real t58 = 4 * fs[i00] - 4 * fs[i01] + t12 + fs_xy[i00] - t13 - 4 * fs[i10] + t1 - t47 + 4 * fs[i11] + fs_xy[i10] - t42 + t43 - t51 + fs_xy[i01] + t52 + fs_xy[i11];

    ks[0] = fs[i00];
    ks[1] = fs_y[i00];
    ks[2] = -t1 - t2 + 3 * fs[i01] - fs_y[i01];
    ks[3] = t5 + fs_y[i00] - 2 * fs[i01] + fs_y[i01];
    ks[4] = fs_x[i00];
    ks[5] = fs_xy[i00];
    ks[6] = -t8 - t9 + t10 - fs_xy[i01];
    ks[7] = t12 + fs_xy[i00] - t13 + fs_xy[i01];
    ks[8] = -t12 - t2 + 3 * fs[i10] - fs_x[i10];
    ks[9] = -t8 - t17 + t18 - fs_xy[i10];
    ks[10] = t35;
    ks[11] = t44;
    ks[12] = t5 + fs_x[i00] - 2 * fs[i10] + fs_x[i10];
    ks[13] = t1 + fs_xy[i00] - t47 + fs_xy[i10];
    ks[14] = t53;
    ks[15] = t58;
}
