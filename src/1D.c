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
#include "../include/1D.h"

static const int NKS    = 4;

/* Private functions */
static void cellcoeffs(const real fs[2], const real fs_x[2], real ks[NKS]);
static real poly(const real ks[NKS], int edx, real x);

/* "Constructor" */
Interpolation interpolation1D(const Range *rx, const real fs[], 
                              const real fs_x[]) {

    real *ds_x  = malloc(rx->len*sizeof(real));
    int le      = rx->len - 1;

    for(int i = 0; i < rx->len; i++)  {
        ds_x[i] = rx->dx*fs_x[i];
    }

    Interpolation ip;
    ip.x0s   = malloc(sizeof(real)); 
    ip.dxs   = malloc(sizeof(real)); 
    ip.lens  = malloc(sizeof(int)); 
    ip.ks    = malloc(NKS*le*sizeof(real));

    ip.x0s[0]   = rx->x0;
    ip.dxs[0]   = rx->dx;
    ip.lens[0]  = le;


    for(int i = 0; i < le; i++)  {
        cellcoeffs(&fs[i], &ds_x[i], &ip.ks[i*4]);
    }

    free(ds_x);

    return ip;
}

real interpolate1D(const Interpolation* ip, real x)    {
    int i;
    real xp;

    xp_i(&xp, &i, x, ip->x0s[0], ip->dxs[0]);

    if(i>= 0 && i < ip->lens[0])   {
        return poly(&ip->ks[NKS*i], 0, xp);
    }
    else    {
        return NAN;
    }
}
real diff1D_x(const Interpolation* ip, real x)    {
    int i;
    real xp;

    xp_i(&xp, &i, x, ip->x0s[0], ip->dxs[0]);

    if(i>= 0 && i < ip->lens[0])   {
        return poly(&ip->ks[NKS*i], 1, xp)/ip->dxs[0];
    }
    else    {
        return NAN;
    }
}

real diff1D(const Interpolation* ip, int edx, real x)    {
    if(edx > 3)  {
        return 0.;
    }
    else if(edx < 0)  {
        return NAN;
    }
    else {
        int i;
        real xp;

        xp_i(&xp, &i, x, ip->x0s[0], ip->dxs[0]);

        if(i>= 0 && i < ip->lens[0])   {
            return poly(&ip->ks[NKS*i], edx, xp)/ipow(ip->dxs[0], edx);
        }
        else    {
            return NAN;
        }
    }
}

void map1D(const Interpolation* ip, int len, 
           const real xs[], real out[])  {
    for(int i = 0; i < len; i++) {
        out[i]  = interpolate1D(ip, xs[i]);
    }
}

void map1D_x(const Interpolation* ip, int len, 
           const real xs[], real out[])  {
    for(int i = 0; i < len; i++) {
        out[i]  = diff1D_x(ip, xs[i]);
    }
}

void map1D_diff(const Interpolation* ip, int edx, int len, 
           const real xs[], real out[])  {
    for(int i = 0; i < len; i++) {
        out[i]  = diff1D(ip, edx, xs[i]);
    }
}

/* Private functions */

/* Create polynome.
 * Any modern compiler should unroll the loop. */
real poly(const real ks[NKS], int edx, real x)    {
    real r    = 0;
    for(int i = 0; i < 4; i++)  {
        r   += ks[i]*xi(x, edx, i); 
    }
    return r;
}


void cellcoeffs(const real fs[2], const real fs_x[2], real ks[NKS]) {
    ks[0] = fs[0];
    ks[1] = fs_x[0];
    ks[2] = -2*fs_x[0] - 3*fs[0] + 3*fs[1] - fs_x[1];
    ks[3] = 2*fs[0] + fs_x[0] - 2*fs[1] + fs_x[1];
}
