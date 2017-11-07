/* Copyright (C) 2017 Andreas Füglistaler <andreas.fueglistaler@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef FILE_HELPERS_H
#define FILE_HELPERS_H

#include <math.h>
#include <stdio.h>

#include "../include/loci.h"

static void xp_i(real *xp, int *i, real x, real x0, real dx)  {
    real xr = x - x0;  
    real xt = xr/dx;  
    real id = floor(xt); 

    *xp     = xt - id;           
    *i      = (int)id;          
}

static real ipow(real x, int n)    {
    int i;
    real r = 1;
    for(i = 0; i < n; i++)  {
        r *= x;
    }
    return r;
}

static real xi(real x, int edx, int i)    {
    static const int as[4][4]   =  {
        {1, 1, 1, 1},
        {0, 1, 2, 3},
        {0, 0, 2, 6},
        {0, 0, 0, 6}
    };
    return as[edx][i]*ipow(x, i-edx);
}

#endif
