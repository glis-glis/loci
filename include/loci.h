/* loci
 * Local multivariate cubic interpolations in up to 4 dimensions
 *
 * Copyright (C) 2017 Andreas Füglistaler <andreas.fueglistaler@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef FILE_LOCI_H
#define FILE_LOCI_H

/* If user defines precision */
#ifdef REAL
typedef REAL real;
/* Else use double precision */
#else
typedef double real;
#endif

/* Return sizeof(real) */
int loci_size(void);

/* Range from `x0` to `x0`+(`len`-1)*`dx`, with step `dx`and length `len`. */
typedef struct {
    real x0;
    real dx;
    int  len;
} Range;

/* Interpolation in up to 4 dimensions dim on a grid with interpolation 
 * coeffisients ks. The grid parameters `x0s`[dim], `dxs`[dim] and `lens`[dim]
 * have the same definition as in the structure `Range`. */
typedef struct {
    real *x0s;
    real *dxs;
    int  *lens;
    real *ks;
} Interpolation;

/* Destructor of interpolation `ip`, used after last use of a Interpolation 
 * object. */
void loci_free(Interpolation* ip);

#include "1D.h"
#include "2D.h"
#include "3D.h"
#include "4D.h"
#include "macros.h"

/* Compile-time multiple dispatch macros. */

/* interpolation(...) with NARGS = number of arguments
 *      NARGS= 3: interpolation1D(ARGS)
 *      NARGS= 6: interpolation2D(ARGS)
 *      NARGS=11: interpolation3D(ARGS)
 *      NARGS=20: interpolation4D(ARGS). */
#define interpolation(...)  \
        _CONC(_INTERPOLATION, _NARGS(__VA_ARGS__))(__VA_ARGS__)

/* interpolate(...), with NARGS = number of arguments
 *      NARGS=2: interpolate1D(ARGS)
 *      NARGS=3: interpolate2D(ARGS)
 *      NARGS=4: interpolate3D(ARGS)
 *      NARGS=5: interpolate4D(ARGS)
 *
 *      Same for diff_x, diff_y, diff_z, diff_w. */
#define interpolate(...) \
        _CONC(_INTERPOLATE, _NARGS(__VA_ARGS__))(__VA_ARGS__)

#define diff_x(...) \
        _CONC(_DIFF_X, _NARGS(__VA_ARGS__))(__VA_ARGS__)

#define diff_y(...) \
        _CONC(_DIFF_Y, _NARGS(__VA_ARGS__))(__VA_ARGS__)

#define diff_z(...) \
        _CONC(_DIFF_Z, _NARGS(__VA_ARGS__))(__VA_ARGS__)

#define diff_w diff4D_w

/* diff(...), with NARGS = number of arguments
 *      NARGS=3: diff1D(ARGS)
 *      NARGS=5: diff2D(ARGS)
 *      NARGS=7: diff3D(ARGS)
 *      NARGS=9: diff4D(ARGS). */
#define diff(...) \
        _CONC(_DIFF, _NARGS(__VA_ARGS__))(__VA_ARGS__)

/* map(...), with NARGS = number of arguments
 *      NARGS=4: map1D(ARGS)
 *      NARGS=5: map2D(ARGS)
 *      NARGS=6: map3D(ARGS)
 *      NARGS=7: map4D(ARGS) 
 *
 *      Same for map_x, map_y, map_z, map_w. */
#define map(...) \
        _CONC(_MAP, _NARGS(__VA_ARGS__))(__VA_ARGS__)

#define map_x(...) \
        _CONC(_MAP_X, _NARGS(__VA_ARGS__))(__VA_ARGS__)

#define map_y(...) \
        _CONC(_MAP_Y, _NARGS(__VA_ARGS__))(__VA_ARGS__)

#define map_z(...) \
        _CONC(_MAP_Z, _NARGS(__VA_ARGS__))(__VA_ARGS__)

#define map_w map4D_w

/* map_diff(...), with NARGS = number of arguments
 *      NARGS= 5: map_diff1D(ARGS)
 *      NARGS= 7: map_diff2D(ARGS)
 *      NARGS= 9: map_diff3D(ARGS)
 *      NARGS=11: map_diff4D(ARGS). */
#define map_diff(...) \
        _CONC(_MAP_DIFF, _NARGS(__VA_ARGS__))(__VA_ARGS__)

#endif
