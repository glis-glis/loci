/* loci
 * Local cubic interpolations in up to 4 dimensions
 *
 * Copyright (C) 2017 Andreas Füglistaler <andreas.fueglistaler@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef FILE_2D_H
#define FILE_2D_H

#include "loci.h"

/* Create interpolation in ranges `rx`, `ry` for values 
 * `fs`= f(`rx`, `ry`) and derivative values `f_x` = df/dx, `f_y` and 
 * `f_xy` = d^2f/(dx dy). 
 *
 * `fs` is of length rx->len*ry->len. Its values are on every gridpoint:
 * for(int i = 0; i < rx->len; i++) {
 *      for(int j = 0; j < ry->len; j++) {
 *          int ij  = i*ry->len + j;
 *          real xi = rx->x0 + rx->dx*i;
 *          real yj = ry->x0 + ry->dx*j;
 *          fs[ij]  = f(xi, yj);
 *      }
 *  } */
Interpolation interpolation2D(const Range *rx, const Range *ry,
                              const real fs[], const real fs_x[],
                              const real fs_y[], const real fs_xy[]);

/* Interpolate with interpolation `ip` at point (`x`, `y`). */
real interpolate2D(const Interpolation *ip, real x, real y);

/* Interpolate first order derivative in x with interpolation `ip` 
 * at point (`x`, `y`). */
real diff2D_x(const Interpolation *ip, real x, real y);

/* Interpolate first order derivative in y with interpolation `ip` 
 * at point (`x`, `y`). */
real diff2D_y(const Interpolation *ip, real x, real y);

/* Interpolate (`edx`, `edy`)-order derivative in (x, y) with interpolation 
 * `ip` at point (`x`, `y`). */
real diff2D(const Interpolation *ip, int edx, int edy, real x, real y);

/* Map interpolation `ip` onto the input-arrays `xs`, `ys` with length `len` 
 * and write result to output-array `out`. 
 *
 * The mapping is linear:
 * out[i] = interpolate(&ip, xs[i], ys[i]) .*/
void map2D(const Interpolation *ip, int len, const real xs[], 
           const real ys[], real out[]);

/* Map first order derivative in x of interpolation `ip` onto input-array `xs`, 
 * `ys` with length `len` and write result to output-array `out`. */
void map2D_x(const Interpolation *ip, int len, const real xs[], 
             const real ys[], real out[]);

/* Map first order derivative in y of interpolation `ip` onto input-array `xs`, 
 * `ys` with length `len` and write result to output-array `out`. */
void map2D_y(const Interpolation *ip, int len, const real xs[], 
             const real ys[], real out[]);

/* Map (`edx`, `edy`)-order derivative in (x, y) of interpolation `ip` onto 
 * input-array `xs`, `ys` with length `len` and write result to 
 * output-array `out`. */
void map2D_diff(const Interpolation *ip, int edx, int edy, int len,
                const real xs[], const real ys[], real out[]);

#endif
