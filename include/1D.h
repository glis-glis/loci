/* loci
 * Local cubic interpolations in up to 4 dimensions
 *
 * Copyright (C) 2017 Andreas Füglistaler <andreas.fueglistaler@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef FILE_1D_H
#define FILE_1D_H

#include "loci.h"

/* Create interpolation in range `rx` for values `fs`= f(`rx`) and 
 * derivative values `f_x` = df/dx(`rx`). 
 *
 * Given a range r with x0, dx and len
 * g(r) = {g(x0), g(x0 + dx), ..., g(x0+(len-1)*dx). */
Interpolation interpolation1D(const Range *rx, const real fs[],
                              const real fs_x[]);

/* Interpolate with interpolation `ip` at point `x`. */
real interpolate1D(const Interpolation *ip, real x);

/* Interpolate first order derivative in x with interpolation `ip` 
 * at point `x`. */
real diff1D_x(const Interpolation *ip, real x);

/* Interpolate `edx`-order derivative in x with interpolation `ip` 
 * at point `x`. */
real diff1D(const Interpolation *ip, int edx, real x);

/* Map interpolation `ip` onto input-array `xs` with length `len` and write 
 * result to output-array `out`*/
void map1D(const Interpolation *ip, int len, const real xs[], real out[]);

/* Map first order derivative in x of interpolation `ip` onto input-array `xs` 
 * with length `len` and write result to output-array `out`*/
void map1D_x(const Interpolation *ip, int len, const real xs[], real out[]);

/* Map `edx`-order derivative in x of interpolation `ip` onto input-array `xs` 
 * with length `len` and write result to output-array `out`*/
void map1D_diff(const Interpolation *ip, int edx, int len, 
                const real xs[], real out[]);

#endif
