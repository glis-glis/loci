/* loci
 * Local cubic interpolations in up to 4 dimensions
 *
 * Copyright (C) 2017 Andreas Füglistaler <andreas.fueglistaler@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef FILE_3D_H
#define FILE_3D_H

#include "loci.h"

/* Create interpolation in ranges `rx`, `ry`, `rz` for values 
 * `fs`= f(`rx`, `ry`, `rz`) and derivative values `f_x` = df/dx, `f_y`,
 * `f_z`, `f_xy`, `f_xz`, `f_yz`, and `f_xyz` = d^3f/(dx dy dz).
 *
 * `fs` is of length rx->len*ry->len. Its values are on every gridpoint:
 *  for(int i = 0; i < rx->len;i++) {
 *   for(int j = 0; j < ry->len; j++) {
 *    for(int k = 0; k < rz->len; k++) {
 *      int ijk   = i*ry->len*rz->len + j*rz->len + k;
 *      real xi   = rx->x0 + rx->dx*i;
 *      real yj   = ry->x0 + ry->dx*j;
 *      real zk   = rz->x0 + rz->dx*k;
 *      fs[ijk]   = f(xi, yj, zk);
 *    }
 *   }
 *  } */
Interpolation interpolation3D(const Range *rx, const Range *ry,
                              const Range *rz, const real fs[],
                              const real fs_x[], const real fs_y[],
                              const real fs_z[], const real fs_xy[],
                              const real fs_xz[], const real fs_yz[],
                              const real fs_xyz[]);

/* Interpolate with interpolation `ip` at point (`x`, `y`, `z`). */
real interpolate3D(const Interpolation *ip, real x, real y, real z);

/* Interpolate first order derivative in x with interpolation `ip` 
 * at point (`x`, `y`, `z`). */
real diff3D_x(const Interpolation *ip, real x, real y, real z);

/* Interpolate first order derivative in y with interpolation `ip` 
 * at point (`x`, `y`, `z`). */
real diff3D_y(const Interpolation *ip, real x, real y, real z);

/* Interpolate first order derivative in z with interpolation `ip` 
 * at point (`x`, `y`, `z`). */
real diff3D_z(const Interpolation *ip, real x, real y, real z);

/* Interpolate (`edx`, `edy`, `edz`)-order derivative in (x, y, z) with 
 * interpolation `ip` at point (`x`, `y`, `z`). */
real diff3D(const Interpolation *ip, int edx, int edy, int edz,
            real x, real y, real z);
              

/* Map interpolation `ip` onto the input-arrays `xs`, `ys`, `zs` with length 
 * `len` and write result to output-array `out`. 
 *
 * The mapping is linear:
 * out[i] = interpolate(&ip, xs[i], ys[i], zs[i]) .*/
void map3D(const Interpolation *ip, int len, const real xs[],
           const real ys[], const real zs[], real out[]);

/* Map first order derivative in x of interpolation `ip` onto input-array `xs`, 
 * `ys`, `zs` with length `len` and write result to output-array `out`. */
void map3D_x(const Interpolation *ip, int len, const real xs[],
           const real ys[], const real zs[], real out[]);

/* Map first order derivative in y of interpolation `ip` onto input-array `xs`, 
 * `ys`, `zs` with length `len` and write result to output-array `out`. */
void map3D_y(const Interpolation *ip, int len, const real xs[],
           const real ys[], const real zs[], real out[]);

/* Map first order derivative in z of interpolation `ip` onto input-array `xs`, 
 * `ys`, `zs` with length `len` and write result to output-array `out`. */
void map3D_z(const Interpolation *ip, int len, const real xs[],
           const real ys[], const real zs[], real out[]);

/* Map (`edx`, `edy`, `edz`)-order derivative in (x, y, z) of interpolation `ip`
 * onto input-array `xs`, `ys`, `zs` with length `len` and write result to 
 * output-array `out`. */
void map3D_diff(const Interpolation *ip, int edx, int edy, int edz, int len,
                const real xs[], const real ys[], const real zs[], real out[]);

#endif
