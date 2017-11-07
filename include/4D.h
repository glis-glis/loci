/* loci
 * Local cubic interpolations in up to 4 dimensions
 *
 * Copyright (C) 2017 Andreas Füglistaler <andreas.fueglistaler@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef FILE_4D_H
#define FILE_4D_H

#include "loci.h"

/* Create interpolation in ranges `rx`, `ry`, `rz`, `rw` for values 
 * `fs`= f(`rx`, `ry`, `rz`, `rw`) and derivative values `f_x` = df/dx, `f_y`,
 * `f_z`, `f_w`, `f_xy`, `f_xz`, `f_xw`, `f_yz`, `f_yw`, `f_zw`, `f_xyz`,
 * `f_xyz`, `f_xyw`,`f_xzw`, `f_yzw`, `f_xyzw`= d^4f/(dx dy dz dw).
 *
 * `fs` is of length rx->len*ry->len. Its values are on every gridpoint:
 * for(int i = 0; i < rx->len;i++) {
 *  for(int j = 0; j < ry->len; j++) {
 *   for(int k = 0; k < rz->len; k++) {
 *    for(int l = 0; l < rw->len; l++) {
 *      int ijkl  = i*ry->len*rz->len*rw->len + j*rz->len*rw->len + 
 *                  k*rw->len + l;
 *      real xi   = rx->x0 + rx->dx*i;
 *      real yj   = ry->x0 + ry->dx*j;
 *      real zk   = rz->x0 + rz->dx*k;
 *      real wl   = rw->x0 + rw->dx*l;
 *      fs[ijkl]  = f(xi, yj, zk, wl);
 *    }
 *   }
 *  }
 * }
*/
Interpolation interpolation4D(const Range *rx, const Range *ry, 
                              const Range *rz, const Range *rw, 
                              const real fs[], const real fs_x[], 
                              const real fs_y[], const real fs_z[],
                              const real fs_w[], const real fs_xy[], 
                              const real fs_xz[], const real fs_xw[], 
                              const real fs_yz[], const real fs_yw[], 
                              const real fs_zw[], const real fs_xyz[], 
                              const real fs_xyw[], const real fs_xzw[], 
                              const real fs_yzw[], const real fs_xyzw[]);

/* Interpolate with interpolation `ip` at point (`x`, `y`, `z`, `w`). */
real interpolate4D(const Interpolation *ip, real x, real y, real z, real w);

/* Interpolate first order derivative in x with interpolation `ip` 
 * at point (`x`, `y`, `z`, `w`). */
real diff4D_x(const Interpolation *ip, real x, real y, real z, real w);

/* Interpolate first order derivative in y with interpolation `ip` 
 * at point (`x`, `y`, `z`, `w`). */
real diff4D_y(const Interpolation *ip, real x, real y, real z, real w);

/* Interpolate first order derivative in z with interpolation `ip` 
 * at point (`x`, `y`, `z`, `w`). */
real diff4D_z(const Interpolation *ip, real x, real y, real z, real w);

/* Interpolate first order derivative in w with interpolation `ip` 
 * at point (`x`, `y`, `z`, `w`). */
real diff4D_w(const Interpolation *ip, real x, real y, real z, real w);

/* Interpolate (`edx`, `edy`, `edz`, `edw`)-order derivative in (x, y, z, w) 
 * with interpolation `ip` at point (`x`, `y`, `z`, `w`). */
real diff4D(const Interpolation *ip, int edx, int edy, int edz, int edw,
            real x, real y, real z, real w);

/* Map interpolation `ip` onto the input-arrays `xs`, `ys`, `zs`, `ws` with 
 * length `len` and write result to output-array `out`. 
 *
 * The mapping is linear:
 * out[i] = interpolate(&ip, xs[i], ys[i], zs[i], ws[i]) .*/
void map4D(const Interpolation *ip, int len, const real xs[], const real ys[], 
           const real zs[], const real ws[], real out[]);

/* Map first order derivative in x of interpolation `ip` onto input-array `xs`, 
 * `ys`, `zs`, `ws` with length `len` and write result to output-array `out`. */
void map4D_x(const Interpolation *ip, int len, const real xs[], 
             const real ys[], const real zs[], const real ws[], real out[]);

/* Map first order derivative in y of interpolation `ip` onto input-array `xs`, 
 * `ys`, `zs`, `ws` with length `len` and write result to output-array `out`. */
void map4D_y(const Interpolation *ip, int len, const real xs[],
           const real ys[], const real zs[], const real ws[], real out[]);

/* Map first order derivative in z of interpolation `ip` onto input-array `xs`, 
 * `ys`, `zs`, `ws` with length `len` and write result to output-array `out`. */
void map4D_z(const Interpolation *ip, int len, const real xs[],
           const real ys[], const real zs[], const real ws[], real out[]);

/* Map first order derivative in w of interpolation `ip` onto input-array `xs`, 
 * `ys`, `zs`, `ws` with length `len` and write result to output-array `out`. */
void map4D_w(const Interpolation *ip, int len, const real xs[],
           const real ys[], const real zs[], const real ws[], real out[]);

/* Map (`edx`, `edy`, `edz`, `edw`)-order derivative in (x, y, z, w) of 
 * interpolation `ip` onto input-array `xs`, `ys`, `zs`, `ws` with length `len`
 * and write result to output-array `out`. */
void map4D_diff(const Interpolation *ip, int edx, int edy, int edz, int edw,
                int len, const real xs[], const real ys[], 
                const real zs[], const real ws[], real out[]);

#endif
