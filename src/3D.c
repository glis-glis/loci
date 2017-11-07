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
#include "../include/3D.h"

static const int NKS    = 64;

/* Declare private functions */
static void cellcoeffs(int ly, int lz, const real fs[], const real fs_x[],
                       const real fs_y[], const real fs_z[], 
                       const real fs_xy[], const real fs_xz[], 
                       const real fs_yz[], const real fs_xyz[], real ks[NKS]);

static real poly(const real ks[NKS], int edx, int edy, int edz, 
                 real x, real y, real z);


/* "Constructor" */
Interpolation interpolation3D(const Range *rx, const Range *ry, 
                              const Range *rz, const real fs[], 
                              const real fs_x[], const real fs_y[], 
                              const real fs_z[], const real fs_xy[],
                              const real fs_xz[], const real fs_yz[],
                              const real fs_xyz[])   {

    int lx          = rx->len;
    int ly          = ry->len;
    int lz          = rz->len;
    int l           = lx*ly*lz;

    real *ds_x      = malloc(l*sizeof(real));
    real *ds_y      = malloc(l*sizeof(real));
    real *ds_z      = malloc(l*sizeof(real));
    real *ds_xy     = malloc(l*sizeof(real));
    real *ds_xz     = malloc(l*sizeof(real));
    real *ds_yz     = malloc(l*sizeof(real));
    real *ds_xyz    = malloc(l*sizeof(real));

    for(int i = 0; i < lx; i++)  {
     for(int j = 0; j < ly; j++)  {
      for(int k = 0; k < lz; k++)  {
        int ijk     = i*ly*lz + j*lz + k;
        ds_x[ijk]   = rx->dx*fs_x[ijk];
        ds_y[ijk]   = ry->dx*fs_y[ijk];
        ds_z[ijk]   = rz->dx*fs_z[ijk];

        ds_xy[ijk]  = rx->dx*ry->dx*fs_xy[ijk];
        ds_xz[ijk]  = rx->dx*rz->dx*fs_xz[ijk];
        ds_yz[ijk]  = ry->dx*rz->dx*fs_yz[ijk];
        ds_xyz[ijk] = rx->dx*ry->dx*rz->dx*fs_xyz[ijk];
      }
     }
    }

    Interpolation ip;
    ip.x0s  = malloc(3*sizeof(real)); 
    ip.dxs  = malloc(3*sizeof(real)); 
    ip.lens = malloc(3*sizeof(int)); 

    ip.x0s[0]   = rx->x0;
    ip.x0s[1]   = ry->x0;
    ip.x0s[2]   = rz->x0;
    ip.dxs[0]   = rx->dx;
    ip.dxs[1]   = ry->dx;
    ip.dxs[2]   = rz->dx;
    ip.lens[0]  = rx->len - 1;
    ip.lens[1]  = ry->len - 1;
    ip.lens[2]  = rz->len - 1;

    ip.ks  = malloc(NKS*ip.lens[0]*ip.lens[1]*ip.lens[2]*sizeof(real));

    for(int i = 0; i < ip.lens[0]; i++)  {
     for(int j = 0; j < ip.lens[1]; j++)  {
      for(int k = 0; k < ip.lens[2]; k++)  {
        int ijkk= NKS*(i*ip.lens[1]*ip.lens[2] + j*ip.lens[2] + k);
        int ijkf= i*ly*lz + j*lz + k;
        cellcoeffs(ly, lz, &fs[ijkf], &ds_x[ijkf], &ds_y[ijkf], &ds_z[ijkf], 
                   &ds_xy[ijkf], &ds_xz[ijkf], &ds_yz[ijkf], &ds_xyz[ijkf], 
                   &ip.ks[ijkk]);
      }
     }
    }

    free(ds_x);
    free(ds_y);
    free(ds_z);
    free(ds_xy);
    free(ds_xz);
    free(ds_yz);
    free(ds_xyz);

    return ip;
}

real interpolate3D(const Interpolation* ip, real x, real y, real z) {
    int i, j, k;
    real xp, yp, zp;

    xp_i(&xp, &i, x, ip->x0s[0], ip->dxs[0]);
    xp_i(&yp, &j, y, ip->x0s[1], ip->dxs[1]);
    xp_i(&zp, &k, z, ip->x0s[2], ip->dxs[2]);

    if(i>= 0 && i < ip->lens[0] && j>= 0 && j < ip->lens[1] && 
        k>= 0 && k < ip->lens[2])   {
        int ijk = NKS*(i*ip->lens[1]*ip->lens[2] + j*ip->lens[2] + k);
        return poly(&ip->ks[ijk], 0, 0, 0, xp, yp, zp);
    }
    else    {
        return NAN;
    }
}

real diff3D_x(const Interpolation* ip, real x, real y, real z) {
    int i, j, k;
    real xp, yp, zp;

    xp_i(&xp, &i, x, ip->x0s[0], ip->dxs[0]);
    xp_i(&yp, &j, y, ip->x0s[1], ip->dxs[1]);
    xp_i(&zp, &k, z, ip->x0s[2], ip->dxs[2]);

    if(i>= 0 && i < ip->lens[0] && j>= 0 && j < ip->lens[1] && 
        k>= 0 && k < ip->lens[2])   {
        int ijk = NKS*(i*ip->lens[1]*ip->lens[2] + j*ip->lens[2] + k);
        return poly(&ip->ks[ijk], 1, 0, 0, xp, yp, zp)/ip->dxs[0];
    }
    else    {
        return NAN;
    }
}

real diff3D_y(const Interpolation* ip, real x, real y, real z) {
    int i, j, k;
    real xp, yp, zp;

    xp_i(&xp, &i, x, ip->x0s[0], ip->dxs[0]);
    xp_i(&yp, &j, y, ip->x0s[1], ip->dxs[1]);
    xp_i(&zp, &k, z, ip->x0s[2], ip->dxs[2]);

    if(i>= 0 && i < ip->lens[0] && j>= 0 && j < ip->lens[1] && 
        k>= 0 && k < ip->lens[2])   {
        int ijk = NKS*(i*ip->lens[1]*ip->lens[2] + j*ip->lens[2] + k);
        return poly(&ip->ks[ijk], 0, 1, 0, xp, yp, zp)/ip->dxs[1];
    }
    else    {
        return NAN;
    }
}

real diff3D_z(const Interpolation* ip, real x, real y, real z) {
    int i, j, k;
    real xp, yp, zp;

    xp_i(&xp, &i, x, ip->x0s[0], ip->dxs[0]);
    xp_i(&yp, &j, y, ip->x0s[1], ip->dxs[1]);
    xp_i(&zp, &k, z, ip->x0s[2], ip->dxs[2]);

    if(i>= 0 && i < ip->lens[0] && j>= 0 && j < ip->lens[1] && 
        k>= 0 && k < ip->lens[2])   {
        int ijk = NKS*(i*ip->lens[1]*ip->lens[2] + j*ip->lens[2] + k);
        return poly(&ip->ks[ijk], 0, 0, 1, xp, yp, zp)/ip->dxs[2];
    }
    else    {
        return NAN;
    }
}

real diff3D(const Interpolation* ip, int edx, int edy, int edz, 
              real x, real y, real z)    {
    int i, j, k;
    real xp, yp, zp;

    xp_i(&xp, &i, x, ip->x0s[0], ip->dxs[0]);
    xp_i(&yp, &j, y, ip->x0s[1], ip->dxs[1]);
    xp_i(&zp, &k, z, ip->x0s[2], ip->dxs[2]);

    if(i>= 0 && i < ip->lens[0] && j>= 0 && j < ip->lens[1] && 
        k>= 0 && k < ip->lens[2])   {
        int ijk = NKS*(i*ip->lens[1]*ip->lens[2] + j*ip->lens[2] + k);
        return poly(&ip->ks[ijk], edx, edy, edz, xp, yp, zp)/
            ipow(ip->dxs[0], edx)/ipow(ip->dxs[1], edy)/ipow(ip->dxs[2], edz);
    }
    else    {
        return NAN;
    }
}

void map3D(const Interpolation* ip, int len, const real xs[], 
           const real ys[], const real zs[], real out[])  {
    int i;
    for(i = 0; i < len; i++) {
        out[i]  = interpolate3D(ip, xs[i], ys[i], zs[i]);
    }
}

void map3D_x(const Interpolation* ip, int len, const real xs[], 
           const real ys[], const real zs[], real out[])  {
    int i;
    for(i = 0; i < len; i++) {
        out[i]  = diff3D_x(ip, xs[i], ys[i], zs[i]);
    }
}


void map3D_y(const Interpolation* ip, int len, const real xs[], 
           const real ys[], const real zs[], real out[])  {
    int i;
    for(i = 0; i < len; i++) {
        out[i]  = diff3D_y(ip, xs[i], ys[i], zs[i]);
    }
}

void map3D_z(const Interpolation* ip, int len, const real xs[], 
           const real ys[], const real zs[], real out[])  {
    int i;
    for(i = 0; i < len; i++) {
        out[i]  = diff3D_z(ip, xs[i], ys[i], zs[i]);
    }
}

void map3D_diff(const Interpolation *ip, int edx, int edy, int edz, int len,
                const real xs[], const real ys[], const real zs[], 
                real out[]) {
    int i;
    for(i = 0; i < len; i++) {
        out[i]  = diff3D(ip, edx, edy, edz, xs[i], ys[i], zs[i]);
    }
}

/* Private functions */

/* Create polynome.
 * Any modern compiler should unroll the loop. */
real poly(const real ks[NKS], int edx, int edy, int edz, 
                   real x, real y, real z)    {
    real r    = 0;
    for(int i = 0; i < 4; i++)  {
         real xx    = xi(x, edx, i);
         for(int j = 0; j < 4; j++)  {
            real xy    = xx*xi(y, edy, j);
            for(int k = 0; k < 4; k++)  {
                int ijk = i*16 + j*4 + k;
                r      += ks[ijk]*xy*xi(z, edz, k); 
            }
        }
    }
    return r;
}

void cellcoeffs(int ly, int lz, const real fs[], const real fs_x[],
                const real fs_y[], const real fs_z[], 
                const real fs_xy[], const real fs_xz[], 
                const real fs_yz[], const real fs_xyz[], real ks[NKS])   {

    int i000    = 0;
    int i001    = 1;
    int i010    = lz;
    int i011    = lz + 1;

    int i100    = ly*lz;
    int i101    = ly*lz + 1;
    int i110    = ly*lz + lz;
    int i111    = ly*lz + lz + 1;

    /* Maple generated code, do not touch */

    real t1 = 2 * fs_z[i000];
    real t2 = 3 * fs[i000];
    real t5 = 2 * fs[i000];
    real t8 = 2 * fs_yz[i000];
    real t9 = 3 * fs_y[i000];
    real t10 = 3 * fs_y[i001];
    real t12 = 2 * fs_y[i000];
    real t13 = 2 * fs_y[i001];
    real t17 = 3 * fs_z[i000];
    real t18 = 3 * fs_z[i010];
    real t20 = 6 * fs_y[i000];
    real t21 = 6 * fs_y[i001];
    real t22 = 3 * fs_y[i010];
    real t23 = 3 * fs_y[i011];
    real t24 = 4 * fs_yz[i000];
    real t25 = 2 * fs_yz[i001];
    real t26 = 2 * fs_yz[i010];
    real t27 = 6 * fs_z[i000];
    real t28 = 3 * fs_z[i001];
    real t29 = 6 * fs_z[i010];
    real t30 = 3 * fs_z[i011];
    real t31 = 9 * fs[i000];
    real t32 = 9 * fs[i001];
    real t33 = 9 * fs[i010];
    real t35 = t20 - t21 + t22 - t23 + t24 + t25 + t26 + fs_yz[i011] + t27 + t28 - t29 - t30 + t31 - t32 - t33 + 9 * fs[i011];
    real t36 = 4 * fs_y[i000];
    real t37 = 4 * fs_y[i001];
    real t38 = 6 * fs[i000];
    real t39 = 6 * fs[i001];
    real t40 = 6 * fs[i010];
    real t41 = 6 * fs[i011];
    real t42 = 2 * fs_y[i011];
    real t43 = 2 * fs_y[i010];
    real t44 = -t8 - t36 + t37 - t38 + t39 + t40 - t17 + t18 - t41 - fs_yz[i010] + t42 - t43 + t30 - t25 - t28 - fs_yz[i011];
    real t47 = 2 * fs_z[i010];
    real t49 = 4 * fs_z[i000];
    real t50 = 4 * fs_z[i010];
    real t51 = 2 * fs_z[i011];
    real t52 = 2 * fs_z[i001];
    real t53 = -t49 - t8 + t50 - t38 + t39 - t9 + t10 + t40 - t41 - t26 + t23 - t22 + t51 - fs_yz[i001] - t52 - fs_yz[i011];
    real t54 = 4 * fs[i000];
    real t55 = 4 * fs[i001];
    real t56 = 4 * fs[i010];
    real t58 = t54 - t55 + t12 + fs_yz[i000] - t13 - t56 + t1 - t47 + 4 * fs[i011] + fs_yz[i010] - t42 + t43 - t51 + fs_yz[i001] + t52 + fs_yz[i011];
    real t59 = 2 * fs_xz[i000];
    real t60 = 3 * fs_x[i000];
    real t61 = 3 * fs_x[i001];
    real t63 = 2 * fs_x[i000];
    real t64 = 2 * fs_x[i001];
    real t66 = 2 * fs_xyz[i000];
    real t67 = 3 * fs_xy[i000];
    real t68 = 3 * fs_xy[i001];
    real t70 = 2 * fs_xy[i000];
    real t71 = 2 * fs_xy[i001];
    real t73 = 3 * fs_x[i010];
    real t75 = 3 * fs_xz[i000];
    real t76 = 3 * fs_xz[i010];
    real t78 = 6 * fs_xz[i000];
    real t79 = 3 * fs_xz[i001];
    real t80 = 6 * fs_xz[i010];
    real t81 = 3 * fs_xz[i011];
    real t82 = 9 * fs_x[i000];
    real t83 = 9 * fs_x[i001];
    real t84 = 9 * fs_x[i010];
    real t85 = 9 * fs_x[i011];
    real t86 = 6 * fs_xy[i000];
    real t87 = 6 * fs_xy[i001];
    real t88 = 3 * fs_xy[i010];
    real t89 = 3 * fs_xy[i011];
    real t90 = 4 * fs_xyz[i000];
    real t91 = 2 * fs_xyz[i001];
    real t92 = 2 * fs_xyz[i010];
    real t93 = t78 + t79 - t80 - t81 + t82 - t83 - t84 + t85 + t86 - t87 + t88 - t89 + t90 + t91 + t92 + fs_xyz[i011];
    real t94 = 4 * fs_xy[i000];
    real t95 = 4 * fs_xy[i001];
    real t96 = 6 * fs_x[i000];
    real t97 = 6 * fs_x[i001];
    real t98 = 6 * fs_x[i010];
    real t99 = 6 * fs_x[i011];
    real t100 = 2 * fs_xy[i011];
    real t101 = 2 * fs_xy[i010];
    real t102 = -t66 - t94 + t95 - t96 + t97 + t98 - t75 + t76 - t99 - fs_xyz[i010] + t100 - t101 + t81 - t91 - t79 - fs_xyz[i011];
    real t103 = 2 * fs_x[i010];
    real t105 = 2 * fs_xz[i010];
    real t107 = 4 * fs_xz[i000];
    real t108 = 4 * fs_xz[i010];
    real t109 = 2 * fs_xz[i011];
    real t110 = 2 * fs_xz[i001];
    real t111 = -t107 - t66 + t108 - t96 + t97 - t67 + t68 + t98 - t99 - t92 + t89 - t88 + t109 - fs_xyz[i001] - t110 - fs_xyz[i011];
    real t112 = 4 * fs_x[i000];
    real t113 = 4 * fs_x[i001];
    real t114 = 4 * fs_x[i010];
    real t115 = 4 * fs_x[i011];
    real t116 = t112 - t113 + t70 + fs_xyz[i000] - t71 - t114 + t59 - t105 + t115 + fs_xyz[i010] - t100 + t101 - t109 + fs_xyz[i001] + t110 + fs_xyz[i011];
    real t119 = 3 * fs_z[i100];
    real t121 = 2 * fs_xz[i100];
    real t122 = 6 * fs_z[i100];
    real t123 = 3 * fs_z[i101];
    real t124 = 9 * fs[i100];
    real t126 = 3 * fs_x[i100];
    real t127 = 3 * fs_x[i101];
    real t128 = t107 + t110 + t121 + fs_xz[i101] + t27 + t28 - t122 - t123 + t31 - t32 - t124 + 9 * fs[i101] + t96 - t97 + t126 - t127;
    real t129 = 6 * fs[i100];
    real t130 = 6 * fs[i101];
    real t131 = 2 * fs_x[i101];
    real t132 = 2 * fs_x[i100];
    real t133 = -t59 - t112 + t113 - t38 + t39 + t129 - t17 + t119 - t130 - fs_xz[i100] + t131 - t132 + t123 - t110 - t28 - fs_xz[i101];
    real t134 = 3 * fs_y[i100];
    real t136 = 3 * fs_yz[i000];
    real t137 = 3 * fs_yz[i100];
    real t139 = 9 * fs_y[i000];
    real t140 = 9 * fs_y[i001];
    real t141 = 9 * fs_y[i100];
    real t142 = 9 * fs_y[i101];
    real t143 = 6 * fs_yz[i000];
    real t144 = 3 * fs_yz[i001];
    real t145 = 6 * fs_yz[i100];
    real t146 = 3 * fs_yz[i101];
    real t147 = 3 * fs_xy[i100];
    real t148 = 3 * fs_xy[i101];
    real t149 = 2 * fs_xyz[i100];
    real t150 = t139 - t140 - t141 + t142 + t143 + t144 - t145 - t146 + t86 - t87 + t147 - t148 + t90 + t91 + t149 + fs_xyz[i101];
    real t151 = 6 * fs_y[i100];
    real t152 = 6 * fs_y[i101];
    real t153 = 2 * fs_xy[i101];
    real t154 = 2 * fs_xy[i100];
    real t155 = -t66 - t94 + t95 - t20 + t21 + t151 - t136 + t137 - t152 - fs_xyz[i100] + t153 - t154 + t146 - t91 - t144 - fs_xyz[i101];
    real t156 = 3 * fs_y[i110];
    real t158 = 3 * fs_x[i110];
    real t159 = t20 + t22 - t151 - t156 + t31 - t33 - t124 + 9 * fs[i110] + t96 - t98 + t126 - t158 + t94 + t101 + t154 + fs_xy[i110];
    real t160 = 3 * fs_xz[i100];
    real t161 = 3 * fs_xz[i110];
    real t162 = 3 * fs_yz[i010];
    real t163 = 3 * fs_yz[i110];
    real t164 = 9 * fs_z[i000];
    real t165 = 9 * fs_z[i010];
    real t166 = 9 * fs_z[i100];
    real t167 = 9 * fs_z[i110];
    real t168 = t78 - t80 + t160 - t161 + t143 + t162 - t145 - t163 + t164 - t165 - t166 + t167 + t90 + t92 + t149 + fs_xyz[i110];
    real t170 = 9 * fs_x[i100];
    real t171 = 6 * fs_xy[i011];
    real t172 = 2 * fs_xyz[i011];
    real t174 = 4 * fs_xyz[i001];
    real t176 = 6 * fs_xy[i010];
    real t178 = 4 * fs_xyz[i010];
    real t179 = 6 * fs_xz[i011];
    real t184 = 6 * fs_xz[i001];
    real t185 = 27 * fs[i100] - t170 + t171 - t172 + 12 * fs_xy[i001] - t174 + 18 * fs_x[i010] - t176 + 12 * fs_xz[i010] - t178 + t179 - 18 * fs_x[i011] - 18 * fs_x[i000] - 12 * fs_xz[i000] + 18 * fs_x[i001] - t184;
    real t189 = 9 * fs_y[i010];
    real t191 = 6 * fs_yz[i010];
    real t192 = 9 * fs_y[i011];
    real t193 = 3 * fs_yz[i011];
    real t194 = 9 * fs_z[i011];
    real t196 = 9 * fs_x[i111];
    real t197 = 9 * fs_y[i111];
    real t198 = 9 * fs_z[i111];
    real t199 = 3 * fs_xy[i111];
    real t200 = 3 * fs_yz[i111];
    real t201 = 3 * fs_xz[i111];
    real t202 = -12 * fs_xy[i000] - 8 * fs_xyz[i000] + 27 * fs[i010] - t189 + 18 * fs_z[i010] - t191 + t192 - t193 + t194 - 27 * fs[i011] - t196 - t197 - t198 + t199 + t200 + t201;
    real t204 = 6 * fs_yz[i110];
    real t206 = 2 * fs_xyz[i110];
    real t208 = 9 * fs_y[i110];
    real t210 = 9 * fs_x[i110];
    real t211 = 3 * fs_xy[i110];
    real t212 = 6 * fs_xz[i110];
    real t214 = 6 * fs_xy[i100];
    real t216 = 4 * fs_xyz[i100];
    real t218 = 6 * fs_yz[i101];
    real t219 = -fs_xyz[i111] + t204 - 18 * fs_z[i110] - t206 + 27 * fs[i111] + t208 - 27 * fs[i110] + t210 - t211 + t212 + 18 * fs_y[i100] - t214 + 12 * fs_yz[i100] - t216 - 18 * fs_y[i101] + t218;
    real t220 = 6 * fs_xy[i101];
    real t221 = 2 * fs_xyz[i101];
    real t223 = 6 * fs_xz[i100];
    real t224 = 3 * fs_xz[i101];
    real t225 = 9 * fs_z[i101];
    real t227 = 9 * fs_x[i101];
    real t231 = 9 * fs_z[i001];
    real t235 = 6 * fs_yz[i001];
    real t236 = t220 - t221 + 18 * fs_z[i100] - t223 - t224 + t225 - 27 * fs[i101] + t227 - 27 * fs[i000] - 18 * fs_z[i000] + 27 * fs[i001] - t231 - 18 * fs_y[i000] - 12 * fs_yz[i000] + 18 * fs_y[i001] - t235;
    real t239 = 18 * fs[i100];
    real t240 = 6 * fs_x[i100];
    real t241 = 4 * fs_xy[i011];
    real t243 = 12 * fs_x[i010];
    real t244 = 4 * fs_xy[i010];
    real t245 = 12 * fs_x[i011];
    real t246 = 12 * fs_x[i000];
    real t247 = 12 * fs_x[i001];
    real t248 = -t239 + t240 - t241 + t172 - 8 * fs_xy[i001] + t174 - t243 + t244 - t80 + t92 - t179 + t245 + t246 + t78 - t247 + t184;
    real t250 = 18 * fs[i010];
    real t251 = 6 * fs_y[i010];
    real t252 = 6 * fs_y[i011];
    real t253 = 18 * fs[i011];
    real t254 = 6 * fs_x[i111];
    real t255 = 6 * fs_y[i111];
    real t256 = 2 * fs_xy[i111];
    real t257 = 8 * fs_xy[i000] + t90 - t250 + t251 - t165 + t162 - t252 + t193 - t194 + t253 + t254 + t255 + t198 - t256 - t200 - t201;
    real t259 = 18 * fs[i111];
    real t260 = 6 * fs_y[i110];
    real t261 = 18 * fs[i110];
    real t262 = 6 * fs_x[i110];
    real t263 = 2 * fs_xy[i110];
    real t264 = 12 * fs_y[i100];
    real t265 = 4 * fs_xy[i100];
    real t266 = 12 * fs_y[i101];
    real t267 = fs_xyz[i111] - t163 + t167 + fs_xyz[i110] - t259 - t260 + t261 - t262 + t263 - t161 - t264 + t265 - t145 + t149 + t266 - t218;
    real t268 = 4 * fs_xy[i101];
    real t269 = 18 * fs[i101];
    real t270 = 6 * fs_x[i101];
    real t271 = 18 * fs[i000];
    real t272 = 18 * fs[i001];
    real t273 = 12 * fs_y[i000];
    real t274 = 12 * fs_y[i001];
    real t275 = -t268 + t221 - t166 + t160 + t224 - t225 + t269 - t270 + t271 + t164 - t272 + t231 + t273 + t143 - t274 + t235;
    real t278 = 6 * fs[i110];
    real t279 = 2 * fs_x[i110];
    real t280 = t129 - t9 - t70 + t134 - t38 + t40 - t112 + t114 - t278 - t101 + t156 - t22 + t279 - t132 - fs_xy[i100] - fs_xy[i110];
    real t281 = 6 * fs_z[i110];
    real t282 = 2 * fs_xz[i110];
    real t283 = -t107 + t122 - t136 - t66 + t137 - t27 + t29 + t108 - t281 - t92 + t163 - t162 + t282 - fs_xyz[i100] - t121 - fs_xyz[i110];
    real t285 = 4 * fs_xz[i011];
    real t287 = 4 * fs_xz[i001];
    real t288 = -t239 + t240 - t171 + t172 - t87 + t91 - t243 + t176 - 8 * fs_xz[i010] + t178 - t285 + t245 + t246 + 8 * fs_xz[i000] - t247 + t287;
    real t289 = 12 * fs_z[i010];
    real t290 = 6 * fs_z[i011];
    real t291 = 6 * fs_z[i111];
    real t292 = 2 * fs_xz[i111];
    real t293 = t86 + t90 - t250 + t189 - t289 + t191 - t192 + t193 - t290 + t253 + t254 + t197 + t291 - t199 - t200 - t292;
    real t295 = 12 * fs_z[i110];
    real t296 = 4 * fs_xz[i110];
    real t297 = fs_xyz[i111] - t204 + t295 + t206 - t259 - t208 + t261 - t262 + t211 - t296 - t141 + t147 - t145 + t149 + t142 - t146;
    real t298 = 12 * fs_z[i100];
    real t299 = 4 * fs_xz[i100];
    real t300 = 2 * fs_xz[i101];
    real t301 = 6 * fs_z[i101];
    real t302 = 12 * fs_z[i000];
    real t303 = 6 * fs_z[i001];
    real t304 = -t148 + fs_xyz[i101] - t298 + t299 + t300 - t301 + t269 - t270 + t271 + t302 - t272 + t303 + t139 + t143 - t140 + t144;
    real t307 = 12 * fs[i100];
    real t308 = 4 * fs_x[i100];
    real t313 = t307 - t308 + t241 - t172 + t95 - t91 + 8 * fs_x[i010] - t244 + t108 - t92 + t285 - 8 * fs_x[i011] - 8 * fs_x[i000] - t107 + 8 * fs_x[i001] - t287;
    real t314 = 12 * fs[i010];
    real t315 = 12 * fs[i011];
    real t316 = 4 * fs_x[i111];
    real t317 = -t94 - t66 + t314 - t251 + t29 - t162 + t252 - t193 + t290 - t315 - t316 - t255 - t291 + t256 + t200 + t292;
    real t319 = 12 * fs[i111];
    real t320 = 12 * fs[i110];
    real t321 = 4 * fs_x[i110];
    real t322 = -fs_xyz[i111] + t163 - t281 - fs_xyz[i110] + t319 + t260 - t320 + t321 - t263 + t282 + t151 - t154 + t137 - fs_xyz[i100] - t152 + t146;
    real t323 = 12 * fs[i101];
    real t324 = 4 * fs_x[i101];
    real t325 = 12 * fs[i000];
    real t326 = 12 * fs[i001];
    real t327 = t153 - fs_xyz[i101] + t122 - t121 - t300 + t301 - t323 + t324 - t325 - t27 + t326 - t303 - t20 - t136 + t21 - t144;
    real t332 = 2 * fs_z[i100];
    real t334 = 4 * fs_z[i100];
    real t335 = 2 * fs_z[i101];
    real t336 = -t49 - t59 + t334 - t38 + t39 - t60 + t61 + t129 - t130 - t121 + t127 - t126 + t335 - fs_xz[i001] - t52 - fs_xz[i101];
    real t337 = 4 * fs[i100];
    real t339 = t54 - t55 + t63 + fs_xz[i000] - t64 - t337 + t1 - t332 + 4 * fs[i101] + fs_xz[i100] - t131 + t132 - t335 + fs_xz[i001] + t52 + fs_xz[i101];
    real t340 = 2 * fs_y[i100];
    real t342 = 2 * fs_yz[i100];
    real t344 = 4 * fs_yz[i100];
    real t345 = 2 * fs_yz[i101];
    real t346 = -t24 - t66 + t344 - t20 + t21 - t67 + t68 + t151 - t152 - t149 + t148 - t147 + t345 - fs_xyz[i001] - t25 - fs_xyz[i101];
    real t347 = 4 * fs_y[i100];
    real t348 = 4 * fs_y[i101];
    real t349 = t36 - t37 + t70 + fs_xyz[i000] - t71 - t347 + t8 - t342 + t348 + fs_xyz[i100] - t153 + t154 - t345 + fs_xyz[i001] + t25 + fs_xyz[i101];
    real t350 = 2 * fs_y[i110];
    real t351 = -t36 - t70 + t347 + t129 - t38 + t40 - t60 + t73 - t278 - fs_xy[i010] + t350 - t43 + t158 - t126 - t154 - fs_xy[i110];
    real t352 = 2 * fs_yz[i110];
    real t353 = -t24 - t66 + t344 - t75 + t122 - t27 + t29 + t76 - t281 - fs_xyz[i010] + t352 - t26 + t161 - t149 - t160 - fs_xyz[i110];
    real t354 = -t239 + t170 - t89 + fs_xyz[i011] - t87 + t91 - t84 + t88 - t80 + t92 - t81 + t85 + t82 + t78 - t83 + t79;
    real t355 = 4 * fs_yz[i010];
    real t356 = 2 * fs_yz[i011];
    real t357 = 2 * fs_yz[i111];
    real t358 = t86 + t90 - t250 + t251 - t289 + t355 - t252 + t356 - t290 + t253 + t196 + t255 + t291 - t199 - t357 - t201;
    real t360 = 4 * fs_yz[i110];
    real t362 = 4 * fs_yz[i101];
    real t363 = fs_xyz[i111] - t360 + t295 + t206 - t259 - t260 + t261 - t210 + t211 - t212 - t264 + t214 - 8 * fs_yz[i100] + t216 + t266 - t362;
    real t365 = 4 * fs_yz[i001];
    real t366 = -t220 + t221 - t298 + t223 + t224 - t301 + t269 - t227 + t271 + t302 - t272 + t303 + t273 + 8 * fs_yz[i000] - t274 + t365;
    real t369 = t307 - t240 + t100 - fs_xyz[i011] + t95 - t91 + t98 - t101 + t76 - fs_xyz[i010] + t81 - t99 - t96 - t75 + t97 - t79;
    real t370 = 4 * fs_y[i010];
    real t371 = 4 * fs_y[i011];
    real t372 = 4 * fs_y[i111];
    real t373 = -t94 - t66 + t314 - t370 + t29 - t26 + t371 - t356 + t290 - t315 - t254 - t372 - t291 + t256 + t357 + t201;
    real t375 = 4 * fs_y[i110];
    real t378 = -fs_xyz[i111] + t352 - t281 - fs_xyz[i110] + t319 + t375 - t320 + t262 - t263 + t161 + 8 * fs_y[i100] - t265 + t344 - t149 - 8 * fs_y[i101] + t362;
    real t381 = t268 - t221 + t122 - t160 - t224 + t301 - t323 + t270 - t325 - t27 + t326 - t303 - 8 * fs_y[i000] - t24 + 8 * fs_y[i001] - t365;
    real t385 = -t337 + t12 + fs_xy[i000] - t340 + t54 - t56 + t63 - t103 + 4 * fs[i110] + fs_xy[i010] - t350 + t43 - t279 + t132 + fs_xy[i100] + fs_xy[i110];
    real t386 = 4 * fs_z[i110];
    real t387 = t59 - t334 + t8 + fs_xyz[i000] - t342 + t49 - t50 - t105 + t386 + fs_xyz[i010] - t352 + t26 - t282 + fs_xyz[i100] + t121 + fs_xyz[i110];
    real t388 = t307 - t240 + t89 - fs_xyz[i011] + t68 - fs_xyz[i001] + t98 - t88 + t108 - t92 + t109 - t99 - t96 - t107 + t97 - t110;
    real t390 = 4 * fs_z[i011];
    real t391 = 4 * fs_z[i111];
    real t392 = -t67 - t66 + t314 - t251 + 8 * fs_z[i010] - t355 + t252 - t356 + t390 - t315 - t254 - t255 - t391 + t199 + t357 + t292;
    real t395 = -fs_xyz[i111] + t360 - 8 * fs_z[i110] - t206 + t319 + t260 - t320 + t262 - t211 + t296 + t151 - t147 + t344 - t149 - t152 + t345;
    real t397 = 4 * fs_z[i101];
    real t399 = 4 * fs_z[i001];
    real t400 = t148 - fs_xyz[i101] + 8 * fs_z[i100] - t299 - t300 + t397 - t323 + t270 - t325 - 8 * fs_z[i000] + t326 - t399 - t20 - t24 + t21 - t25;
    real t404 = -8 * fs[i100] + t308 - t100 + fs_xyz[i011] - t71 + fs_xyz[i001] - t114 + t101 - t105 + fs_xyz[i010] - t109 + t115 + t112 + t59 - t113 + t110;
    real t407 = t70 + fs_xyz[i000] - 8 * fs[i010] + t370 - t50 + t26 - t371 + t356 - t390 + 8 * fs[i011] + t316 + t372 + t391 - t256 - t357 - t292;
    real t411 = fs_xyz[i111] - t352 + t386 + fs_xyz[i110] - 8 * fs[i111] - t375 + 8 * fs[i110] - t321 + t263 - t282 - t347 + t154 - t342 + fs_xyz[i100] + t348 - t345;
    real t415 = -t153 + fs_xyz[i101] - t334 + t121 + t300 - t397 + 8 * fs[i101] - t324 + 8 * fs[i000] + t49 - 8 * fs[i001] + t399 + t36 + t8 - t37 + t25;

    ks[0] = fs[i000];
    ks[1] = fs_z[i000];
    ks[2] = -t1 - t2 + 3 * fs[i001] - fs_z[i001];
    ks[3] = t5 + fs_z[i000] - 2 * fs[i001] + fs_z[i001];
    ks[4] = fs_y[i000];
    ks[5] = fs_yz[i000];
    ks[6] = -t8 - t9 + t10 - fs_yz[i001];
    ks[7] = t12 + fs_yz[i000] - t13 + fs_yz[i001];
    ks[8] = -t12 - t2 + 3 * fs[i010] - fs_y[i010];
    ks[9] = -t8 - t17 + t18 - fs_yz[i010];
    ks[10] = t35;
    ks[11] = t44;
    ks[12] = t5 + fs_y[i000] - 2 * fs[i010] + fs_y[i010];
    ks[13] = t1 + fs_yz[i000] - t47 + fs_yz[i010];
    ks[14] = t53;
    ks[15] = t58;
    ks[16] = fs_x[i000];
    ks[17] = fs_xz[i000];
    ks[18] = -t59 - t60 + t61 - fs_xz[i001];
    ks[19] = t63 + fs_xz[i000] - t64 + fs_xz[i001];
    ks[20] = fs_xy[i000];
    ks[21] = fs_xyz[i000];
    ks[22] = -t66 - t67 + t68 - fs_xyz[i001];
    ks[23] = t70 + fs_xyz[i000] - t71 + fs_xyz[i001];
    ks[24] = -t70 - t60 + t73 - fs_xy[i010];
    ks[25] = -t66 - t75 + t76 - fs_xyz[i010];
    ks[26] = t93;
    ks[27] = t102;
    ks[28] = t63 + fs_xy[i000] - t103 + fs_xy[i010];
    ks[29] = t59 + fs_xyz[i000] - t105 + fs_xyz[i010];
    ks[30] = t111;
    ks[31] = t116;
    ks[32] = -t63 - t2 + 3 * fs[i100] - fs_x[i100];
    ks[33] = -t59 - t17 + t119 - fs_xz[i100];
    ks[34] = t128;
    ks[35] = t133;
    ks[36] = -t70 - t9 + t134 - fs_xy[i100];
    ks[37] = -t66 - t136 + t137 - fs_xyz[i100];
    ks[38] = t150;
    ks[39] = t155;
    ks[40] = t159;
    ks[41] = t168;
    ks[42] = t185 + t202 + t219 + t236;
    ks[43] = t248 + t257 + t267 + t275;
    ks[44] = t280;
    ks[45] = t283;
    ks[46] = t288 + t293 + t297 + t304;
    ks[47] = t313 + t317 + t322 + t327;
    ks[48] = t5 + fs_x[i000] - 2 * fs[i100] + fs_x[i100];
    ks[49] = t1 + fs_xz[i000] - t332 + fs_xz[i100];
    ks[50] = t336;
    ks[51] = t339;
    ks[52] = t12 + fs_xy[i000] - t340 + fs_xy[i100];
    ks[53] = t8 + fs_xyz[i000] - t342 + fs_xyz[i100];
    ks[54] = t346;
    ks[55] = t349;
    ks[56] = t351;
    ks[57] = t353;
    ks[58] = t354 + t358 + t363 + t366;
    ks[59] = t369 + t373 + t378 + t381;
    ks[60] = t385;
    ks[61] = t387;
    ks[62] = t388 + t392 + t395 + t400;
    ks[63] = t404 + t407 + t411 + t415;
}
