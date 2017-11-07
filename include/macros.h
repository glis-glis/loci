/* Copyright (C) 2017 Andreas Füglistaler <andreas.fueglistaler@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef FILE_MACROS_H
#define FILE_MACROS_H

#define _NARGS_HELPER( _1,  _2,  _3,  _4,  _5,  _6,  _7,  _8,  _9, _10, \
                      _11, _12, _13, _14, _15, _16, _17, _18, _19, _20, \
                      N, ...) N
#define _NARGS(...)                                                      \
    _NARGS_HELPER(__VA_ARGS__, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11,  \
                  10, 9, 8, 7, 6, 5, 4, 3, 2, 1)

#define _CONC(X, Y) _DCONC(X, Y)
#define _DCONC(X, Y) X ## Y

#define _INTERPOLATION3 interpolation1D
#define _INTERPOLATION6 interpolation2D
#define _INTERPOLATION11 interpolation3D
#define _INTERPOLATION20 interpolation4D

#define _INTERPOLATE2 interpolate1D
#define _INTERPOLATE3 interpolate2D
#define _INTERPOLATE4 interpolate3D
#define _INTERPOLATE5 interpolate4D

#define _DIFF_X2 diff1D_x
#define _DIFF_X3 diff2D_x
#define _DIFF_X4 diff3D_x
#define _DIFF_X5 diff4D_x

#define _DIFF_Y3 diff2D_y
#define _DIFF_Y4 diff3D_y
#define _DIFF_Y5 diff4D_y

#define _DIFF_Z4 diff3D_z
#define _DIFF_Z5 diff4D_z

#define _DIFF3 diff1D
#define _DIFF5 diff2D
#define _DIFF7 diff3D
#define _DIFF9 diff4D

#define _MAP4 map1D
#define _MAP5 map2D
#define _MAP6 map3D
#define _MAP7 map4D

#define _MAP_X4 map1D_x
#define _MAP_X5 map2D_x
#define _MAP_X6 map3D_x
#define _MAP_X7 map4D_x

#define _MAP_Y5 map2D_y
#define _MAP_Y6 map3D_y
#define _MAP_Y7 map4D_y

#define _MAP_Z6 map3D_z
#define _MAP_Z7 map4D_z

#define _MAP_DIFF5 map1D_diff
#define _MAP_DIFF7 map2D_diff
#define _MAP_DIFF9 map3D_diff
#define _MAP_DIFF11 map4D_diff

#endif
