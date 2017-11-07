/* Copyright (C) 2017 Andreas Füglistaler <andreas.fueglistaler@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef FILE_TESTHELPERS_H
#define FILE_TESTHELPERS_H

#include <math.h>
#include <stdio.h>

/* Output macro, print _L tabs, followed by printf arguments. */
#define OUT(_L, ...) {                   \
        int _i;                         \
        for(_i = 0; _i < _L; _i++)  {   \
            printf("\t");               \
        }                               \
        printf(__VA_ARGS__);            \
        printf("\n");                   \
    } 

/* maximum */
#define MAX(X, Y) ((X) > (Y) ? (X) : (Y))
/* square */
#define SQ(X) ((X)*(X))
/* cube */
#define CU(X) ((X)*SQ(X))
/* quatruple */
#define QU(X) (SQ(X)*SQ(X))

/* calculate the norm of `xs` - `ys` with array length `len`. */
static double errnorm(int len, const double xs[], const double ys[]) {
    int i;
    double err  = 0;
    for(i = 0; i < len; i++) {
        double e    = xs[i] - ys[i];
        err += e*e;
    }
    return sqrt(err);
}

/* calculate the norm of  the array `xs` with length `len`. */
static double norm(int len, const double xs[])  {
    int i;
    double n    = 0;
    for(i = 0; i < len; i++) {
        n   += xs[i]*xs[i];
    }
    return sqrt(n);
}

#endif
