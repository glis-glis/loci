/* Copyright (C) 2017 Andreas Füglistaler <andreas.fueglistaler@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <stdio.h>

#include "runtests.h"
#include "../include/loci.h"

int main(void)  {

    puts("loci: Local Interpolations!\n");

    basic1D();
    exp1D();
    boundary1D();

    basic2D();
    exp2D();
    boundary2D();

    basic3D();
    exp3D();
    boundary3D();

    basic4D();
    exp4D();
    boundary4D();
    return 0;
}
