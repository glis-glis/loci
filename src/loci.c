/* Copyright (C) 2017 Andreas Füglistaler <andreas.fueglistaler@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <stdlib.h>

#include "../include/loci.h"

/* "Destructor" */
void loci_free(Interpolation *ip)    {
    free(ip->x0s);
    free(ip->dxs);
    free(ip->lens);
    free(ip->ks);
}

/* Return sizeof(real) */
int loci_size(void) {
    return sizeof(real);
}
