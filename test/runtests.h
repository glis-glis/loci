/* Copyright (C) 2017 Andreas Füglistaler <andreas.fueglistaler@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef FILE_RUNTESTS_H
#define FILE_RUNTESTS_H

static const double tols[5]   = {1e-10, 1e-3, 1e-2, 1e-1, 1};

void basic1D(void);
void exp1D(void);
void boundary1D(void);

void basic2D(void);
void exp2D(void);
void boundary2D(void);

void basic3D(void);
void exp3D(void);
void boundary3D(void);

void basic4D(void);
void exp4D(void);
void boundary4D(void);

#endif
