# -*- coding: utf-8 -*-
# loci
# Local cubic interpolations in up to 4 dimensions
#
# Copyright (C) 2017 Andreas Füglistaler <andreas.fueglistaler@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import ctypes as ct


class Range(ct.Structure):
    """Range from `x0` to `x0`+(`len`-1)*`dx`, with step `dx`and length `len`.
    """
    _fields_    = [("x0", ct.c_double),
                   ("dx", ct.c_double),
                   ("len", ct.c_int)]

class CInterpolation(ct.Structure):
    """Interpolation in up to 4 dimensions dim on a grid with interpolation 
    coeffisients ks. The grid parameters `x0s`[dim], `dxs`[dim] and `lens`[dim]
    have the same definition as in the class `Range`.
    """
    _fields_    = [("x0s", ct.POINTER(ct.c_double)),
                   ("dxs", ct.POINTER(ct.c_double)),
                   ("lens", ct.POINTER(ct.c_int)),
                   ("ks", ct.POINTER(ct.c_double))]

from Interpolation1D import Interpolation1D
from Interpolation2D import Interpolation2D
from Interpolation3D import Interpolation3D
from Interpolation4D import Interpolation4D

def Interpolation(*args):
    """interpolation(...) with NARGS = len(args)
        NARGS= 3: interpolation1D(args)
        NARGS= 6: interpolation2D(args)
        NARGS=11: interpolation3D(args)
        NARGS=20: interpolation4D(args). 
    """
    l   = len(args)
    if l == 3:
        return Interpolation1D(*args)
    elif l == 6:
        return Interpolation2D(*args)
    elif l == 11:
        return Interpolation3D(*args)
    elif l == 20:
        return Interpolation4D(*args)
    else:
        raise RuntimeError("Wrong number of arguments")
