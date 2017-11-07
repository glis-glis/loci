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
import numpy as np
from . import CInterpolation

class Interpolation1D:
    """Interpolation in 1 dimension. """ 

    def __init__(self, rx, f, f_x):
        """Create interpolation in range `rx` for values `fs`= f(`rx`) and 
        derivative values `f_x` = df/dx(`rx`). 
         
        Given a range r with x0, dx and len
        g(r) = {g(x0), g(x0 + dx), ..., g(x0+(len-1)*dx). 
        """
        self._lib                          = ct.CDLL("libloci.so")

        lib                          = self._lib
        lib.interpolation1D.restype  = CInterpolation
        lib.interpolate1D.restype    = ct.c_double
        lib.diff1D_x.restype         = ct.c_double
        lib.diff1D.restype           = ct.c_double
        lib.map1D.restype            = None
        lib.map1D_x.restype          = None
        lib.map1D_diff.restype       = None
        
        if callable(f): 
            #expects f AND f_x to be callable
            xs      = np.r_[rx.x0:rx.x0 + rx.len*rx.dx:rx.dx]
            fs      = np.array([f(xi) for xi in xs])
            fs_x    = np.array([f_x(xi) for xi in xs])
        else:
            fs      = f
            fs_x    = f_x

        self._ip = lib.interpolation1D(ct.pointer(rx),
                     fs.ctypes.data_as(ct.POINTER(ct.c_double)),
                     fs_x.ctypes.data_as(ct.POINTER(ct.c_double)))

    def interpolate(self, x):
        """Interpolate at point `x`. """
        return self._lib.interpolate1D(ct.pointer(self._ip), ct.c_double(x))

    def diff_x(self, x):
        """Interpolate first order derivative in x at point `x`. """
        return self._lib.diff1D_x(ct.pointer(self._ip), ct.c_double(x))

    def diff(self, edx, x):
        """Interpolate `edx`-order derivative in x at point `x`. """
        return self._lib.diff1D(ct.pointer(self._ip), ct.c_int(edx), 
                                ct.c_double(x))

    def map(self, xs):
        """Map interpolation onto input-array `xs` and return result. """
        l   = len(xs)
        out = np.empty(l)
        self._lib.map1D(ct.pointer(self._ip), ct.c_int(l), 
            xs.ctypes.data_as(ct.POINTER(ct.c_double)),
            out.ctypes.data_as(ct.POINTER(ct.c_double)))
        return out

    def map_x(self, xs):
        """Map first order derivative in x onto input-array `xs` 
        and return result.
        """
        l   = len(xs)
        out = np.empty(l)
        self._lib.map1D_x(ct.pointer(self._ip), ct.c_int(l), 
            xs.ctypes.data_as(ct.POINTER(ct.c_double)),
            out.ctypes.data_as(ct.POINTER(ct.c_double)))
        return out

    def map_diff(self, edx, xs):
        """Map `edx`-order derivative in x onto input-array `xs` 
        and return result.
        """
        l   = len(xs)
        out = np.empty(l)
        self._lib.map1D_diff(ct.pointer(self._ip), ct.c_int(edx), ct.c_int(l), 
            xs.ctypes.data_as(ct.POINTER(ct.c_double)),
            out.ctypes.data_as(ct.POINTER(ct.c_double)))
        return out

    def __del__(self):
        self._lib.loci_free(ct.pointer(self._ip))
