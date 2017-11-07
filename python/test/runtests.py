# -*- coding: utf-8 -*-
# loci
# Local cubic interpolations in up to 4 dimensions
#
# Copyright (C) 2017 Andreas Füglistaler <andreas.fueglistaler@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from loci import Interpolation
from loci import Range
from numpy import *
from numpy.linalg import norm
from scipy import *

from test1D import basic1D, exp1D, boundary1D
from test2D import basic2D, exp2D, boundary2D
from test3D import basic3D, exp3D, boundary3D
from test4D import basic4D, exp4D, boundary4D

basic1D()
exp1D()
boundary1D()

basic2D()
exp2D()
boundary2D()

basic3D()
exp3D()
boundary3D()

basic4D()
exp4D()
boundary4D()
