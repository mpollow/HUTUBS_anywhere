#!/usr/bin/env python

import numpy as np

import coordinates as co
import sphharm as sh

NMAX = 5
cart = co.gaussian(NMAX)
np.savetxt(f'gaussian{NMAX}.txt', cart)

