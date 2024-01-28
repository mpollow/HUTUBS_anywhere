#!/usr/bin/env python

import pytest

import numpy as np
import os
import glob

import sys
sys.path.append('..')

from HUTUBS_anywhere import *


PATH = '../HUTUBS_HRIRs'

def check_HRTF_match(SH_filename):
    SOFA_filename = SH_filename.replace('SHcoefficients', 'HRIRs').replace('.mat', '.sofa')
    hrir, sr_sofa, pos = get_HRIRs(SOFA_filename)
    hrtf = np.fft.rfft(hrir)

    theta, phi = pos2tp(pos)
    hrtf_SH, sr_mat = get_HRTFs_fromSH(SH_filename, theta, phi)

    assert sr_sofa == sr_mat
    assert np.allclose(hrtf, hrtf_SH, rtol=np.inf, atol=1e-8)
    # test only for absolute tolerance

def test_all_HRIRs(path=PATH):
    for filepath in glob.glob(os.path.join(PATH, '*.mat')):
        check_HRTF_match(filepath)

def test_SH_basis():
    nmax = 5
    cart = np.loadtxt(f'./coordinates_gaussian{nmax}.txt')
    a = loadmat(f'./SH_AKtools_gaussian{nmax}.mat')
    SHc_AK = a['Ycomplex']
    SHr_AK = a['Yreal']
    
    x, y, z = cart.T
    theta = np.arctan2((x**2 + y**2)**0.5, z)
    phi = np.arctan2(y, x)

    SHc = SH_complex(nmax, theta, phi)
    SHr = SH_real(nmax, theta, phi)

    assert np.allclose(SHc, SHc_AK)
    assert np.allclose(SHr, SHr_AK)


if __name__=="__main__":
    retcode = pytest.main()

