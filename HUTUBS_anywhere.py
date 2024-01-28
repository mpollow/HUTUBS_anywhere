#!/usr/bin/env python

import numpy as np
import math
import h5py
from scipy.io import loadmat
from scipy.special import lpmv


def acn2nm(acn):
    """
    Compute (n, m) index from Ambisonics channel numbering (ACN)
    >>> acn2nm(5)
    (2, -1)
    """
    t = type(acn)
    acn = np.array(acn)
    n = (acn**0.5).astype(int)
    m = (np.round(acn - n**2 - n)).astype(int)
    if t is list:
        n, m = list(n), list(m)
    return n, m

def pos2tp(pos):
    # converts position (az, el, r) to theta and phi
    az, el, _ = pos.T
    theta = (90 - el) / 180 * np.pi
    phi = az / 180 * np.pi
    return theta, phi

def SH_complex(nmax, theta, phi):
    nSH = (nmax + 1)**2
    SH = np.zeros((len(theta), nSH), dtype=np.complex128)
    for acn in np.arange(nSH):
        n, m = acn2nm(acn)
        fac0 = math.factorial(n - m)
        fac1 = math.factorial(n + m)
        N = ( (2*n + 1) / (4 * np.pi) * fac0 / fac1 )**0.5
        P = lpmv(m, n, np.cos(theta))
        X = np.exp(1j * m * phi)
        SH[:, acn] = N * P * X
    return SH

def SH_real(nmax, theta, phi):
    nSH = (nmax + 1)**2
    SH = np.zeros((len(theta), nSH), dtype=np.float64)
    for acn in np.arange(nSH):
        n, m = acn2nm(acn)
        fac0 = math.factorial(n - m)
        fac1 = math.factorial(n + m)
        N = ( (2*n + 1) / (4 * np.pi) * fac0 / fac1 )**0.5
        P = lpmv(m, n, np.cos(theta)) * (-1.0)**m
        if m < 0:
            X = np.sin(abs(m) * phi) * (-1.0)**(m-1)
        else:
            X = np.cos(abs(m) * phi)
        if m != 0:
            X *= 2.0**0.5
        SH[:, acn] = N * P * X
    return SH

def get_HRIRs(SOFA_filename):
    h5 = h5py.File(SOFA_filename)
    hrir = h5['Data.IR'][:]
    sr = h5['Data.SamplingRate'][0]
    pos = h5['SourcePosition'][:]
    return hrir, sr, pos


def get_HRTFs_fromSH(SH_filename, theta, phi):
    a = loadmat(SH_filename)
    
    hrtfSH0, hrtfSH1, _nmax, n, m, _freqs, _sr, unit, comment = a['HRIR'][0][0]

    try:
        # TOA compensation only used for HRTF measurements, not for simulations
        toa0sh = a['TOA'][0][0][0]
        toa1sh = a['TOA'][0][0][1]
        toa_compensation = True
    except:
        toa_compensation = False
    
    # get rid of the extra dimensions
    nmax = _nmax[0][0]
    sr = _sr[0][0]
    freqs = _freqs[0]
   
    # verify number of SH coefficients and frequency sampling
    assert (nmax+1)**2 == hrtfSH0.shape[0]
    assert np.all(np.linspace(0, sr/2, len(freqs)) == freqs)

    sh = SH_complex(nmax, theta, phi)
    
    hrtf0 = sh @ hrtfSH0
    hrtf1 = sh @ hrtfSH1

    if toa_compensation:
        sh_real = SH_real(nmax, theta, phi)
        t0 = (sh_real @ toa0sh) / sr
        t1 = (sh_real @ toa1sh) / sr

        omega = 2 * np.pi * freqs
        hrtf0 = hrtf0 * np.exp(-1j * omega * t0)
        hrtf1 = hrtf1 * np.exp(-1j * omega * t1)

    hrtf = np.stack((hrtf0, hrtf1), axis=1)
    #hrir = np.fft.irfft(hrtf)

    return hrtf, sr

