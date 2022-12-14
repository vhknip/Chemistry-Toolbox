# -*- coding: utf-8 -*-
"""
Spyder Editor

@author: vincent.nip
"""
import numpy as np
def vector_norm(vector: list) -> list:
    """
    

    Parameters
    ----------
    vector : list
        This function reads in a vector as a list.

    Returns
    -------
    list
        This function will return a normalized vector as a list.
        If a 0-vector is passed in, the original vector is returned.

    """
    magn = np.linalg.norm(vector)
    if magn != 0:
        normalized = [x / magn for x in vector]
        return normalized
    # deals with 0-vector passed in
    else:
        return vector

def vector_map(vector: list, basis: set, target_basis: list) -> list:
    """
    

    Parameters
    ----------
    vector : list
        Vector as a list that user wants to map to a different basis
    basis : list
        Basis as a list that the original vector space is spanned by
    target_basis : list
        Target basis as a list that user wants to map vector to

    Returns
    -------
    mapped_vector: list
        Original vector exclusively mapped on to new specified basis

    """
    
    mapped_vector = [0] * len(target_basis)
    for x in range(0,len(target_basis)):
        for y in range(0,len(basis)):
            if basis[y] == target_basis[x]:
                mapped_vector[x] = vector[y]
    return mapped_vector
                
    
    
def match(sample_sig: list,
          sample_mz: list,
          ref_sig: list,
          ref_mz: list,
          mode: str = 'direct',
          mult: float = 100) -> float:
    """
    

    Parameters
    ----------
    sample_sig : list
        Signal values of each m/z value detected in sample compound as list
    sample_mz : list
        m/z values corresponding to each signal in sample compound as list
    ref_sig : list
        Signal values of each m/z value detected in reference compound as list
    ref_mz : list
        m/z values corresponding to each signal in reference compound as list
    mode : {'direct', 'forward', 'reverse'}
        modes for matching that include 'direct','forward', and 'reverse'
        direct: normal match factor considering all values
        forward: match factor considering only non-zero values found in sample
            used for low concentrations / low ion abundance samples
        reverse: match factor considering only non-zero values found in reference
            used for high noise spectra
    mult: float
        Multiplier used for match factor score.
        Agilent (Unknowns Analysis and OpenLAB CDS) use 100 (default)
        Thermo Fisher (Chromeleon, Xcalibur, Tracefinder, NIST Lib) use 1000 

    Returns
    -------
    match_factor : float
        Quantification of similarity between both spectra
        Essentially a dot product between two spectra
            Need to ensure spectra are evalauted on correct vector space

    """
    # construct common basis
    common_basis = set(sample_mz + ref_mz)
    
    # mape to desired basis
    if mode == 'direct':
        basis = list(set(sample_mz + ref_mz))
    elif mode == 'forward':
        basis = sample_mz
        for signal in range(len(sample_sig)):
            if sample_sig[signal] == 0:
                ref_sig[signal] = 0
    elif mode == 'reverse':
        basis = ref_mz
        for signal in range(len(ref_sig)):
            if ref_sig[signal] == 0:
                sample_sig[signal] = 0
    else:
        print('mode is defined as direct, forwaard, or reverse')

    # map both spectra to common basis
    mapped_sample = vector_map(sample_sig, sample_mz, basis)
    mapped_ref = vector_map(ref_sig, ref_mz, basis)
    
    # normalize basis
    norm_sample = vector_norm(mapped_sample)
    norm_ref = vector_norm(mapped_ref)
    
    return np.dot(norm_sample,norm_ref) * mult

