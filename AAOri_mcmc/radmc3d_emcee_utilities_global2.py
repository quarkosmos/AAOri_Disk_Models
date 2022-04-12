import numpy as np
import matplotlib.pyplot as plt
import emcee
import corner


from astropy.io import ascii
import os
import re
import os.path
import radmc3dPy
from radmc3dPy import *
import radmc3dPy.natconst as nc

import aaori_emc_radmc3dtool_global2 as rtool

import config



def log_likelihood(theta, rt_run_dir, p2w):
    ymodel = rtool.run_radmc3d(rt_run_dir, p2w)
    return -0.5 * np.sum(((config.y - ymodel)/config.yerr) ** 2)

# The function for specifying input paramters: prior
# for example
def log_prior(theta):
    if config.testgroup == 'FDd2':
        hrdisk, rdisk, plh, inclination, mfrac_d1 = theta
        #totalmfrac = mfrac_d1 + mfrac_d2
        if 0.09 <= hrdisk <= 0.2 and \
            200 <= rdisk <= 1000 and \
            1/20 <= plh <= 3/7 and \
            0 <= inclination <= 80 and \
            0.1 <= mfrac_d1 <= 0.9:
            return 0.0
        else:
            return -np.inf

    if config.testgroup == 'FDd3':
        hrdisk, rdisk, plh, inclination, mfrac_lg, sm1_frac = theta
            #totalmfrac = mfrac_d1 + mfrac_d2
        if 0.09 <= hrdisk <= 0.2 and \
            200 <= rdisk <= 1000 and \
            1/20 <= plh <= 3/7 and \
            0 <= inclination <= 80 and \
            0.1 <= mfrac_lg <= 0.9 and \
            0.1 <= sm1_frac <= 0.9:
            return 0.0
        else:
            return -np.inf

    if config.testgroup == 'TDd3':
        hrdisk, rdisk, plh, inclination, mfrac_lg, sm1_frac, gap_rin, gap_rout, p_in, eta_in_sm1, eta_in_sm2, eta_in_lg, eta_gap_sm1, eta_gap_sm2, eta_gap_lg = theta
        if 0.09 <= hrdisk <= 0.2 and \
            200 <= rdisk <= 1000 and \
            1/20 <= plh <= 3/7 and \
            0 <= inclination <= 80 and \
            0.1 <= mfrac_lg <= 0.9 and \
            0.1 <= sm1_frac <= 0.9 and \
            0.25 <= gap_rin <=80 and \
            0.25 <= gap_rout <=80 and \
            gap_rin < gap_rout and \
            -3 <= p_in <= 3 and \
            1e-2 <= eta_in_sm1 <= 1 and\
            1e-2 <= eta_in_sm2 <= 1 and \
            1e-2 <= eta_in_lg <= 1 and \
            1e-8 <= eta_gap_sm1 <= 1e-3 and\
            1e-8 <= eta_gap_sm2 <= 1e-3 and \
            1e-8 <= eta_gap_lg <= 1e-3:
            return 0.0
        else:
            return -np.inf



def log_probability_global(theta):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    else:
        #global rt_run_dir
        p2w=rtool.params2dict(theta)
        rt_run_dir = rtool.write_paraminp(p2w, config.topdir+config.subdir)
        return lp + log_likelihood(theta, rt_run_dir, p2w)
