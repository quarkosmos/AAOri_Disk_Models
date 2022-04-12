
import numpy as np
import matplotlib.pyplot as plt
import emcee
import corner

from astropy.io import ascii
import time
import os
import re
import os.path
import radmc3dPy
from radmc3dPy import *
import radmc3dPy.natconst as nc

import aaori_emc_radmc3dtool_global2 as rtool
import radmc3d_emcee_utilities_global2 as emcutil

import config

log_probability_global = emcutil.log_probability_global

# Don't forget to clear it in case the file already exists
#h5filename = topdir+"AAOri_emc/test1_rt_emcee/FD_test.h5"
if os.path.exists(config.h5filename):
    os.remove(config.h5filename)
backend = emcee.backends.HDFBackend(config.h5filename)
backend.reset(config.nwalkers, config.ndim)

sampler = emcee.EnsembleSampler(config.nwalkers, config.ndim, log_probability_global, backend=backend)
start = time.time()
sampler.run_mcmc(config.pos, config.niter)
end = time.time()
serial_data_time = end - start
print("Serial took {0:.1f} seconds".format(serial_data_time))

#\uc5ec\uae30\uae4c\uc9c0 \ub3cc\ub9ac\uace0, \ub098\uc11c.. \uc774\ud6c4\uc5d0 ndim vs step \uadf8\ub9bc\uc744 \ud655\uc778\ud55c \ud6c4 discard step\uc744 \uacb0\uc815\uc744 \ud574\uc57c \ud558\uaca0\ub2e4.
fig, axes = plt.subplots(config.ndim, figsize=(10, 15), sharex=True)
samples = sampler.get_chain()

for i in range(config.ndim):
    ax = axes[i]
    ax.plot(samples[:, :, i], "k", alpha=0.3)
    ax.set_xlim(0, len(samples))
    ax.set_ylabel(config.labels[i])
    ax.yaxis.set_label_coords(-0.1, 0.5)

axes[-1].set_xlabel("step number");


# \uc704\uc5d0\uc11c discard step \uc815\ud55c \ud6c4\uc5d0 thin, discard \uac12 \uacb0\uc815\ud558\uc5ec \ub2e4\uc74c \ud568\uc218\uc5d0 \ubc18\uc601\ud558\uae30.
flat_samples = sampler.get_chain(flat=True)
print(flat_samples.shape)

fig = corner.corner(flat_samples,show_titles=True,labels=config.labels,plot_datapoints=True,quantiles=[0.16, 0.5, 0.84])


##### Run MCMC from the last ended sampler by calling HDF5 file
new_backend = emcee.backends.HDFBackend(config.h5filename)
print("Initial size: {0}".format(new_backend.iteration))
new_sampler = emcee.EnsembleSampler(config.nwalkers, config.ndim, log_probability_global, backend=new_backend)
start = time.time()
new_sampler.run_mcmc(None, 48)
end = time.time()
serial_data_time = end - start
print("Serial took {0:.1f} seconds".format(serial_data_time))
print("Final size: {0}".format(new_backend.iteration))
