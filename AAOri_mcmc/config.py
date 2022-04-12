#config.py
import numpy as np

topdir = '/Users/khkim/Work2/MyDiskModelPy/'
obsdatatxt = topdir+'AAOri_emc/obsdata/AAOri_K4_spec_sample4emcee.txt'
odata = np.loadtxt(obsdatatxt, skiprows= 0 + 1+ 2, usecols=(0,1,2,3,4))
wobs = odata[:,0]
fnobs = odata[:,1] #fn
fnerr = odata[:,2]
fobs = odata[:,3] #nfn
ferr = odata[:,4]


x=wobs
y=fobs
yerr = ferr


# select testing group
# Set up the backend
#testgroup = 'FDd2'
#subdir='AAOri_emc/test_FDd2_rt_emcee/'
#h5filename = topdir+subdir+"FDd2.h5"
#fixed_param_file = topdir+'AAOri_emc/default_inps/problem_params_fixed_FD_2dusts.inp'
#dustkap = ('dr_MW55, WD-mm')
#niter = 2
#
#testgroup = 'TDd2'
#subdir='TDd2_rt_emcee/'
#h5filename = topdir+subdir+"TDd2.h5"
#dustkap = ('dr_MW55, WD-mm')
#niter = 100
#
#testgroup = 'FDd3'
#subdir='FDd3_rt_emcee/'
#h5filename = topdir+subdir+"FDd3.h5"
#dustkap = ('dr_MW55, mgalsi2_01-1um, WD-mm')
#fixed_param_file = topdir+'AAOri_emc/default_inps/problem_params_fixed_FD_2dusts.inp'
#niter = 100
#
testgroup = 'TDd3'
subdir='AAOri_emc/test-multi/'
h5filename = topdir+subdir+"TDd3.h5"
dustkap = ('dr_MW55, mgalsi2_01-1um, WD-mm')
fixed_param_file = topdir+'AAOri_emc/default_inps/problem_params_fixed_TD_3dusts.inp'
niter = 2
#global current_params4run

if testgroup == 'FDd2':
    initial_params4run = {'hrdisk': 0.12,\
        'rdisk': 600, \
            'plh': 0.143,\
                'inclination': 40,\
                    'mfrac_d1': 0.5,\
                    }
    labels = [r'h$_r$', r'r$_{out}$', r'$\psi$', r'inc', r'f$_{small}$']
    for k, v in initial_params4run.items():
        exec('%s = %s' % (k,v))
    #theta_i = [hrdisk, rdisk, plh, inclination, mfrac_d1, mfrac_d2]
    theta_i = [hrdisk, rdisk, plh, inclination, mfrac_d1]

    ndim = len(theta_i)
    nwalkers = ndim*2 + 1
    #pos = theta_i + 1e-3 * np.random.randn(nwalkers, ndim)

    po0 = theta_i[0] + 1e-2 * np.random.randn(nwalkers)
    po1 = theta_i[1] + 5e1 * np.random.randn(nwalkers)
    po2 = theta_i[2] + 1e-2 * np.random.randn(nwalkers)
    po3 = theta_i[3] + 1e1 * np.random.randn(nwalkers)
    po4 = theta_i[4] + 1e-1 * np.random.randn(nwalkers)

    #pos = np.array((nwalkers,ndim)
    posa = []
    for j in range(nwalkers):
        posn=[po0[j], po1[j], po2[j], po3[j], po4[j]]
        posa.append(posn)
    pos = np.array(posa)

#Full disk and 3 dusts test를 위한 initial parameter dictionary
if testgroup == 'FDd3':
    initial_params4run = {'hrdisk': 0.12,\
        'rdisk': 600, \
            'plh': 0.143,\
                'inclination': 40,\
                    'mfrac_lg': 0.5,\
                    'sm1_frac': 0.1,\
                }
    labels = [r'h$_r$', r'r$_{out}$', r'$\psi$', r'inc', r'f$_{large}$', r'sm$_{1}$/sm$_{Total}$']
    for k, v in initial_params4run.items():
        exec('%s = %s' % (k,v))

    theta_i = [hrdisk, rdisk, plh, inclination, mfrac_lg, sm1_frac]

    ndim = len(theta_i)
    nwalkers = ndim*2 + 1
    #pos = theta_i + 1e-3 * np.random.randn(nwalkers, ndim)

    po0 = theta_i[0] + 1e-2 * np.random.randn(nwalkers)
    po1 = theta_i[1] + 5e1 * np.random.randn(nwalkers)
    po2 = theta_i[2] + 1e-2 * np.random.randn(nwalkers)
    po3 = theta_i[3] + 1e1 * np.random.randn(nwalkers)
    po4 = theta_i[4] + 1e-1 * np.random.randn(nwalkers)
    po5 = theta_i[5] + 1e-1 * np.random.randn(nwalkers)

    #pos = np.array((nwalkers,ndim)
    posa = []
    for j in range(nwalkers):
        posn=[po0[j], po1[j], po2[j], po3[j], po4[j], po5[j]]
        posa.append(posn)
    pos = np.array(posa)

#Full disk and 3 dusts test를 위한 initial parameter dictionary
if testgroup == 'TDd3':
    initial_params4run = {'hrdisk': 0.12,\
                'rdisk': 600, \
                'plh': 0.143,\
                'inclination': 40,\
                'mfrac_lg': 0.5,\
                'sm1_frac': 0.5,\
                'gap_rin':10, \
                'gap_rout':60, \
                'p_in':-1, \
                'eta_in_sm1': 1e-1, \
                'eta_in_sm2': 1e-1, \
                'eta_in_lg': 1e-1, \
                'eta_gap_sm1': 1e-6, \
                'eta_gap_sm2': 1e-6, \
                'eta_gap_lg': 1e-6
            }
    labels = [r'h$_r$', r'r$_{out}$', r'$\psi$', r'inc', r'f$_{large}$', r'sm$_{1}$/sm$_{Total}$', r'r$_{gap\_in}$', r'r$_{gap\_out}$', r'p$_{in}$',r'$\eta_{in\_sm1}$', r'$\eta_{in\_sm2}$', r'$\eta_{in\_lg}$', r'$\eta_{gap\_sm1}$', r'$\eta_{gap\_sm2}$', r'$\eta_{gap\_lg}$']
    theta_i = np.array([])
    for k, v in initial_params4run.items():
        exec('%s = %s' % (k,v))
    theta_i = [hrdisk, rdisk, plh, inclination, mfrac_lg, sm1_frac, gap_rin, gap_rout, p_in, eta_in_sm1, eta_in_sm2, eta_in_lg, eta_gap_sm1, eta_gap_sm2, eta_gap_lg]

    ndim = len(theta_i)
    nwalkers = ndim*2 + 1
    #pos = theta_i + 1e-3 * np.random.randn(nwalkers, ndim)

    po0 = theta_i[0] + 1e-2 * np.random.randn(nwalkers)
    po1 = theta_i[1] + 5e1 * np.random.randn(nwalkers)
    po2 = theta_i[2] + 1e-2 * np.random.randn(nwalkers)
    po3 = theta_i[3] + 1e1 * np.random.randn(nwalkers)
    po4 = theta_i[4] + 1e-1 * np.random.randn(nwalkers)
    po5 = theta_i[5] + 1e-1 * np.random.randn(nwalkers)
    po6 = theta_i[6] + 1 * np.random.randn(nwalkers)
    po7 = theta_i[7] + 5 * np.random.randn(nwalkers)
    po8 = theta_i[8] + 5e-1 * np.random.randn(nwalkers)
    po9 = theta_i[9] + 5e-5 * np.random.randn(nwalkers)
    po10 = theta_i[10] + 5e-5 * np.random.randn(nwalkers)
    po11 = theta_i[11] + 5e-5 * np.random.randn(nwalkers)
    po12 = theta_i[12] + 1e-5 * np.random.randn(nwalkers)
    po13 = theta_i[13] + 1e-5 * np.random.randn(nwalkers)
    po14 = theta_i[14] + 1e-5 * np.random.randn(nwalkers)

    #pos = np.array((nwalkers,ndim)
    posa = []
    for j in range(nwalkers):
        posn=[po0[j], po1[j], po2[j], po3[j], po4[j], po5[j], po6[j], po7[j], po8[j], po9[j], po10[j], po11[j], po12[j], po13[j], po14[j]]
        posa.append(posn)
    pos = np.array(posa)
