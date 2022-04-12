"""
This is a code including several Functions to facilitate batch running models

testSEDcsq() : calculate chisquare between data and model fluxes and return csq and redcsq
# genModelDir() : make directory
genProbParmInp() : read input parameter table sets and make a directory for each model setup. then make a problem_params.inp file in each model directory
batch_run_radmc3d() : go to each model directory, and run radmc3dPy
genDefaultFigures(): SED, surface density, density distribution


"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as plb
from astropy.io import ascii
from scipy.stats import chisquare
import os
import csv
import re
import os.path
import uuid

import radmc3dPy
from radmc3dPy import *
import radmc3dPy.natconst as nc

from matplotlib.backends.backend_pdf import PdfPages


import config

def params2dict(theta):
    """
    basically emcee will generate parameter values depending on the parameters defined.
    starting from initial_params4run in the preamp of the main code, the current_params4run will be updated with variable parameters
    generated for a specific nwalk and nstep.
    Then, current_params4run dic will be used for current_params2write dictionry.

    """
    #global params2write
    #global current_params4run
    current_params4run = config.initial_params4run.copy()
    if config.testgroup == 'FDd2':
        a = list(current_params4run.keys())
        for i in range(len(a)):
            current_params4run[a[i]] = theta[i]

        params2write = current_params4run.copy()
        params2write['mfrac'] =[current_params4run.get('mfrac_d1'), 1-current_params4run.get('mfrac_d1')]
        del params2write['mfrac_d1']


    if config.testgroup == 'FDd3':
        a = list(current_params4run.keys())
        for i in range(len(a)):
            current_params4run[a[i]] = theta[i]

        params2write = current_params4run.copy()
        smfrac = 1-current_params4run.get('mfrac_lg')
        mfrac_sm1 = current_params4run.get('sm1_frac')*smfrac
        mfrac_sm2 = (1-current_params4run.get('sm1_frac'))*smfrac
        params2write['mfrac'] =[mfrac_sm1, mfrac_sm2, current_params4run.get('mfrac_lg')]
        del params2write['mfrac_lg']
        del params2write['sm1_frac']

    if config.testgroup == 'TDd3':
        a = list(current_params4run.keys())
        for i in range(len(a)):
            current_params4run[a[i]] = theta[i]

        params2write = current_params4run.copy()
        smfrac = 1-current_params4run.get('mfrac_lg')
        mfrac_sm1 = current_params4run.get('sm1_frac')*smfrac
        mfrac_sm2 = (1-current_params4run.get('sm1_frac'))*smfrac
        hrdisk_in = current_params4run.get('hrdisk')*(0.01)**current_params4run.get('plh')
        params2write['mfrac'] =[mfrac_sm1, mfrac_sm2, current_params4run.get('mfrac_lg')]
        del params2write['mfrac_lg']
        del params2write['sm1_frac']
        params2write['ind_drfact_idust'] = [current_params4run.get('eta_in_sm1'), current_params4run.get('eta_in_sm2'), current_params4run.get('eta_in_lg')]
        params2write['gap_drfact_idust'] = [current_params4run.get('eta_gap_sm1'), current_params4run.get('eta_gap_sm2'), current_params4run.get('eta_gap_lg')]
        params2write['ind_drpower'] = [current_params4run.get('p_in'), current_params4run.get('p_in'), current_params4run.get('p_in')]
        params2write['hrdisk_out'] = current_params4run.get('hrdisk')
        params2write['hrdisk_in'] = hrdisk_in
        del params2write['eta_in_sm1']
        del params2write['eta_in_sm2']
        del params2write['eta_in_lg']
        del params2write['eta_gap_sm1']
        del params2write['eta_gap_sm2']
        del params2write['eta_gap_lg']
        del params2write['p_in']
        print(params2write)

    return params2write



def write_paraminp(p2w, dirpath):
    """
    this function is to make a directory to run a model through radmc3dPy.
    It receive parameters generated from emcee EnsembleSampler and make a problem_params.inp file in the model run directory.
    Then, it will returen the directory name, so let the program know where other files (opacity, density, wavelength etc..)
    should be written to run a disk model via radmc3dPy

    the input dic_param : params2write
    """
    #rundirname = uuid.uuid4().hex
    #dirpath = config.topdir+config.subdir
    ndir = len(next(os.walk(dirpath))[1])
    num_credir = ndir+1
    destdir = dirpath+'m_'+str(num_credir)+'/'
    #destdir = dirpath+rundirname+'/'
    destparamfile = destdir+'problem_params.inp'
    os.system('mkdir '+destdir)
    os.system('cp '+config.fixed_param_file+' '+destparamfile)
    d_key = list(p2w.keys())
    d_val = list(p2w.values())
    f = open(destparamfile, "a")
    for  i in range(len(d_key)):
        if d_key[i] == 'rdisk':
            s = d_key[i]+'                   = '+str(d_val[i])+'*au'
        elif d_key[i] == 'gap_rin':
            s = d_key[i]+'                   = ['+str(d_val[i])+'*au]'
        elif d_key[i] == 'gap_rout':
            s = d_key[i]+'                   = ['+str(d_val[i])+'*au]'
        else:
            s = d_key[i]+'                   = '+str(d_val[i])
        f.write(s)
        f.write('\n')

    if config.testgroup=='TDd3':
        f.write('plh_in                     ='+str(p2w.get('plh')))
        f.write('\n')
        f.write('plh_out                     ='+str(p2w.get('plh')))
        f.write('\n')

    f.close()

    return destdir



#def delete_modelrundir(destdir):





def run_radmc3d(rt_run_dir, p2w):
    incangl  = p2w.get('inclination')
    os.chdir(rt_run_dir)
    dustelements = re.split(", ", config.dustkap)
    for i in range(0, len(dustelements)):
        dustopacfilename = 'dustkappa_'+dustelements[i]+'.inp'
        os.system('cp -v '+config.topdir+'AAOri_emc/default_inps/'+dustopacfilename+' .')


    rt_model = radmc3dPy.setup.problemSetupDust(model='ppdisk_aaori', verbose='False', binary=False)
    if os.path.exists('radmc3d.inp'):
        with open('radmc3d.inp', 'r+') as file:
            for line in file:
                pass
            file.write('mc_scat_maxtauabs = 15.d0')
    os.system('cp -v '+config.topdir+'AAOri_emc/default_inps/stars.inp.forFDtest.inp ./stars.inp')
    os.system('cp -v '+config.topdir+'AAOri_emc/default_inps/wavelength_micron.inp.forFDtest.inp ./wavelength_micron.inp')
    os.system('radmc3d mctherm countwrite 100000')
    if incangl == 0:
        os.system('radmc3d sed countwrite 100000')
    else:
        os.system('radmc3d sed incl '+str(incangl)+' countwrite 100000')

    modeltxt = 'spectrum.out'
    mdata = np.loadtxt(modeltxt, skiprows= 0 + 1+ 2)
    wmodel = mdata[1:,0] #do not use the first wavelength & its flux
    fmodel = mdata[1:,1]*(3e14)/wmodel/414/414

    genDefaultFigures(rt_run_dir)

    os.system('rm dustkappa_*.inp')
    os.system('rm dustopac.inp')
    os.system('rm stars.inp')
    os.system('rm wavelength_micron.inp')


    return fmodel



def genDefaultFigures(rt_run_dir):

    from astropy.modeling.models import BlackBody
    from astropy import units as u
    from astropy.visualization import quantity_support

    #obsdatatxt = topdir+'AAOri_emc/obsdata/AAOri_K4_spec_sample4emcee.txt'
    obsdatatxt = config.topdir+'AAOri_emc/obsdata/AAOri_K4_spec_sample4emcee.txt'
    odata = np.loadtxt(obsdatatxt, skiprows= 0 + 1+ 2, usecols=(0,1,2,3,4))
    wobs = odata[:,0]
    fnobs = odata[:,1] #fn
    fnerr = odata[:,2]
    fobs = odata[:,3] #nfn
    ferr = odata[:,4]

### for the blackbody photosphere
    jidx = [i for i, value in enumerate(wobs) if value == 1.25]
    hidx = [i for i, value in enumerate(wobs) if value == 1.65]
    fobs_J = fobs[jidx]
    fobs_H = fobs[hidx]
    lfl_H = fobs_H

    bb = BlackBody(temperature=4600.2*u.K)
    wav = np.arange(0.1, 2000, 0.1) * u.micron
    flux = bb(wav)
    lfl = flux*wav.to(u.Hz, equivalencies=u.spectral())

    flux_J = bb(1.25*u.micron)
    wJ = 1.25 * u.micron
    lfl_J = flux_J * wJ.to(u.Hz, equivalencies=u.spectral())


    bbscale = fobs_J*(u.erg/u.cm/u.cm/u.s)/lfl_J
    lfl_norm = lfl * bbscale
############

    with PdfPages(rt_run_dir+'sed_sigma_rho.pdf') as pdf:

        dustelements = re.split(", ", config.dustkap)

        os.chdir(rt_run_dir)

            # 2 x 2 subplot frame
        fig = plt.figure(figsize=(11,8))
        f1 = fig.add_subplot(2,2,1)
        f2 = fig.add_subplot(2,2,2)
        f3 = fig.add_subplot(2,2,3)
        f4 = fig.add_subplot(2,2,4)


            # plot SED
        modeltxt = 'spectrum.out'
        mdata = np.loadtxt(modeltxt, skiprows= 0 + 1+ 2)
        wmodel = mdata[1:,0]
        fmodel = mdata[1:,1]*(3e14)/wmodel/414/414
        #fig1
        f1.plot(wav, lfl_norm, 'k', linewidth=0.5,  label='photosphere (blackbody)')
        f1.plot(wmodel,fmodel, 'r')
        f1.scatter(wmodel, fmodel, c='y', s=10, marker="s")
        f1.scatter(wobs,fobs, s=10)
        f1.errorbar(wobs,fobs, yerr=ferr, linestyle='none', ecolor='black', elinewidth=0.5)
        f1.set_xscale('log')
        f1.set_yscale('log')
        f1.set_xlabel('wavelength '+ r'($\mu$m)')
        f1.set_ylabel(r'$\nu$$f_\nu$ (erg/s/cm$^{2}$)')
        f1.set_xlim([0.1, 2000])
        f1.set_ylim([8e-15, 2e-8])
        #f1.axis([0.1, 2000, 8e-15, 2e-8])
        f1.set_title('SED')#dircode)
        #f1.legend([r'$\chi^2_{all}$:'+(str(csvalue[0])[:4])])
        #f1.text(0.2, 1e-12, r'$\chi^2_{all}$:'+(str(csvalue[0])[:4]) , fontsize=9)
        #f1.text(0.2, 2.5e-13, r'$\chi^2_{10}$:'+(str(csvalue[2])[:4]) , fontsize=9)
        #f1.text(0.2, 8e-14, r'$\chi^2_{other}$:'+(str(csvalue[4])[:4]) , fontsize=9)

        f3.plot(wmodel,fmodel, 'r')
        f3.scatter(wmodel, fmodel, c='y', s=10, marker="s")
        f3.scatter(wobs,fobs, s=10)
        f3.errorbar(wobs,fobs, yerr=ferr, linestyle='none', ecolor='black', elinewidth=0.5)
        f3.set_xscale('linear')
        f3.set_yscale('log')
        f3.set_xlabel('wavelength '+ r'($\mu$m)')
        f3.set_ylabel(r'$\nu$$f_\nu$ (erg/s/cm$^{2}$)')
        f3.set_xlim([5, 15])
        f3.set_ylim([8e-11, 9e-10])
        #f3.axis([7, 20, 8e-11, 2e-9])
        f3.set_title('10'+ r'($\mu$m)'+' zoom up')


        # plot surface density
        data = analyze.readData(ddens=True, binary=False)


        #---- surface density vs radius
        # plot dust surface density profile calculated from getSigmaDust method
        s=data.getSigmaDust()
        rtemp=np.array(s[0])
        r=rtemp[1:]
        ytemp = s[1][:,0]
        #axisrange=[0.1, 410, 1e-10, 0.9]
        dustlegend = ['all']
        f2.plot(r, ytemp, 'b')
        f2.set_xscale('log')
        f2.set_yscale('log')
        f2.set_xlabel('r (AU)')
        f2.set_ylabel(r'$\Sigma$ (g/cm$^{2}$)')
        #plt.axis(axisrange)
        #plt.ylim(9e-11, 1e1)
        f2.set_title('surface density')
        for i in range(0, len(dustelements)):
            s=data.getSigmaDust(idust=i)
            rtemp=np.array(s[0])
            r=rtemp[1:]
            ytemp = s[1][:,0]
            f2.plot(r, ytemp)
            dustlegend.append(dustelements[i])
            # if i == 0:
            #     dustlegend.append(dustelements[i][2:-1])
            # elif i == len(dustelements)-1:
            #     dustlegend.append(dustelements[i][1:-2])
            # else:
            #     dustlegend.append(dustelements[i][1:-1])
        f2.legend(dustlegend)


        # plot density contours
        vmin=(1e-33)
        vmax=(1e-14)

        analyze.plotSlice2D(data, var='ddens', plane='xy', log=True, linunit='au', cmap=plb.cm.jet, vmin=vmin, vmax=vmax)
        data.getTau(wav=10)
        analyze.plotSlice2D(data, var='taux', plane='xy', log=True, linunit='au', contours=True, clev=[1.0], clcol='w', cllabel=True)
        data.getTau(wav=850)
        analyze.plotSlice2D(data, var='taux', plane='xy', log=True, linunit='au', contours=True, clev=[1.0], clcol='g', cllabel=True)
        plb.xscale('log')
        #Tdata = analyze.readData(dtemp=True)
        #analyze.plotSlice2D(Tdata, var='dtemp', plane='xy', ispec=0, log=True, linunit='au', contours=True, clev=[1400.0], clcol='k', cllabel=True)
        f4.set_title('dust density (all)')


        #fig.suptitle(dircode)
        fig.subplots_adjust(hspace=0.35, wspace=0.35)
        pdf.savefig(orientation='landscape')
        plt.close()
