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

import radmc3dPy
from radmc3dPy import *
import radmc3dPy.natconst as nc

from matplotlib.backends.backend_pdf import PdfPages



def testSEDcsq(specdir, varpar):  #chisq 게산에서 첫번째 자료는 제외 (model값이 항상 너무 뚝 떨어)
    obsdatatxt_chisq = '/Users/khkim/Work2/MyDiskModelPy/AAOri_emc/obsdata/AAOri_K4_spec_sample4emcee.txt'
    #obsdatatxt_chisq = '/Users/khkim/Work2/MyDiskModelPy/AAOri_emc/obsdata/AAOri_K4_spec_sample4emcee_forfullsed.txt'
    odata_chisq = np.loadtxt(obsdatatxt_chisq, skiprows= 0 + 1+ 2, usecols=(0,1,2,3,4))
    wobs_chisq = odata_chisq[1:,0]
    fnobs_chisq = odata_chisq[1:,1] #fn
    fobs_chisq = odata_chisq[1:,3] #nfn
    ferr_chisq = odata_chisq[1:,4]

    modeltxt = specdir+'/spectrum.out'
    mdata = np.loadtxt(modeltxt, skiprows= 4)
    wmodel = mdata[1:,0]
    fmodel = mdata[1:,1]*(3e14)/wmodel/414/414
    fnmodel = fmodel*wmodel/(3e14)*(1e23)

    #chisquare(fnobs_chisq, fnmodel, ddof=[17])
    #chisquare(fnobs_chisq, fnmodel)

    csq_ind = (fnobs_chisq-fnmodel)**2/fnmodel
    csq_all = sum(csq_ind)
    Nall = len(csq_ind)
    #print(csq)
    print(len(csq_ind))
    redcsq= csq_all/(len(csq_ind)-varpar)

    #--------------
    #find index for wavelength between 8 and 12 microns
    mk = [i for i,e in enumerate(wmodel) if e >= 8. and e<=12.0]
    #print(len(mk))
    #wmk10 = wmodel[mk[0]:mk[-1]]
    #fnm10 = fnmodel[mk[0]:mk[-1]]
    fnm10 = [fnmodel[i] for i in mk]
    #print(wmk10)
    print(len(fnm10))

    ok = [i for i,e in enumerate(wobs_chisq) if e >= 8. and e<=12.0]
    #wo10 = wobs_chisq[ok[0]:ok[-1]]
    #fno10 = fnobs_chisq[ok[0]:ok[-1]]
    fno10 = [fnobs_chisq[i] for i in ok]

    chisq_y10_ind = (np.array(fno10)-np.array(fnm10))**2/np.array(fnm10)
    chisq_y10 = sum(chisq_y10_ind)
    N_y10 = len(chisq_y10_ind)

    #--------------
    #find index for wavelength other than between 8 and 12 microns
    mt = [i for i,e in enumerate(wmodel) if e < 8. or e > 12.]
    #wmt10 = wmodel[mt[0]:mk[-1]]
    fnm10n = [fnmodel[i] for i in mt]
    print(len(fnm10n))

    ot = [i for i,e in enumerate(wobs_chisq) if e < 8. or e > 12.]
    #wo10 = wobs_chisq[ok[0]:ok[-1]]
    fno10n = [fnobs_chisq[i] for i in ot]

    chisq_n10_ind = (np.array(fno10n)-np.array(fnm10n))**2/np.array(fnm10n)
    chisq_n10 = sum(chisq_n10_ind)
    N_n10 = len(chisq_n10_ind)




#chisq_all   Nall    chisq_y10   N_y10    chisq_n10    N_n10
    return csq_all, Nall, chisq_y10, N_y10, chisq_n10, N_n10


#def genModelDir(dirname):
#    os.system('mkdir '+dirname)


def genProbParmInp(parm_table, dirname):
    #parm_table = "/Users/khkim/Work2/MyDiskModelPy/AAOri/param_table.inp/param_table_gap.inp"
    #dirname = "/Users/khkim/Work2/MyDiskModelPy/AAOri/rtmodels/"
    parm_table = parm_table

    num_cols = 0

    with open(parm_table, 'r') as pt:
        reader = csv.reader(pt, delimiter=' ', skipinitialspace=True)
        first_row = next(reader)
        num_cols = len(first_row)

    with open(parm_table, 'r') as pt:
        #reader = csv.reader(pt, delimiter=' ', skipinitialspace=True)
        linecount = 0
        for l in pt:
            if not l.strip().startswith('#'):
                linecount += 1

        print(linecount, num_cols)

    for icolumn in range(2, num_cols):
        gparm = np.loadtxt(parm_table, dtype=str, usecols=(1, icolumn))
        dircode = gparm[0, 1]
        destdir = dirname+dircode
        destparmfile = destdir+'/problem_params.inp'
        os.system('mkdir '+destdir)
        #os.system('cp /Users/khkim/Work2/MyDiskModelPy/AAOri/param_table.inp/problem_params_fixed.inp '+destparmfile)
        os.system('cp /Users/khkim/Work2/MyDiskModelPy/AAOri_emc/default_inps/problem_params_fixed_TDFD_emcee.inp '+destparmfile)
        f = open(destparmfile, "a")
        for  l in range(1, linecount):
            s = gparm[l,0]+'                   = '+gparm[l,1]
            f.write(s)
            f.write('\n')
        f.close()


def batch_run_radmc3d(parm_table, dirname):
        #parm_table = "/Users/khkim/Work2/MyDiskModelPy/AAOri/param_table.inp/param_table_gap.inp"
        #dirname = "/Users/khkim/Work2/MyDiskModelPy/AAOri/rtmodels/"
        parm_table = parm_table
        num_cols = 0
        with open(parm_table, 'r') as pt:
            reader = csv.reader(pt, delimiter=' ', skipinitialspace=True)
            first_row = next(reader)
            num_cols = len(first_row)

        with open(parm_table, 'r') as pt:
            linecount = 0
            #reader = csv.reader(pt, delimiter=' ', skipinitialspace=True)
            for l in pt:
                if not l.strip().startswith('#'):
                    linecount += 1
            print(linecount)



        # cp dust opacity input file to the destination directory
        for icolumn in range(2, num_cols):
            gparm = np.loadtxt(parm_table, dtype=str, usecols=(1, icolumn))
            dircode = gparm[0, 1]
            dustkap = gparm[1,1]
            destdir = dirname+dircode
            incangl = gparm[linecount-1, 1]
            #print('incangle=', incangl)
            os.chdir(destdir)
            dustelements = re.split(",", dustkap)
            for i in range(0, len(dustelements)):
                if i == 0:
                    dustopacfilename = 'dustkappa_'+dustelements[i][2:-1]+'.inp'
                elif i == len(dustelements)-1:
                    dustopacfilename = 'dustkappa_'+dustelements[i][1:-2]+'.inp'
                else:
                    dustopacfilename = 'dustkappa_'+dustelements[i][1:-1]+'.inp'
                os.system('cp -v /Users/khkim/Work2/MyDiskModelPy/AAOri_emc/default_inps/'+dustopacfilename+' .')

            # run a model
            model = radmc3dPy.setup.problemSetupDust(model='ppdisk_aaori', verbose='False', binary=False)

            # Copy the user defined input files for the wavelength ranges and stellar property
            #os.system('cp -v /Users/khkim/Work2/MyDiskModelPy/AAOri_emc/default_inps/stars.inp.forfullsed.inp ./stars.inp')
            #os.system('cp -v /Users/khkim/Work2/MyDiskModelPy/AAOri_emc/default_inps/wavelength_micron.inp.forfullsed.inp ./wavelength_micron.inp')
            os.system('cp -v /Users/khkim/Work2/MyDiskModelPy/AAOri_emc/default_inps/stars.inp.forFDtest.inp ./stars.inp')
            os.system('cp -v /Users/khkim/Work2/MyDiskModelPy/AAOri_emc/default_inps/wavelength_micron.inp.forFDtest.inp ./wavelength_micron.inp')


            # Calculate the dust temperature
            os.system('radmc3d mctherm')

            # ------------------------------------------------------------------------------------------------------------
            # Run SED
            # ------------------------------------------------------------------------------------------------------------
            if float(incangl) == 0.0:
                os.system('radmc3d sed')
            else:
                os.system('radmc3d sed incl '+incangl)



def genDefaultFigures(parm_table, dirname, outputdir, filename, outtxtfile):
    # plot SED and compare obs and model data. plot surface density profile (overplot each dust speceise)
    # plot density distribution
    #parm_table = "/Users/khkim/Work2/MyDiskModelPy/AAOri/param_table.inp/param_table_gap.inp"
    #dirname = "/Users/khkim/Work2/MyDiskModelPy/AAOri/rtmodels/"
    #outputdir = "/Users/khkim/Work2/MyDiskModelPy/AAOri/rtmodels/out_fig_csq/"
    #outtxtfile = "chisqtestresults_g.txt" or "chisqtestresults_g.txt"

    parm_table = parm_table
    with open(parm_table, 'r') as pt:
        reader = csv.reader(pt, delimiter=' ', skipinitialspace=True)
        first_row = next(reader)
        num_cols = len(first_row)
        #print(num_cols)

    with PdfPages(outputdir+filename+'.pdf') as pdf:
        # the observational data
        obsdatatxt = '/Users/khkim/Work2/MyDiskModelPy/AAOri_emc/obsdata/AAOri_K4_spec_sample4emcee.txt'
        odata = np.loadtxt(obsdatatxt, skiprows= 0 + 1+ 2, usecols=(0,3,4))
        wobs = odata[:,0]
        fobs = odata[:,1]
        ferr = odata[:,2]

        if not os.path.isfile(outputdir+outtxtfile):
            setupOutFiles(outputdir+outtxtfile)
        #else:
        #print("here")

        # for loop to get data files in the models directories
        for icolumn in range(2, num_cols):
            gparm = np.loadtxt(parm_table, dtype=str, usecols=(1, icolumn))
            dircode = gparm[0, 1]
            dustkap = gparm[1,1]
            dustelements = re.split(",", dustkap)

            destdir = dirname+dircode
            os.chdir(destdir)
            print(destdir)
            # 2 x 2 subplot frame
            fig = plt.figure(figsize=(11,8))
            f1 = fig.add_subplot(2,2,1)
            f2 = fig.add_subplot(2,2,2)
            f3 = fig.add_subplot(2,2,3)
            f4 = fig.add_subplot(2,2,4)


            # plot SED
            modeltxt = 'spectrum.out'
            mdata = np.loadtxt(modeltxt, skiprows= 0 + 1+ 2)
            wmodel = mdata[:,0]
            fmodel = mdata[:,1]*(3e14)/wmodel/414/414

            csvalue = testSEDcsq(destdir,0)

            # write the chisq test chisqtestresults
            file_object = open(outputdir+outtxtfile, 'a')
            file_object.write("\n")
            file_object.write(dircode+', '+(str(csvalue[0])[:6])+', '+(str(csvalue[1]))+', '+(str(csvalue[2])[:6])+', '+(str(csvalue[3]))+', '+(str(csvalue[4])[:6])+', '+(str(csvalue[5])))
            file_object.close()

            #fig1
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
            f1.text(0.2, 1e-12, r'$\chi^2_{all}$:'+(str(csvalue[0])[:4]) , fontsize=9)
            f1.text(0.2, 2.5e-13, r'$\chi^2_{10}$:'+(str(csvalue[2])[:4]) , fontsize=9)
            f1.text(0.2, 8e-14, r'$\chi^2_{other}$:'+(str(csvalue[4])[:4]) , fontsize=9)

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
                if i == 0:
                    dustlegend.append(dustelements[i][2:-1])
                elif i == len(dustelements)-1:
                    dustlegend.append(dustelements[i][1:-2])
                else:
                    dustlegend.append(dustelements[i][1:-1])
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


            fig.suptitle(dircode)
            fig.subplots_adjust(hspace=0.35, wspace=0.35)
            pdf.savefig(papertype='a4', orientation='landscape')
            plt.close()



def setupOutFiles(outfilename):
    with open(outfilename, 'w') as wfile:
        wfile.write('%s\n' % '#AA Ori RADMC3D protoplanetary disk model fitting chisquare values')
        wfile.write('%s\n' % 'dirindex, chisq_all, Nall, chisq_y10, N_y10, chisq_n10, N_n10')
        wfile.close()
