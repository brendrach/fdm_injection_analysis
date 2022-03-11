## Import files
from __future__ import print_function
import sys
import libstempo 
import libstempo.plot as LP, libstempo.toasim as LT

import numpy as np
import glob, os, json
import pickle
import random
import csv
import math

import enterprise
from enterprise.pulsar import Pulsar
import enterprise.signals.parameter as parameter
from enterprise.signals import utils
from enterprise.signals import signal_base
from enterprise.signals import selections
from enterprise.signals.selections import Selection
from enterprise.signals import white_signals
from enterprise.signals import gp_signals
from enterprise.signals import deterministic_signals

from PTMCMCSampler.PTMCMCSampler import PTSampler as ptmcmc

import argparse

## Initialize the parser to accept user inputs
parser = argparse.ArgumentParser(description = "Initiate data sets with stochastic GWB injections with given numpy seed (realization) and amplitudes")

## Required arguments:
parser.add_argument("-timpath", required = True, help = "Path to base directory holding timfiles")
parser.add_argument("-parpath", required = True, help = "Path to base directory containing all parfiles")
parser.add_argument("-noisefile", required = True, help = "Location of noisefile")
parser.add_argument("-timsavepath", required = True, help = "Path to directory where we will save the new tim and parfiles")



## Optional arguments:
parser.add_argument("--datacut_length", dest='datacut_length', default=7.0, help="The minimum observing baseline a pulsar must meet in order to be included in the simulation.", type=float)
parser.add_argument("--ephemeris", dest = 'ephem', help = "Choose solar system ephemeris for loading in pulsars; Default: DE438", choices = ['DE430', 'DE435', 'DE436', 'DE438', 'BayesEphem'], default = 'DE438')
parser.add_argument("--gamma", dest = 'gamma',  help = "Specify index of stochastic GWB powerlaw function; Default: 13./3.", type = float, default = 13./3.)
parser.add_argument("--inject_fdm", dest = 'inject_fdm', action = 'store_true', default = False, help = "Whether to inject an FDM signal into our simulated data")
parser.add_argument("--inject_GWB", dest = 'inject_GWB', action = 'store_true', default = False, help = "Whether to inject a GWB signal into our simulated data")
parser.add_argument("--fdm_A", dest = 'fdm_A', help = "The fuzzy dark matter amplitude to inject; Default: None", default=None, type=float)
parser.add_argument("--fdm_f", dest = 'fdm_f', help = "The fuzzy dark matter frequency to inject; Default: None", default=None, type=float)
parser.add_argument("--GWB_A", dest = 'GWB_A', help = "The GWB amplitude to inject; Default: None", default=None, type=float)
parser.add_argument("--average_epoch", dest = 'average_epoch', help = "The epoch to average over", default=None, type=float)

## Load the arguments:
args = parser.parse_args()

print(args.average_epoch)

if args.average_epoch:
    print('YES!')

## Begin writing all arguments to an output file for record keeping
output_params = open(args.timsavepath + "simulation_params.txt","a")
output_params.write(str(args) + '\n')

## Check if the user will be injecting fuzzy dark matter
if args.inject_fdm:

    ## If they do inject fuzzy dark matter, they must specify all of the FDM parameters or the code will exit.
    if args.fdm_A is None or args.fdm_f is None:
        raise ValueError("You have decided to inject a fuzzy dark matter signal but did not specify any or all of the injected values." + \
                     " The values you've provided are fdm_A, fdm_f = " + str(args.fdm_A) + ', ' + str(args.fdm_f) + '.')

## Check if the user will be injecting a gravitational wave background (GWB)
if args.inject_GWB:

    ## If they do, make sure they specified the GWB amplitude
    if args.GWB_A is None:
        raise ValueError("You have decided to inject a gravitational background signal but did not specify the injected amplitude.")


## Get parameter noise dictionary
noise_ng12 = args.noisefile

params = {}
with open(noise_ng12, 'r') as noise:
    params.update(json.load(noise))

def get_noise(psr_name):
    
    efac_keys = []
    efac_vals = []
    for key, value in params.items():
        if psr_name in key and 'efac' in key:
            efac_keys.append(key)
            efac_vals.append(value)
            
           
    equad_keys = []
    equad_vals = []
    for key, value in params.items():
        if psr_name in key and 'equad' in key:
            equad_keys.append(key)
            equad_vals.append(value)
    

    red_noise_gamma_val = [val for key, val in params.items() if psr_name in key and 'red_noise_gamma' in key][0]
    red_noise_amp_val = [val for key, val in params.items() if psr_name in key and 'red_noise_log10_A' in key][0]


    return efac_keys, efac_vals, equad_keys, equad_vals, red_noise_gamma_val, red_noise_amp_val

def epoch_average(psr, dt=1):
    toas, residuals = psr.toas(), psr.residuals()
    
    rtime = dt * np.round(toas * 86400 / dt)
    epochs = np.sort(np.unique(rtime))
    
    ret = []
    for epoch in epochs:
        mask = (rtime == epoch)
        resv, errv = residuals[mask], psr.toaerrs[mask]
        
        avgvar = 1 / np.sum(1/errv**2)
        avgres = np.sum(resv/errv**2) * avgvar

        ret.append((epoch / 86400, avgres, math.sqrt(avgvar)))
    ret = np.array(ret)
    
    return ret[:,0], ret[:,1], ret[:,2]


## Get all of the par/tim files sorted.
parfiles = sorted(glob.glob(args.parpath + '*.par'))
timfiles = sorted(glob.glob(args.timpath + '*.tim'))

psrcut_parfiles = []
psrcut_timfiles = []
for i, parfile in enumerate(parfiles):
    if 'J1713' in parfile and 't2' not in parfile:
        parfiles.remove(parfile)
    else:
        with open(parfile, newline='') as f:
            reader = csv.reader(f)
            for row in reader:
                if 'START' in row[0].split():
                    start = float(row[0].split()[1])
                    
                if 'FINISH' in row[0].split():
                    finish = float(row[0].split()[1])
            
            if (finish-start)/365.25 > args.datacut_length:
                psrcut_parfiles.append(parfiles[i])
                psrcut_timfiles.append(timfiles[i])

## Write the pulsar names to the output file for record keeping.
output_params.write("This dataset includes: \n")
# import each pulsar into tempo2
tempopsr = []

## Loop through the pulsars and make each a tempo2 object.
for i in range(len(psrcut_parfiles)):
    psr = libstempo.tempopulsar(parfile = psrcut_parfiles[i],
                                timfile = psrcut_timfiles[i], maxobs=50000)
             
    print('Creating Tempo object for ' + psr.name)
    output_params.write(psr.name + '\n')

    if args.average_epoch: 
        print("This psr has " + str(len(psr.toas())) + " observations before averaging.")
        toas, res, err = epoch_average(psr, dt=args.average_epoch)
        fp = LT.fakepulsar(psrcut_parfiles[i], toas, err)
        fp.stoas[:] += res / 86400
        print("This psr has " + str(len(fp.toas())) + " observations after averaging.")
        fp.formbats()

        tempopsr.append(fp)

    else:
        tempopsr.append(psr)
    

# remove all noise from the pulsars.
for psr in tempopsr:
    print('Adding rednoise to ' + psr.name)

    efac_keys, efac_vals, equad_keys, equad_vals, red_noise_gamma_val, red_noise_amp_val = get_noise(psr.name)

    if args.average_epoch:
        efac_val = np.mean(efac_vals)
        equad_val = np.mean(equad_vals)

        LT.add_efac(psr, efac=efac_val)
        LT.add_equad(psr, equad=10**equad_val)
        LT.add_rednoise(psr, 10**red_noise_amp_val, red_noise_gamma_val)

    else:
        for i in range(len(equad_vals)):
            equad_vals[i] = 10**equad_vals[i]

        print(equad_vals)

        LT.make_ideal(psr)
        LT.add_efac(psr, efac=efac_vals, flags = efac_keys, flagid = 'f')
        LT.add_equad(psr, equad=equad_vals, flags = equad_keys, flagid = 'f')
        LT.add_rednoise(psr, 10**red_noise_amp_val, red_noise_gamma_val)
    

## Inject fdm.
if args.inject_fdm:
    print("Injecting Fuzzy Dark Matter!")
    fdm_phase_e = random.random() * 2*np.pi
    output_params.write('fdm_phase_e = ' + str(fdm_phase_e) + '\n')
    for psr in tempopsr:
        fdm_phase_p = random.random() * 2*np.pi
        output_params.write(psr.name + ' fdm_phase_p = ' + str(fdm_phase_p) + '\n')
        LT.add_fdm(psr, args.fdm_A, args.fdm_f, fdm_phase_e, fdm_phase_p)

# Inject a GWB.
if args.inject_GWB:
    print("Injecting Gravitational Wave Background")
    LT.createGWB(tempopsr, Amp = args.GWB_A, gam = args.gamma)


for psr in tempopsr:
    print(psr)
    if args.average_epoch:

        print("Cleaning up design matrix")

        for i in range(1,126):
            psr[f'DMX_{i:04d}'].fit = False

        for par in ['JUMP1','FD1','FD2','FD3']:
            psr[par].fit = False

        dmsum = np.sum(psr.designmatrix(updatebats=False,incoffset=True)**2, axis=0)

        pars_to_ignore = []
        for i in range(len(dmsum)-1):
            if dmsum[i+1] == 0:
                pars_to_ignore.append(psr.pars()[i])

        for par in pars_to_ignore:
            psr[par].fit = False

        psr.fit(include_noise=False)
    else:
        psr.fit()

## Save the tim files. 
print("Saving the tim files!")
for psr in tempopsr:

    psr.savetim(args.timsavepath + psr.name + '-sim.tim')
    psr.savepar(args.timsavepath + psr.name + '-sim.par')
    libstempo.purgetim(args.timsavepath + psr.name + '-sim.tim')
    
    
