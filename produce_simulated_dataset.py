## Import files
from __future__ import print_function
import sys
import libstempo 
import libstempo.plot as LP, libstempo.toasim as LT

import numpy as np
import glob, os, json
import pickle
import random

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
parser.add_argument("-timsavepath", required = True, help = "Path to directory where we will save the new tim and parfiles")



## Optional arguments:
parser.add_argument("--datacut_length", dest='datacut_length', default=7.0, help="The minimum observing baseline a pulsar must meet in order to be included in the simulation.", type=float)
parser.add_argument("--read_pickle", dest = 'read_pickle', default = False, action = 'store_true', help = "Flag to read dataset as pickle to save i/o time; Default: False")
parser.add_argument("--ephemeris", dest = 'ephem', help = "Choose solar system ephemeris for loading in pulsars; Default: DE438", choices = ['DE430', 'DE435', 'DE436', 'DE438', 'BayesEphem'], default = 'DE438')
parser.add_argument("--gamma", dest = 'gamma',  help = "Specify index of stochastic GWB powerlaw function; Default: 13./3.", type = float, default = 13./3.)
parser.add_argument("--inject_fdm", dest = 'inject_fdm', action = 'store_true', default = False, help = "Whether to inject an FDM signal into our simulated data")
parser.add_argument("--inject_GWB", dest = 'inject_GWB', action = 'store_true', default = False, help = "Whether to inject a GWB signal into our simulated data")
parser.add_argument("--fdm_A", dest = 'fdm_A', help = "The fuzzy dark matter amplitude to inject; Default: None", default=None, type=float)
parser.add_argument("--fdm_f", dest = 'fdm_f', help = "The fuzzy dark matter frequency to inject; Default: None", default=None, type=float)
parser.add_argument("--GWB_A", dest = 'GWB_A', help = "The GWB amplitude to inject; Default: None", default=None, type=float)

## Load the arguments:
args = parser.parse_args()

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


## Check if the user has a pickle file containing the par and tim files.
if args.read_pickle is True:
    psrcut = pickle.load(open(args.parpath + "psrs.p", "rb"))

    ## Get the names of the pulsars that are in pickled dataset.
    psrcut_parfiles = []
    for psr in psrcut:
        for parfile in parfiles:
            if psr.name in parfile:
                psrcut_parfiles.append(parfile)
                
    psrcut_timfiles = []
    for psr in psrcut:
        for timfile in timfiles:
            if psr.name in timfile:
                psrcut_timfiles.append(timfile)

## If not, get the par and tim from their usual txt format. 
else:
    ## Get all of the par/tim files sorted.
    parfiles = sorted(glob.glob(args.parpath + '*.par'))
    timfiles = sorted(glob.glob(args.timpath + '*.tim'))

    ## Get the tempo2 version of J1713
    for parfile in parfiles:
        if 'J1713' in parfile and 't2' not in parfile:
            parfiles.remove(parfile)

    ## Zip the par/tim files for easy access.
    psrs = []
    for p, t in zip(parfiles, timfiles):
        psrname = parfiles[0].split('/')[-1].split('_')[0]
        psr = Pulsar(p, t, ephem='DE438')
        psrs.append(psr)

    ## Eliminate all pulsars that do not have a baseline greater than our specific cutoff.
    psrcut = []
    for psr in psrs:
        tmax = max(psr.toas)
        tmin = min(psr.toas)
        Tspan = tmax - tmin
        if Tspan / 525600 / 60 > args.datacut_length:
            psrcut.append(psr)
            
    ## Get the names of all of our remaining pulsars.
    psrcut_parfiles = []
    for psr in psrcut:
        for parfile in parfiles:
            if psr.name in parfile:
                psrcut_parfiles.append(parfile)
                
    psrcut_timfiles = []
    for psr in psrcut:
        for timfile in timfiles:
            if psr.name in timfile:
                psrcut_timfiles.append(timfile)

## Write the pulsar names to the output file for record keeping.
output_params.write("This dataset includes: \n")
# import each pulsar into tempo2
tempopsr = []
i=0
## Loop through the pulsars and make each a tempo2 object.
for psr in psrcut:
    print(psr.name)
    output_params.write(psr.name + '\n')
    psr = libstempo.tempopulsar(parfile = psrcut_parfiles[i],
                                timfile = psrcut_timfiles[i], maxobs=50000)
    tempopsr.append(psr)
    i += 1

# remove all noise from the pulsars.
for psr in tempopsr:
    LT.make_ideal(psr)

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
    psr.fit()

## Save the tim files. 
print("Saving the tim files!")
for psr in tempopsr:
    psr.savetim(args.timsavepath + psr.name + '-sim.tim')
    psr.savepar(args.timsavepath + psr.name + '-sim.par')
    libstempo.purgetim(args.timsavepath + psr.name + '-sim.tim')
    