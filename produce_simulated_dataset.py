from __future__ import print_function
import sys
import libstempo 
import libstempo.plot as LP, libstempo.toasim as LT

import numpy as np
import glob, os, json
import pickle

import matplotlib.pyplot as plt
import corner

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

#Initialize the parser to accept user inputs
parser = argparse.ArgumentParser(description = "Initiate data sets with stochastic GWB injections with given numpy seed (realization) and amplitudes")

#Required arguments:
parser.add_argument("-outdir", required = True, help = "Base directory to store output chain and parameter files")
parser.add_argument("-timpath", required = True, help = "Path to base directory holding timfiles")
parser.add_argument("-parpath", required = True, help = "Path to directory containing all parfiles")


#Optional arguments:
parser.add_argument("--datacut_length", dest='datacut_length', default=7.0, help="The minimum observing baseline a pulsar must meet in order to be included in the simulation.")
parser.add_argument("--read_pickle", dest = 'read_pickle', action = 'store', default = True, help = "Flag to read dataset as pickle to save i/o time; Default: True")
parser.add_argument("--ephemeris", dest = 'ephem', help = "Choose solar system ephemeris for loading in pulsars; Default: DE438", choices = ['DE430', 'DE435', 'DE436', 'DE438', 'BayesEphem'], default = 'DE4368')
parser.add_argument("--gamma", dest = 'gamma',  help = "Specify index of stochastic GWB powerlaw function; Default: 13./3.", type = float, default = 13./3.)
parser.add_argument("--save_pickle", dest = 'save_pickle', action = 'store', default = False, help = "Flag to save dataset as pickle to save i/o time; Default: False")
parser.add_argument("--inject_fdm", dest = 'inject_fdm', default = False, help = "Whether to inject an FDM signal into our simulated data")
parser.add_argument("--inject_GWB", dest = 'inject_GWB', default = False, help = "Whether to inject a GWB signal into our simulated data")
parser.add_argument("--fdm_A", dest = 'fdm_A', help = "The fuzzy dark matter amplitude to inject; Default: None", default=None)
parser.add_argument("--fdm_f", dest = 'fdm_f', help = "The fuzzy dark matter frequency to inject; Default: None", default=None)
parser.add_argument("--fdm_phase_p", dest = 'fdm_phase_p', help = "The fuzzy dark matter pulsar term phase to inject; Default: None", default=None)
parser.add_argument("--fdm_phase_e", dest = 'fdm_phase_e', help = "The fuzzy dark matter earth term phase to inject; Default: None", default=None)
parser.add_argument("--GWB_A", dest = 'GWB_A', help = "The GWB amplitude to inject; Default: None", default=None)

#Load the arguments:
args = parser.parse_args()

print(type(args.datacut_length))

## Check if the user will be injecting fuzzy dark matter
if args.inject_fdm:

    ## If they do inject fuzzy dark matter, they must specify all of the FDM parameters or the code will exit.
    if args.fdm_A is None or args.fdm_f is None or args.fdm_phase_p is None or args.fdm_phase_e is None:
        raise ValueError("You have decided to inject a fuzzy dark matter signal but did not specify any or all of the injected values." + \
                     " The values you've provided are fdm_A, fdm_f, fdm_phase_p, fdm_phase_e = " + str(fdm_A) + ', ' + str(fdm_f) + ', ' \
                     + str(fdm_phase_p) + ', ' + str(fdm_phase_e))

if args.inject_GWB:

    if args.GWB_A is None:
        raise ValueError("You have decided to inject a gravitational background signal but did not specify the injected amplitude.")

outdir = args.outdir

if args.read_pickle is True:
    print(args.read_pickle)
    psrcut = pickle.load(open("psrs.p", "rb"))

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


else:
    parfiles = sorted(glob.glob(args.parpath + '*.par'))
    timfiles = sorted(glob.glob(args.timpath + '*.tim'))

    for parfile in parfiles:
        if 'J1713' in parfile and 't2' not in parfile:
            parfiles.remove(parfile)

    psrs = []
    for p, t in zip(parfiles, timfiles):
        psrname = parfiles[0].split('/')[-1].split('_')[0]
        psr = Pulsar(p, t, ephem='DE438')
        psrs.append(psr)

    psrcut = []
    for psr in psrs:
        tmax = max(psr.toas)
        tmin = min(psr.toas)
        Tspan = tmax - tmin
        if Tspan / 525600 / 60 > float(args.datacut_length):
            psrcut.append(psr)
        # print(psr.name)
    print(len(psrcut))

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

# import each pulsar into tempo2
tempopsr = []
i=0
for psr in psrcut:
    print(psr.name)
    psr = libstempo.tempopulsar(parfile = psrcut_parfiles[i],
                                timfile = psrcut_timfiles[i], maxobs=50000)
    tempopsr.append(psr)
    i += 1

# remove all noise from the pulsars
# add in EFAC noise
for psr in tempopsr:
    print(psr.name)
    LT.make_ideal(psr)
    #LT.add_efac(psr, efac=1.0)

if args.inject_fdm:
    print("Injecting Fuzzy Dark Matter!")
    for psr in tempopsr:
        LT.add_fdm(psr, args.fdm_A, args.fdm_f, args.fdm_phase_e, args.fdm_phase_p)
        psr.fit()

if args.inject_GWB:
    print("Injecting Gravitational Wave Background")
    for psr in tempopsr:
        LT.createGWB(psr, Amp = args.GWB_A, gam = args.gamma)


        
