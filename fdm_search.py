from __future__ import division

import os, glob, json, pickle
import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg as sl
import sys


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
import enterprise.constants as const

import corner
from PTMCMCSampler.PTMCMCSampler import PTSampler as ptmcmc
import enterprise_extensions
from enterprise_extensions import models, model_utils, hypermodel
from enterprise_extensions import sampler as samp

import argparse

## Initialize the parser to accept user inputs
parser = argparse.ArgumentParser(description = "Initiate data sets with stochastic GWB injections with given numpy seed (realization) and amplitudes")

## Required arguments:
parser.add_argument("-timpath", required = True, help = "Path to base directory holding timfiles")
parser.add_argument("-parpath", required = True, help = "Path to base directory containing all parfiles")
parser.add_argument("-outdir", required = True, help = "Path to directory where we will save the output")
parser.add_argument("-noisepath", required = True, help = "Path to directory where the noise parameters are stored")



## Load the arguments:
args = parser.parse_args()

parfiles = sorted(glob.glob(args.parpath + '/*.par'))
timfiles = sorted(glob.glob(args.timpath + '/*.tim'))

psrlist = []
for timfile in timfiles:
    temp_split = timfile.split('/')[-1].split('-')
    if len(temp_split) == 3:
        psrlist.append(temp_split[0] + '-' + temp_split[1])
    else:
        psrlist.append(temp_split[0])

print(psrlist)
   
psr_parfiles = []
psr_timfiles = []

for fname in parfiles:
    for psr in psrlist:
        if psr in fname:
            if 'J1713+0747_NANOGrav_12yv2.gls.par' in fname:
                print(fname)
                pass
            else:
                psr_parfiles.append(fname)
            
for fname in timfiles:
    for psr in psrlist:
        if psr in fname:
            psr_timfiles.append(fname)

print(psr_parfiles)
print(psr_timfiles)
            

psrs = []
for p, t in zip(psr_parfiles, psr_timfiles):
    psr = Pulsar(p, t, ephem='DE438', clk='BIPM(2018)')
    psrs.append(psr)

## Get parameter noise dictionary
noise_ng12 = os.path.join(args.noisepath, 'channelized_12p5yr_v3_full_noisedict.json')

params = {}
with open(noise_ng12, 'r') as fp:
    params.update(json.load(fp))

print(params)


pta =  models.model_fdm(psrs, noisedict=params, white_vary=False, tm_svd=False,
            Tmin_fdm=None, Tmax_fdm=None, gw_psd='powerlaw',
            red_psd='powerlaw', components=5, n_rnfreqs = 30, 
            n_gwbfreqs=5, gamma_common=None, delta_common=None,
            dm_var=False, dm_psd='powerlaw', dm_annual=False,
            upper_limit=False, bayesephem=False, wideband=False,
            pshift=False, pseed=None, model_CRN=True,
            amp_upper=-11, amp_lower=-18, 
            freq_upper=-7, freq_lower=-9,
            use_fixed_freq=False, fixed_freq=-8)


x0 = np.hstack([p.sample() for p in pta.params])
ndim = len(x0)

# initial jump covariance matrix
cov = np.diag(np.ones(ndim) * 0.01**2) # helps to tune MCMC proposal distribution

groups = samp.get_parameter_groups(pta)

# sampler object
sampler = ptmcmc(ndim, pta.get_lnlikelihood, pta.get_lnprior, cov,
                 outDir=args.outdir, groups=groups,
                 resume=False)

# sampler for N steps
N = int(5e6)

jp = samp.JumpProposal(pta)
sampler.addProposalToCycle(jp.draw_from_prior, 5)
sampler.addProposalToCycle(jp.draw_from_fdm_prior, 30)

sampler.addProposalToCycle(jp.draw_from_par_prior(['fdm_log10_A']), 20) #draw from uniform Mc
sampler.addProposalToCycle(jp.draw_from_par_prior(['fdm_log10_f']), 50)
sampler.addProposalToCycle(jp.draw_from_par_prior(['fdm_phase_e']), 20)

np.savetxt(args.outdir + "/pars.txt", list(map(str, pta.param_names)), fmt="%s")
np.savetxt(
    args.outdir + "/priors.txt",
    list(map(lambda x: str(x.__repr__()), pta.params)),
    fmt="%s",
)


# SCAM = Single Component Adaptive Metropolis
# AM = Adaptive Metropolis
# DE = Differential Evolution
## You can keep all these set at default values
#sampler.sample(x0, N, SCAMweight=30, AMweight=15, DEweight=50,)
#sampler = super_model.setup_sampler(resume=False, outdir=outDir, sample_nmodel=False)
sampler.sample(x0, N, SCAMweight=25, AMweight=40, DEweight=0,  writeHotChains = True,)

#file.close()
