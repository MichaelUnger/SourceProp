#!/usr/bin/env python
# -*- coding: utf-8 -*-

import emcee
import triangle
import numpy as np
import ROOT
from rootext import ToNumpy

ROOT.gSystem.Load("libProp.so")

nwalkers = 50
nsteps = 500
options = "fitFiles/PaperDefault.txt"
output = "mcmc.root"
interface = ROOT.prop.MCMCInterface(nwalkers, options, output)

def logProb(par):
  n = len(par)
  arg = ROOT.std.vector("double")(n)
  for i in xrange(n):
      arg[i] = par[i]
  return interface.GetLogProb(arg)

startVals = ToNumpy(interface.GetFreeParStartValues())
ndim = startVals.size
pos = [startVals + 1e-3*startVals*np.random.randn(ndim) for i in range(nwalkers)]
sampler = emcee.EnsembleSampler(nwalkers, ndim, logProb)

print("Running MCMC...")
sampler.run_mcmc(pos, nsteps, rstate0 = np.random.get_state(), storechain = False)
print("Done.")

interface.CloseFile()

print("Mean acceptance fraction: {0:.3f}"
      .format(np.mean(sampler.acceptance_fraction)))
