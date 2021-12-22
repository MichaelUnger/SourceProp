import numpy as np
from crpropa import *

nEvents = 500000
doNeutrinos = True
doPhotons = False
spectralIndex = -1.
minEnergyObs = 1*TeV
minEnergy  = 1e19*eV
maxEnergy = 5e20*eV

distances = np.loadtxt('agnDist.txt')

infrared = IRB_Gilmore12

#Random_seedThreads(123)

# simulation setup
H0 = 67.8 * 1000 * meter / second / Mpc
# account for funny conversion in setParameters(double h, double oM)
setCosmologyParameters(H0 * Mpc / 1e5, 0.308)

for d in distances:
    filename = "sbSim/agn2_" + str(int(d*10)) + ".root"
    print "==========================================================="
    print filename
    sim = ModuleList()
    sim.add(SimplePropagation())
    sim.add(Redshift())
    sim.add(PhotoPionProduction(CMB, doPhotons, doNeutrinos))
    sim.add(PhotoPionProduction(infrared, doPhotons, doNeutrinos))
    sim.add(PhotoDisintegration(CMB))
    sim.add(PhotoDisintegration(infrared))
    sim.add(NuclearDecay(doPhotons, doNeutrinos))
    sim.add(ElectronPairProduction(CMB))
    sim.add(ElectronPairProduction(infrared))
    sim.add(OutputROOTEvent(filename))
    sim.add(MinimumEnergy(minEnergyObs))

    # observer
    obs = Observer()
    obs.add(ObserverPoint())
    sim.add(obs)

    # source
    source = Source()
    source.add(SourceUniform1D(d*Mpc, (d*1.0001*Mpc)))
    source.add(SourceRedshift1D())

    composition = SimpleComposition(minEnergy, maxEnergy, spectralIndex)
    composition.add(1, 1)  # H
    composition.add(4, 2)  # He-4
    composition.add(14, 7)  # N-14
    composition.add(28, 14)  # Si-28
    composition.add(56, 26)  # Fe-56
    source.add(composition)

    # run simulation
    sim.setShowProgress(True)
    sim.run(source, nEvents, True)
