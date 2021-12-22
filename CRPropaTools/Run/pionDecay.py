from crpropa import *

nEvents = 10000000
doNeutrinos = True
doPhotons = False
spectralIndex = -1.
minEnergyObs = 1*TeV
minEnergy  = 1e17*eV
maxEnergy = 1e22*eV
minDistance = 0.5*Mpc
maxDistance = 12*Gpc


#Random_seedThreads(123)

# simulation setup
H0 = 67.8 * 1000 * meter / second / Mpc
# account for funny conversion in setParameters(double h, double oM)
setCosmologyParameters(H0 * Mpc / 1e5, 0.308)
sim = ModuleList()
sim.add(SimplePropagation())
sim.add(Redshift())
sim.add(PionDecay(doPhotons, doNeutrinos))
sim.add(OutputROOTEvent('crpropa.root'))
sim.add(MinimumEnergy(minEnergyObs))

# observer
obs = Observer()
obs.add(ObserverPoint())
#obs.add(ObserverOutput1D('events_sim1D.txt'))
sim.add(obs)

# source
source = Source()
source.add(SourceUniform1D(minDistance, maxDistance))
source.add(SourceRedshift1D())

#http://en.wikipedia.org/wiki/Table_of_nuclides_%28complete%29
composition = SimpleComposition(minEnergy, maxEnergy, spectralIndex)
composition.add(211)  # pi+
composition.add(-211)  # pi-
source.add(composition)

# run simulation
sim.setShowProgress(True)
sim.run(source, nEvents, True)
