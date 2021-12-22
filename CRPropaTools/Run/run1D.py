from crpropa import *

nEvents = 500000
doNeutrinos = True
doPhotons = False
spectralIndex = -1.
minEnergyObs = 1*TeV
minEnergy  = 1e17*eV
maxEnergy = 1e22*eV
minDistance = 0.5*Mpc
maxDistance = 12*Gpc

infrared = IRB_Stecker05

#Random_seedThreads(123)

# simulation setup
H0 = 67.8 * 1000 * meter / second / Mpc
# account for funny conversion in setParameters(double h, double oM)
setCosmologyParameters(H0 * Mpc / 1e5, 0.308)
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
composition.add(2112)  # n
composition.add(12)  # nu_e
# PSB chain
composition.add(1, 1)  # H
composition.add(2, 1)  # D
composition.add(3, 2)  # He-3
composition.add(4, 2)  # He-4
composition.add(5, 3)  # Li-5
composition.add(6, 3)  # Li-6
composition.add(7, 3)  # Li-7
composition.add(8, 4)  # Be-8
composition.add(9, 4)  # Be-9
composition.add(10, 5)  # B-10
composition.add(11, 5)  # B-11
composition.add(12, 6)  # C-12
composition.add(13, 6)  # C-13
composition.add(14, 7)  # N-14
composition.add(15, 7)  # N-15
composition.add(16, 8)  # O-16
composition.add(17, 8)  # O-17
composition.add(18, 8)  # O-18
composition.add(19, 9)  # F-19
composition.add(20, 10)  # Ne-20
composition.add(21, 10)  # Ne-21
composition.add(22, 10)  # Ne-22
composition.add(23, 11)  # Na-23
composition.add(24, 12)  # Mg-24
composition.add(25, 12)  # Mg-25
composition.add(26, 12)  # Mg-26
composition.add(27, 13)  # Al-27
composition.add(28, 14)  # Si-28
composition.add(29, 14)  # Si-29
composition.add(30, 14)  # Si-30
composition.add(31, 15)  # P-31
composition.add(32, 16)  # S-32
composition.add(33, 16)  # S-33
composition.add(34, 16)  # S-34
composition.add(35, 17)  # Cl-35
composition.add(36, 16)  # S-36
composition.add(37, 17)  # Cl-37
composition.add(38, 18)  # Ar-38
composition.add(39, 19)  # K-39
composition.add(40, 20)  # Ca-40
composition.add(41, 20)  # Ca-41
composition.add(42, 20)  # Ca-42
composition.add(43, 20)  # Ca-43
composition.add(44, 20)  # Ca-44
composition.add(45, 21)  # Sc-45
composition.add(46, 22)  # Ti-46
composition.add(47, 22)  # Ti-47
composition.add(48, 22)  # Ti-48
composition.add(49, 22)  # Ti-49
composition.add(50, 22)  # Ti-50
composition.add(51, 23)  # V-51
composition.add(52, 24)  # Cr-52
composition.add(53, 24)  # Cr-53
composition.add(54, 24)  # Cr-54
composition.add(55, 25)  # Mn-55
composition.add(56, 26)  # Fe-56
source.add(composition)

# run simulation
sim.setShowProgress(True)
sim.run(source, nEvents, True)
