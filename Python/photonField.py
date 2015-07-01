import numpy as np
import scipy as sp
from scipy.special import zeta
from os import path
import sys

datadir = 'SED/'

cm = 1e-2
cm3 = cm**3
eV = 1.60217657e-19  # [J]
erg = 1e-7  # [J]
c0 = 299792458  # [m/s]
h = 6.62606957e-34  # [m^2 kg / s]
kB = 1.3806488e-23  # [m^2 kg / s^2 / K]
T_CMB = 2.72548  # CMB temperature [K]
pc = 3.08567758e16 # [m]

# interfaces
class CMB:
    """
    Cosmic microwave background radiation
    """
    name = 'CMB'
    info = 'CMB'
    redshift = None

    def getDensity(self, eps, z=0):
        """
        Comoving spectral number density dn/deps [1/m^3/J] at given photon energy eps [J] and redshift z.
        Multiply with (1+z)^3 for the physical number density.
        """
        return 8*np.pi / c0**3 / h**3 * eps**2 / (np.exp(eps/(kB*T_CMB)) - 1)

    def getEmin(self, z=0):
        """Minimum effective photon energy in [J]"""
        return 1e-10 * eV

    def getEmax(self, z=0):
        """Maximum effective photon energy in [J]"""
        return 0.1 * eV

def temperatureFromPeak(peakEnergy, sigma):
    b = sigma + 2
    x = (sp.special.lambertw(-np.exp(-b) * b) + b).real
    return peakEnergy / x / kB

def peakFromTemperature(T, sigma):
    b = sigma + 2
    x = (sp.special.lambertw(-np.exp(-b) * b) + b).real
    return x * (kB*T)


class ModifiedBlackBody:
    """
    modified black body
    """
    name = 'MBB'
    info = 'MBB'
    redshift = None

    def __init__(self, T=90, sigma=1.) :
        self.T = T
        self.sigma = sigma
        self.name = 'MBB_' + str(int(round(T,0))) + "_" + str(sigma)
        if sigma > 0:
            self.info = 'MBB, $T=' + str(int(round(T,0))) + "$ K, $\sigma=" + \
                        str(sigma) + "$"
        else:
            self.info = 'BB, $T=' + str(int(round(T,0))) + "$ K"
        integralMBB = T**3 * zeta(sigma + 3., 1) * sp.special.gamma(sigma + 3.)
        peakEnergy = peakFromTemperature(T, sigma)
        TBB = temperatureFromPeak(peakEnergy, 0)
        integralBB = TBB**3 * zeta(3., 1) * sp.special.gamma(3.)
        self.integralRatio = integralBB / integralMBB

    @classmethod
    def fromPeak(cls, peakEnergy, sigma):
        """create from peakEnergy and sigma """
        T = temperatureFromPeak(peakEnergy, sigma)
        return cls(T, sigma)

    def getDensity(self, eps, z=0):
        """
        Spectral number density dn/deps [1/m^3/J] at given photon energy eps [J]
        """
        kT = self.T * kB
        bb = 8*np.pi / (c0 * h)**3 * eps**2 / (np.exp(eps/(kT)) - 1)
        return self.integralRatio * bb * (eps/kT)**self.sigma

    def getEmin(self, z=0):
        """Minimum effective photon energy in [J]"""
        return 1e-10 * eV

    def getEmax(self, z=0):
        """Maximum effective photon energy in [J]"""
        return 100000 * eV


class BrokenPowerLaw:
    """
    Broken Power Law Photon Spectrum
    """
    name = 'SP_0.05_2.0_52'
    info = 'SP_0.05_2.0_52'
    redshift = None

    def __init__(self, eps0 = 0.05, beta = -2, nom = 5, denom = 2) :
        self.eps0 = eps0
        self.beta = float(beta)
        self.alpha = float(nom) / denom
        self.name = 'BPL_' + str(eps0/eV) + "_" + \
                    str(abs(round(beta, 1))) + "_" + str(nom) + "/" + str(denom);
        self.info = 'BPL, $\epsilon_0 = ' + str(eps0/eV) + \
                    '$ eV, $\\alpha=' + str(nom) + "/" + str(denom) + \
                    "$, $\\beta=" + str(round(beta, 1)) + "$"
        TBB = temperatureFromPeak(eps0, 0)
        integralBB = 8*np.pi / (c0 * h)**3 * (TBB*kB)**3 * zeta(3., 1) * \
                     sp.special.gamma(3.)
        integralBPL = eps0 * ( 1. / (self.alpha + 1) - 1. / (self.beta + 1))
        print integralBPL
        self.integralRatio = integralBB / integralBPL


    def getDensity(self, eps, z=0):
        """
        Spectral number density dn/deps [1/m^3/J] at given photon energy eps [J]
        """
        x = eps / self.eps0
        C = self.integralRatio
        alpha = self.alpha
        beta = self.beta
        return ((eps > self.eps0)*x**beta+(eps <= self.eps0)*x**alpha)*C

    def getEmin(self, z=0):
        """Minimum effective photon energy in [J]"""
        return 1e-10 * eV

    def getEmax(self, z=0):
        """Maximum effective photon energy in [J]"""
        return 100000 * eV



class BlackBody:
    """
    black body
    """
    name = 'BB'
    info = 'BB'
    redshift = None

    def __init__(self, T=90) :
        self.T = T
        self.name = 'BB_' + str(T)

    def getDensity(self, eps, z=0):
        """
        """
        kT = self.T * kB

        return 8*np.pi / (c0 * h)**3 * eps**2 / (np.exp(eps/(kT)) - 1)

    def getEmin(self, z=0):
        """Minimum effective photon energy in [J]"""
        return 1e-10 * eV

    def getEmax(self, z=0):
        """Maximum effective photon energy in [J]"""
        return 100000 * eV

class Elbaz2011:
    """
    Elbaz 2011 starburst arXiv:1105.2537
    """
    name = 'Elbaz11'
    info = 'Elbaz11'
    redshift = None

    def __init__(self, prolong = True) :
        lmbda, L = np.genfromtxt(datadir + 'starburst_sed_Elbaz2011.tbl', unpack=True)
        self.eps = (h*c0 / (lmbda * 1e-6))[::-1]
        Lreverse = L[::-1]
        Lsun = 3.846e26 # J/sec
        N = 1 # number of starforming regions
        R = 100*pc # size of starforming regions
        self.dndeps = Lreverse * Lsun /(N*np.pi*R**2*c0*self.eps**2)
        if prolong :
            epsLast = self.eps[-1]
            dndepsLast = self.dndeps[-1]
            epsFirst = self.eps[0]
            dndepsFirst = self.dndeps[0]

            eps = np.logspace(np.log10(epsLast/eV), 2, 200)[1:] * eV
            self.eps = np.concatenate((self.eps, eps))
            self.dndeps = np.concatenate((self.dndeps, dndepsLast*(eps/epsLast)**(-2.5)))

            eps = np.logspace(-5, np.log10(epsFirst/eV), 200)[:-1] * eV
            self.eps = np.concatenate((eps, self.eps))
            self.dndeps = np.concatenate((dndepsFirst*(eps/epsFirst)**(2.5), self.dndeps))


    def getDensity(self, eps, z=0):
        """
        """
        return np.interp(eps, self.eps, self.dndeps, 0, 0)

    def getEmin(self, z=0):
        """Minimum effective photon energy in [J]"""
        return self.eps[0]

    def getEmax(self, z=0):
        """Maximum effective photon energy in [J]"""
        return self.eps[-1]


if __name__ == '__main__':
    from pylab import *
    eps = logspace(-4, -0.4, 200) * eV
    x  = eps / eV
    eps0 = 0.05 * eV
    bb = ModifiedBlackBody.fromPeak(eps0, 0)
    mbb1 = ModifiedBlackBody.fromPeak(eps0, 1)
    mbb2 = ModifiedBlackBody.fromPeak(eps0, 2)
    bpl = BrokenPowerLaw(eps0)
    y1 = bb.getDensity(eps) / (1/eV) / (1/cm3)
    y2 = mbb1.getDensity(eps) / (1/eV) / (1/cm3)
    y3 = mbb2.getDensity(eps) / (1/eV) / (1/cm3)
    y4 = bpl.getDensity(eps) / (1/eV) / (1/cm3)
    figure()
    plot(x, y1, label=bb.info)
    plot(x, y2, label=mbb1.info)
    plot(x, y3, label=mbb2.info)
    plot(x, y4, label=bpl.info)

    legend(loc='upper right')
#    semilogx()
#    loglog()
    ylim(ymin=0)
    ylabel('d$n/$d$\epsilon$ [eV$^{-1}\,$cm$^{-3}$]')
    xlabel('$\epsilon$ [eV]')
    savefig('photon.png')
    savefig('photon.pdf')
