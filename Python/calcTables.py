from numpy import *
import photonField
import interactionRate as iR


eV = 1.60217657e-19


crpDataDir = '/ssd/munger/Mag/CRPropa3New/CRPropa3-data/'


temperatures = [10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 175, 200, 225, 250, 275, 300, 350, 400, 450, 500, 750, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000]
temperatures = [95]

for temperature in temperatures:


    gamma = logspace(6, 14, 201)  # tabulated UHECR Lorentz-factors


    # cross sections for A < 12 (various sources)
    ddir1 = crpDataDir + 'tables/PD_external/'
    isotopes1 = genfromtxt(ddir1 + 'isotopes.txt')
    x = genfromtxt(ddir1+'eps.txt') * eV * 1e6  # [J]
    n = len(x)
    d1sum = genfromtxt(ddir1+'xs_sum.txt',
                       dtype=[('Z',int), ('N',int), ('xs','%if8'%n)])
    d1exc = genfromtxt(ddir1+'xs_excl.txt',
                       dtype=[('Z',int), ('N',int), ('ch',int), ('xs','%if8'%n)])
    eps1 = iR.romb_pad_logspaced(x, 513)  # padding
    xs1sum = array([iR.romb_pad_zero(x, 513) for x in d1sum['xs']])*1e-31
    xs1exc = array([iR.romb_pad_zero(x, 513) for x in d1exc['xs']])*1e-31


    # cross sections for A >= 12 (TALYS)
    ddir2 = crpDataDir + 'tables/PD_Talys1.6_Khan/'
    isotopes2 = genfromtxt(ddir2 + 'isotopes.txt')
    x = genfromtxt(ddir2+'eps.txt') * eV * 1e6  # [J]
    n = len(x)
    d2sum = genfromtxt(ddir2+'xs_sum.txt',
                       dtype=[('Z',int), ('N',int), ('xs','%if8'%n)])
    d2exc = genfromtxt(ddir2+'xs_thin.txt',
                       dtype=[('Z',int), ('N',int), ('ch',int), ('xs','%if8'%n)])
    eps2 = iR.romb_pad_logspaced(x, 513)  # padding
    xs2sum = array([iR.romb_pad_zero(x, 513) for x in d2sum['xs']])*1e-31
    xs2exc = array([iR.romb_pad_zero(x, 513) for x in d2exc['xs']])*1e-31

    fields = [
        photonField.BrokenPowerLaw(temperature/1000., -2, 3, 2)
#        photonField.ModifiedBlackBody(T=363, sigma=0),
#        photonField.ModifiedBlackBody(T=206, sigma=1),
#        photonField.ModifiedBlackBody(T=148, sigma=2),
        ]

    print '-------- PD'
    for field in fields:
        print field.name

        # Calculate total interaction rate
        R1 = array([iR.invMFP_fast(eps1, x, gamma, field) for x in xs1sum])
        R2 = array([iR.invMFP_fast(eps2, x, gamma, field) for x in xs2sum])

        # save
        fname = 'data/pd_%s.txt' % field.name
        output = r_[
            c_[d1sum['Z'], d1sum['N'], R1],
            c_[d2sum['Z'], d2sum['N'], R2]]
        fmt = '%i\t%i' + '\t%g'*201
        hdr = 'Photo-disintegration with the %s\nZ, N, 1/lambda [1/Mpc] for log10(gamma) = 6-14 in 201 steps' % field.info
        savetxt(fname, output, fmt=fmt, header=hdr)


        # Calculate branching ratios
        # for A < 12
        B1 = array([iR.invMFP_fast(eps1, x, gamma, field) for x in xs1exc])
        for (Z, N, A) in isotopes1:
            s = (d1exc['Z'] == Z) * (d1exc['N'] == N)
            B1[s] /= sum(B1[s], axis=0)
            B1[isnan(B1)] = 0  # set to 0 when total cross section is 0

        # for A > 12
        B2 = array([iR.invMFP_fast(eps2, x, gamma, field) for x in xs2exc])
        for (Z, N, A) in isotopes2:
            s = (d2exc['Z'] == Z) * (d2exc['N'] == N)
            B2[s] /= sum(B2[s], axis=0)
            B2[isnan(B2)] = 0  # set to 0 when total cross section is 0

        # save
        fname = 'data/pd_branching_%s.txt'%field.name
        output = r_[
            c_[d1exc['Z'], d1exc['N'], d1exc['ch'], B1],
            c_[d2exc['Z'], d2exc['N'], d2exc['ch'], B2]]
        fmt = '%i\t%i\t%06d' + '\t%g'*201
        hdr = 'Photo-disintegration with the %s\nZ, N, channel, branching ratio for log10(gamma) = 6-14 in 201 steps' % field.info
        savetxt(fname, output, fmt=fmt, header=hdr)


    # proton / neutron cross sections [1/m^2] for tabulated energies [J]
    # truncate to largest length 2^i + 1 for Romberg integration
    d = genfromtxt(crpDataDir + 'tables/PPP/xs_proton.txt', unpack=True)
    eps1 = d[0] * 1e9 * eV  # [J]
    xs1  = d[1] * 1e-34  # [m^2]

    d = genfromtxt(crpDataDir + 'tables/PPP/xs_neutron.txt', unpack=True)
    eps2 = d[0] * 1e9 * eV  # [J]
    xs2  = d[1] * 1e-34  # [m^2]

    lgamma = linspace(6, 16, 251)
    gamma  = 10**lgamma


    print '-------- PP'

    for field in fields:
        print field.name
        r1 = iR.invMFP_fast(eps1[0:2049], xs1[0:2049], gamma, field)
        r2 = iR.invMFP_fast(eps2[0:2049], xs2[0:2049], gamma, field)

        fname = 'data/ppp_%s.txt' % field.name
        data  = c_[lgamma, r1, r2]
        fmt   = '%.2f\t%.6e\t%.6e'
        header = 'Photo-pion interaction rate with the %s\nlog10(gamma)\t1/lambda_proton [1/Mpc]\t1/lambda_neutron [1/Mpc]'%field.info
        savetxt(fname, data, fmt=fmt, header=header)
