# Modified 2 things inside Pyspeckit:
# 1)In measurements.py: 'ptol' from 2 to 8A (interval to find em. lines)
# 2)In optical.py: 'def get_optical_lines()': not to find H3-2 and H4-2

import pyspeckit
import asciitable
import pyfits
import numpy
import numpy as np
from pylab import *
import math
from matplotlib.font_manager import FontProperties
prop = matplotlib.font_manager.FontProperties(size=11)
import matplotlib.pyplot as plt
from os.path import basename, splitext
import csv
from astroquery.irsa_dust import IrsaDust  # IrsaDust needed for the ebv value

# Rest wavelengths of the lines we are fitting - use as initial guesses
Hbeta = 4861.33
OIIIa = 4958.92
OIIIb = 5006.84
OI = 6300.30
NIIa = 6549.86
NIIb = 6585.27
Halpha = 6564.614
SIIa = 6718.29
SIIb = 6732.68
# Offsets between OIII lines
OIIIb_off = OIIIb - OIIIa
OI_off_Ha = OI - Halpha


# Load the SDSS spectra


def loadsdss(file, z):
    hdulist = pyfits.open(file+'.fits')
    VAC = (10**hdulist[1].data.loglam)
    wave = []
    x = 2.735182E-4
    for i in range(0, len(VAC)):
        wave.append(VAC[i]/(1.+x+131.4182/VAC[i]**2+2.76249E8/VAC[i]**4)/(1+z))
    flux = hdulist[1].data.flux*10**-17
    err = hdulist[1].data.ivar*10**-17
    #bunit = hdulist[0].header['bunit']
    c0 = hdulist[0].header['coeff0']
    c1 = hdulist[0].header['coeff1']
    units = 'erg/s/cm^2/Ang'
    xarr = pyspeckit.units.SpectroscopicAxis(wave, units='angstroms')
    spec = pyspeckit.OpticalSpectrum(header=hdulist[0].header, xarr=xarr,
                                     data=flux*1e17, err=err)
    return spec


# Galactic extinction correction


def galextinction(file, spec):  # Take the ebv of the galaxy from IrsaDust
    name = file
    table = IrsaDust.get_query_table(name, section='ebv')
    ebv = table['ext SFD mean'][0]

    spec.deredden(ebv=ebv)  # Deredden in place
    return spec


# Load spectra from the trump folder


def loadtrump(file, z):
    hdulist = pyfits.open(file)
    wave, flux, err = [], [], []
    for i in range(0, len(hdulist[0].data)):
        flux.append(hdulist[0].data[i])
        err.append(hdulist[1].data[i])
        wave.append(hdulist[2].data[i]/(1+z))

    xarr = pyspeckit.units.SpectroscopicAxis(wave, units='angstroms',
                                             refX_units='angstroms')
    spec = pyspeckit.OpticalSpectrum(xarr=xarr, data=np.asarray(flux),
                                     error=np.asarray(err))
    spec.units = 'erg s^{-1} cm^{-2} \\AA^{-1}'
    spec.xarr.units = 'angstroms'
    return spec


# Load spectra from the Kitt folder


def get_wstart(ref, wave_ref, wave_per_pixel):
    return wave_ref - ((ref-1)*wave_per_pixel)


def get_wavelength(start_wave, wave_per_pixel, size):
    return np.array([start_wave + i*wave_per_pixel for i in range(size)])


def check(file):
    hdulist = pyfits.open(file)
    print hdulist[0].header.keys()
    print hdulist[0].header.ascardlist()
    return


def loadkitt(file, z):  # 3 columns are fluxes: 1st, 2nd are improved, 3rd not
    spectra_all_path = basename(file)
    hdulist = pyfits.open(spectra_all_path)
    flux = hdulist[0].data[0][0]
    flux = flux/numpy.mean(flux)  # Flux normalisation
    ref_pixel = hdulist[0].header['CRPIX1']
    coord_ref_pixel = hdulist[0].header['CRVAL1']
    wave_pixel = hdulist[0].header['CD1_1']  # CDELT1
    wstart = get_wstart(ref_pixel, coord_ref_pixel, wave_pixel)
    wave = get_wavelength(wstart, wave_pixel, len(flux))/(1+z)
    f = plt.figure()
    xarr = pyspeckit.units.SpectroscopicAxis(wave, units='angstroms',
                                             refX_units='angstroms')
    spec = pyspeckit.OpticalSpectrum(xarr=xarr, data=flux)
    spec.units = 'erg s^{-1} cm^{-2} \\AA^{-1}'
    spec.xarr.units = 'angstroms'
    return spec


# Function to fit the OI/Halpha/NII/SII lines


def fit_OI(file, extinction, z, sigma):

    # Initialize spectrum object and plot region surrounding Halpha-[NII]
    spec = loadsdss(file, z)
    # spec = loadtrump(file, z)
    # spec = loadkitt(file, z)

    galextinction(extinction, spec)  # Take galactic extinction into account

    spec.plotter(xmin=OI-20, xmax=SIIb + 20)

    # Baseline : Erase the background between the lines
    spec.baseline(xmin=5800, xmax=7000, exclude=[4850, 4870, 4940, 4970, 4995,
                  5020, 6250, 6350, 6450, 6746, 6815, 6884, 7003, 7126, 7506,
                  7674, 8142, 8231], subtract=True, reset_selection=True,
                  highlight_fitregion=True, order=1)

    # Calculate the amplitude of OI, NII, Halpha and SII
    spec.specfit.selectregion(xmin=OI - 10, xmax=OI + 10)
    ampOI = np.max(spec.data[spec.specfit.xmin:spec.specfit.xmax])

    spec.specfit.selectregion(xmin=NIIa - 10, xmax=NIIa + 10)
    ampNIIa = np.max(spec.data[spec.specfit.xmin:spec.specfit.xmax])

    spec.specfit.selectregion(xmin=NIIb - 10, xmax=NIIb + 10)
    ampNIIb = np.max(spec.data[spec.specfit.xmin:spec.specfit.xmax])

    spec.specfit.selectregion(xmin=Halpha - 10, xmax=Halpha + 10)
    ampHa = np.max(spec.data[spec.specfit.xmin:spec.specfit.xmax])

    spec.specfit.selectregion(xmin=SIIa - 10, xmax=SIIa + 10)
    ampSIIa = np.max(spec.data[spec.specfit.xmin:spec.specfit.xmax])

    spec.specfit.selectregion(xmin=SIIb - 10, xmax=SIIb + 10)
    ampSIIb = np.max(spec.data[spec.specfit.xmin:spec.specfit.xmax])

    spec.plotter(xmin=OI-20, xmax=SIIb + 20, ymin=-ampHa*1.5, ymax=ampHa*1.2)

    # Fit the [NII] and [SII] doublets and [OI], allow 2 components for Halpha
    # The widths of all narrow lines are tied to the width of [SIIb]
    # The wavelength is tied to the one of Halpha
    # The amplitude is fixed by calculation
    guesses = [ampOI, OI, sigma, ampNIIa, NIIa, sigma, ampHa, Halpha, sigma,
               ampHa/2., Halpha, 50, ampNIIb, NIIb, sigma, ampSIIa, SIIa,
               sigma, ampSIIb, SIIb, sigma]
    tied = ['', 'p[7] + {0}'.format(OI_off_Ha), 'p[20]', '',
            'p[7] + {0}'.format(NIIa-Halpha), 'p[20]', '', '', 'p[20]', '',
            'p[7]', '', '', 'p[7] + {0}'.format(NIIb-Halpha), 'p[20]', '',
            'p[7] + {0}'.format(SIIa-Halpha), 'p[20]', '',
            'p[7] + {0}'.format(SIIb-Halpha), 'p[20]']
    F, T = False, True
    # Fix the amplitude
    fixed = [T, F, F, T, F, F, T, F, F, F, F, F, T, F, F, T, F, F, T, F, F]

    # Actually do the fit, and plot the gaussians and the residual
    spec.specfit(guesses=guesses, tied=tied, annotate=False)
    spec.specfit.plot_components(add_baseline=True, component_yoffset=-ampHa*1)
    spec.specfit.plotresiduals(axis=spec.plotter.axis, 	clear=False,
                               yoffset=-ampHa*1.2, label=False)
    spec.plotter.refresh()

    # Let's use the measurements class to derive information about the emission
    # lines.  The galaxy's redshift and the flux normalization of the spectrum
    # must be supplied to convert measured fluxes to line luminosities.  If the
    # spectrum we loaded in FITS format, BUNITS would be read and we would not
    # need to supply 'fluxnorm'.
    spec.measure(z, fluxnorm=1e-17)  # Initialize the measurements class

    # Now overplot positions of lines and annotate
    y = spec.plotter.ymax * 0.85  # Location of annotations in y

    for i, line in enumerate(spec.measurements.lines.keys()):

        # If this line is not in our database of lines, don't try to annotate
        if line not in spec.speclines.optical.lines.keys():
            continue

        x = spec.measurements.lines[line]['modelpars'][1]
        # Draw dashed line to mark its position
        spec.plotter.axis.plot([x]*2, [spec.plotter.ymin, spec.plotter.ymax],
                               ls='--', color='k')
        spec.plotter.axis.annotate(spec.speclines.optical.lines[line][-1],
                                   (x, y), rotation=90, ha='right',
                                   va='center')  # Label it

    # Make some nice axis labels
    spec.plotter.axis.set_xlabel(r'Wavelength $(\AA)$')
    spec.plotter.axis.set_ylabel(r'Flux $(10^{-17} \mathrm{erg/s/cm^2/\AA})$')
    spec.plotter.refresh()

    # Save the figure
    name = splitext(basename(file))[0]+'Ha'
    spec.plotter.figure.savefig(name)

    # Print out spectral line information
    print '\n', file
    print ("Line   Flux (erg/s/cm^2)     Amplitude (erg/s/cm^2)    "
           "FWHM (Angstrom)   Luminosity (erg/s)")
    for line in spec.measurements.lines.keys():
        print line, spec.measurements.lines[line]['flux'],\
            spec.measurements.lines[line]['amp'],\
            spec.measurements.lines[line]['fwhm'],\
            spec.measurements.lines[line]['lum']

    # Had we not supplied the objects redshift (or distance), the line
    # luminosities would not have been measured, but integrated fluxes would
    # still be derived.  Also, the measurements class separates the broad and
    # narrow H-alpha components, and identifies which lines are which. Nice!

    spec.specfit.plot_fit()

    return spec


# Function to fit the OIII/Hbeta lines


def fit_OIII(file, extinction, z, sigma):
    # Initialize spectrum object and plot region surrounding Hbeta\OIIb complex
    spec = loadsdss(file, z)
    #spec = loadtrump(file, z)
    #spec = loadkitt(file,z)

    galextinction(extinction, spec)  # Take galactic extinction into account
    print file

    spec.plotter(xmin=Hbeta - 20, xmax=OIIIb + 20)

    spec.baseline(xmin=4840, xmax=5050, exclude=[4850, 4870, 4940, 4970,
                  4995, 5020, 6450, 6746, 6815, 6884, 7003, 7126, 7506,
                  7674, 8142, 8231], subtract=True, reset_selection=True,
                  highlight_fitregion=True, order=1)

    #Calculate the amplitude of OIII and Hbeta
    spec.specfit.selectregion(xmin=Hbeta - 10, xmax=Hbeta + 10)
    ampHb = np.max(spec.data[spec.specfit.xmin:spec.specfit.xmax])

    spec.specfit.selectregion(xmin=OIIIa - 10, xmax=OIIIa + 10)
    ampOIIIa = np.max(spec.data[spec.specfit.xmin:spec.specfit.xmax])

    spec.specfit.selectregion(xmin=OIIIb - 10, xmax=OIIIb + 10)
    ampOIIIb = np.max(spec.data[spec.specfit.xmin:spec.specfit.xmax])

    # Initialize spectrum and plot region surrounding Hbeta-[OIII] complex
    spec.plotter(xmin=Hbeta - 20, xmax=OIIIb + 20, ymin=-ampOIIIb*1.5,
                 ymax=ampOIIIb*1.2)

    # The widths of all narrow lines are tied to the width of [OIIIb]
    # The Wavelength of all is tied with the wavelength of [OIIb]
    # The amplitude is calculated and then fixed
    guesses = [ampHb, Hbeta, sigma, ampOIIIa, OIIIa, sigma, ampOIIIb,
               OIIIb, sigma]
    tied = ['', 'p[7] - {0}'.format(OIIIb-Hbeta), 'p[8]', '',
            'p[7] - {0}'.format(OIIIb-OIIIa), 'p[8]', '', '', '']
    F, T = False, True
    fixed = [T, F, F, T, F, F, T, F, F]

    # Actually do the fit.
    spec.specfit(guesses=guesses, tied=tied, fixed=fixed, annotate=False,
                 exclude=[3275.6, 4700, 5070, 9619.5])
    spec.specfit.plot_components(add_baseline=True,
                                 component_yoffset=-0.8*ampOIIIb)
    spec.specfit.plotresiduals(axis=spec.plotter.axis, clear=False,
                               yoffset=-1.2*ampOIIIb, label=False)
    spec.plotter.refresh()

    # Let's use the measurements class to derive information about the emission
    # lines.  The galaxy's redshift and the flux normalization of the spectrum
    # must be supplied to convert measured fluxes to line luminosities.  If the
    # spectrum we loaded in FITS format, BUNITS would be read and we would not
    # need to supply 'fluxnorm'.
    spec.measure(z, fluxnorm=1e-17)  # Initialize the measurements class

    # Now overplot positions of lines and annotate
    y = spec.plotter.ymax * 0.85    # Location of annotations in y

    for i, line in enumerate(spec.measurements.lines.keys()):

        # If this line is not in our database, don't try to annotate it
        if line not in spec.speclines.optical.lines.keys():
            continue

        #Location of emission line
        x = spec.measurements.lines[line]['modelpars'][1]
        # Draw dashed line to mark its position
        spec.plotter.axis.plot([x]*2, [spec.plotter.ymin, spec.plotter.ymax],
                               ls='--', color='k')
        #Label it
        spec.plotter.axis.annotate(spec.speclines.optical.lines[line][-1],
                                   (x, y), rotation=90, ha='right',
                                   va='center')

    #Make some nice axis labels
    spec.plotter.axis.set_xlabel(r'Wavelength $(\AA)$')
    spec.plotter.axis.set_ylabel(r'Flux $(10^{-17} \mathrm{erg/s/cm^2/\AA})$')
    spec.plotter.refresh()

    # Save the figure
    name = splitext(basename(file))[0]+'Hb'
    spec.plotter.figure.savefig(name)

    # Print out spectral line information
    print ("Line   Flux (erg/s/cm^2)     Amplitude (erg/s/cm^2)    "
           "FWHM (Angstrom)   Luminosity (erg/s)")
    for line in spec.measurements.lines.keys():
        print line, spec.measurements.lines[line]['flux'],\
            spec.measurements.lines[line]['amp'],\
            spec.measurements.lines[line]['fwhm'],\
            spec.measurements.lines[line]['lum']

    # Had we not supplied the objects redshift (or distance), the line
    # luminosities would not have been measured, but integrated fluxes would
    # still be derived.  Also, the measurements class separates the broad and
    # narrow H-alpha components, and identifies which lines are which. Nice !

    return spec


# Function that return all infos about the fit of Halpha region


def grablinesha(spec):
# We check first if Halpha is usable :
    try:
        Halphaf = spec.measurements.lines['H_alpha']['flux']
    except KeyError:
        Halphaf = 1e-50

    if Halphaf <= 1e-50:
        fluxrat = [1e-10, 1e-10, 1e-10]
        Halphaf = 1e-50
        fluxraterr = [0, 0, 0]
        fluxline = [0, 0, 0, 0, 0, 0, 0]
        broad = False

    else:
        # Now we control the error on Halpha
        try:
            errHalpha = (spec.measurements.lines['H_alpha']['modelerrs'][0] /
                         spec.measurements.lines['H_alpha']['modelpars'][0] +
                         0.03)
            print 'errHalpha : ', errHalpha
        except KeyError:
            errHalpha = 0

        if errHalpha >= 1 or errHalpha <= 0:
            fluxrat = [1e-10, 1e-10, 1e-10]
            fluxraterr = [0, 0, 0]
            fluxline = [0, 0, 0, 0, 0, 0, 0]
            broad = False
        else:
            # Now we calculate the rest of the emission line fluxes
            try:
                OIf = spec.measurements.lines['OI']['flux']
            except KeyError:
                OIf = 1e-50
            if OIf <= 1e-50:
                OIf = 1e-50

            try:
                SIIaf = spec.measurements.lines['SIIa']['flux']
            except KeyError:
                SIIaf = 1e-50
            if SIIaf <= 1e-50:
                SIIaf = 1e-50

            try:
                SIIbf = spec.measurements.lines['SIIb']['flux']
            except KeyError:
                SIIbf = 1e-50
            if SIIbf <= 1e-50:
                SIIbf = 1e-50

            try:
                NIIaf = spec.measurements.lines['NIIa']['flux']
            except KeyError:
                NIIaf = 1e-50
            if NIIaf <= 1e-50:
                NIIaf = 1e-50

            try:
                NIIbf = spec.measurements.lines['NIIb']['flux']
            except KeyError:
                NIIbf = 1e-50
            if NIIbf <= 1e-50:
                NIIbf = 1e-50

            SIIf = SIIaf+SIIbf
            fluxrat = [OIf/Halphaf, NIIbf/Halphaf, SIIf/Halphaf]

            try:
                Halpha1f = spec.measurements.lines['H_alpha_1']['flux']
                FWHM = spec.measurements.lines['H_alpha_1']['fwhm']
            except KeyError:
                Halpha1f = 1e-50
                FWHM = 0

            if Halpha1f <= 1e-50:
                Halpha1f = 1e-50
                FWHM = 0
                broad = False
            else:
                if FWHM > 21.897:
                    broad = True
                else:
                    broad = False

            # Now we look at the errors on the emission lines
            lines = ['OI', 'SIIa', 'NIIb']
            pererrs = []

            for i in lines:
                try:
                    pererr = (spec.measurements.lines[i]['modelerrs'][0] /
                              spec.measurements.lines[i]['modelpars'][0] +
                              errHalpha)  # PROBLEM
                    if i == 'SIIa':
                        j = 'SIIb'
                        pererr += (spec.measurements.lines[j]['modelerrs'][0] /
                                   spec.measurements.lines[j]['modelpars'][0])

                    pererrs.append(pererr)
                except KeyError:
                    pererrs.append(0)
            # If the error is higher than flux or negative,
            # then error=0 and fluxrate also.
            # Goes away from the diagnotics plots
            if pererrs[0] >= 1 or pererrs[0] <= 0:
                pererrs[0] = 0
                OIf = 1e-50
                fluxrat[0] = 1e-10
            if pererrs[1] >= 1 or pererrs[1] <= 0:
                pererrs[1] = 0
                SIIf = 1e-50
                fluxrat[1] = 1e-10
            if pererrs[2] > 1 or pererrs[2] < 0:
                pererrs[2] = 0
                NIIbf = 1e-50
                fluxrat[2] = 1e-10
            fluxraterr = [OIf/Halphaf*pererrs[0], NIIbf/Halphaf*pererrs[1],
                          SIIf/Halphaf*pererrs[2]]
            fluxline = [OIf, NIIaf, Halphaf, Halpha1f, NIIbf, SIIaf, SIIbf]

    return fluxrat, fluxraterr, fluxline, broad


def grablineshb(spec):

    # We check first if Hbeta is usable :
    try:
        Hbetaf = spec.measurements.lines['H_beta']['flux']

    except KeyError:
        Hbetaf = 1e-50

    if Hbetaf <= 1e-50:
        fluxrat = 1e-10
        Hbetaf = 1e-50
        fluxraterr = 0
        fluxline = [0, 0, 0]

    else:
        # Now we control the error on Hbeta
        try:
            # modelerrs always 0 !! modelpars : [amplitude,wavelength, sigma]
            errHbeta = (spec.measurements.lines['H_beta']['modelerrs'][0] /
                        spec.measurements.lines['H_beta']['modelpars'][0]+0.03)
            print 'error Hbeta : ', errHbeta

        except KeyError:
            errHbeta = 0

        if errHbeta >= 1 or errHbeta <= 0:
            fluxrat = 1e-10
            fluxraterr = 0
            fluxline = [0, 0, 0]

        else:
            # Now we calculate the rest of the emission line fluxes
            try:
                OIIIaf = spec.measurements.lines['OIIIa']['flux']
            except KeyError:
                OIIIaf = 1e-50
            if OIIIaf <= 1e-50:
                OIIIaf = 1e-50

            try:
                OIIIbf = spec.measurements.lines['OIIIb']['flux']
            except KeyError:
                OIIIbf = 1e-50
            if OIIIbf <= 1e-50:
                OIIIbf = 1e-50

            fluxrat = OIIIbf/Hbetaf
            pererrs = []
            lines = ['OIIIa', 'OIIIb']
            for i in lines:
                try:
                    pererr = (spec.measurements.lines[i]['modelerrs'][0] /
                              spec.measurements.lines[i]['modelpars'][0] +
                              errHbeta+0.03)
                    pererrs.append(pererr)
                except KeyError:
                    pererrs.append(0)
            # If the error is higher than flux or negative,
            # then error=0 and fluxrate also.
            # Goes away from the diagnotics plots
            if pererrs[0] >= 1 or pererrs[0] <= 0:
                pererrs[0] = 0
                OIIIaf = 1e-50
            if pererrs[1] >= 1 or pererrs[1] <= 0:
                pererrs[1] = 0
                OIIIbf = 1e-50
                fluxrat = 1e-10

            fluxraterr = [OIIIbf/Hbetaf*pererrs[1]]
            fluxline = [Hbetaf, OIIIaf, OIIIbf]

    return fluxrat, fluxraterr, fluxline


def plotformat():
    figure()
    axes([0.15, 0.15, 0.95-0.15, 0.95-0.2])  # [left, bottom, width, height]
    setp(gca().get_ymajorticklabels(), fontsize='large')
    setp(gca().get_xmajorticklabels(), fontsize='large')
    return


def n2haplot(x, y, xerr, yerr, label, color):
    plotformat()
    xlabel(r'$\log$ [NII]/H$\alpha$', size=28)
    ylabel(r'$\log$ [OIII]/H$\beta$', size=28)

    n2ha = arange(-1.25, 0, 0.05)
    o3hb = 0.61/((n2ha-0.05))+1.3
    n2ha2 = arange(-2, 0.4, 0.05)
    o3hb2 = 0.61/((n2ha2-0.47))+1.19
    n2ha3 = arange(-1.5, 4)
    o3hb3 = 0.7 - 1.2 * (n2ha3 + 0.4)
    n2ha4 = arange(-1.5, 4)
    o3hb4 = 0.7 - 1.2 * (n2ha4 - 0.4)

    ylim(-1.25, 1.5)
    xlim(-2, 1)

    plot(n2ha, o3hb, ':', color='black')
    plot(n2ha2, o3hb2, color='black')

    # plot(n2ha3,o3hb3,color='blue')
    # plot(n2ha4,o3hb4,color='red')
    # print x,y
    # plot(x,y,'*')

    for i in range(0, len(x)):
        errorbar(math.log10(x[i]), math.log10(y[i]),
                 xerr=[[math.log10(x[i])-math.log10(x[i]-xerr[i])],
                 [math.log10(xerr[i]+x[i])-math.log10(x[i])]],
                 yerr=[[math.log10(y[i])-math.log10(y[i]-yerr[i])],
                 [math.log10(yerr[i]+y[i])-math.log10(y[i])]], marker='o',
                 color=color[i], ms=5, label=label[i])
    figtext(0.3, 0.4, 'HII', size=20)
    figtext(0.75, 0.5, 'AGN', size=20)
    figtext(0.63, 0.2, 'COMP', size=20)

    # lg=legend(loc='upper right',prop=prop,scatterpoints=1,numpoints=1)
    # lg.get_frame().set_linewidth(0)
    line1 = Line2D(range(1), range(1), color="white", marker='o',
                   markerfacecolor="red")
    line2 = Line2D(range(1), range(1), color="white", marker='o',
                   markerfacecolor="Blue")
    plt.legend((line1, line2), ('v<1000km/s', 'v>1000km/s'),
               numpoints=1, loc=3)
    savefig('n2ha.png')

    return


def s2haplot(x, y, xerr, yerr, label, color):

    plotformat()
    xlabel(r'$\log$ [SII]/H$\alpha$', size=28)
    ylabel(r'$\log$ [OIII]/H$\beta$', size=28)

    s2ha = arange(-2, 0.2, 0.05)
    o3hb = 0.72/((s2ha-0.32))+1.3
    s2ha2 = arange(-0.315, 0.25, 0.05)
    o3hb2 = 1.89*s2ha2+0.76

    ylim(-1.2, 1.5)
    xlim(-1.2, 0.8)

    plot(s2ha, o3hb, color='black')
    plot(s2ha2, o3hb2, ':', color='black')

    for i in range(0, len(x)):
        errorbar(math.log10(x[i]), math.log10(y[i]),
                 xerr=[[math.log10(x[i])-math.log10(x[i]-xerr[i])],
                 [math.log10(xerr[i]+x[i])-math.log10(x[i])]],
                 yerr=[[math.log10(y[i])-math.log10(y[i]-yerr[i])],
                 [math.log10(yerr[i]+y[i])-math.log10(y[i])]], marker='o',
                 color=color[i], ms=5, label=label[i])

    figtext(0.3, 0.4, 'HII', size=20)
    figtext(0.4, 0.75, 'Seyfert', size=20)
    figtext(0.7, 0.5, 'LINER', size=20)

    # lg=legend(loc='upper right',prop=prop,scatterpoints=1,numpoints=1)
    # lg.get_frame().set_linewidth(0)
    line1 = Line2D(range(1), range(1), color="white", marker='o',
                   markerfacecolor="red")
    line2 = Line2D(range(1), range(1), color="white", marker='o',
                   markerfacecolor="Blue")
    plt.legend((line1, line2), ('v<1000km/s', 'v>1000km/s'),
               numpoints=1, loc=3)
    savefig('s2ha.png')

    return


def o1haplot(x, y, xerr, yerr, label, color):
    plotformat()
    xlabel(r'$\log$ [OI]/H$\alpha$', size=28)
    ylabel(r'$\log$ [OIII]/H$\beta$', size=28)

    o1ha = arange(-2.5, -0.7, 0.05)
    o3hb = 0.73/((o1ha+0.59))+1.33
    o1ha2 = arange(-1.126, 0, 0.05)
    o3hb2 = 1.18*o1ha2+1.3

    ylim(-1.2, 1.5)
    xlim(-2.2, 0.0)

    plot(o1ha, o3hb, color='black')
    plot(o1ha2, o3hb2, ':', color='black')
    for i in range(0, len(x)):
        errorbar(math.log10(x[i]), math.log10(y[i]),
                 xerr=[[math.log10(x[i])-math.log10(x[i]-xerr[i])],
                 [math.log10(xerr[i]+x[i])-math.log10(x[i])]],
                 yerr=[[math.log10(y[i])-math.log10(y[i]-yerr[i])],
                 [math.log10(yerr[i]+y[i])-math.log10(y[i])]], marker='o',
                 color=color[i], ms=5, label=label[i])

    figtext(0.3, 0.4, 'HII', size=20)
    figtext(0.5, 0.7, 'Seyfert', size=20)
    figtext(0.7, 0.4, 'LINER', size=20)

    # lg=legend(loc='upper right',prop=prop,scatterpoints=1,numpoints=1)
    # lg.get_frame().set_linewidth(0)
    line1 = Line2D(range(1), range(1), color="white", marker='o',
                   markerfacecolor="red")
    line2 = Line2D(range(1), range(1), color="white", marker='o',
                   markerfacecolor="Blue")
    plt.legend((line1, line2), ('v<1000km/s', 'v>1000km/s'),
               numpoints=1, loc=3)
    savefig('o1ha.png')

    return


# MAIN..................................................................

# 1) Open the .fits file and make the fit !

# Read the data file (name, redshift and name for extinction)
sigma = 3
names = []
redshift = []
extinctionname = []
with open('data', 'r') as infile:
    csv_reader = csv.reader(infile, delimiter='\t')
    for line in csv_reader:
        names.append(line[0])
        redshift.append(line[1])
        extinctionname.append(line[2])
infile.close()
redshift = map(float, redshift)
nbr = len(redshift)  # nbr of spectrum
spec = [None]*2*nbr  # create the spec array

# Do the fit
for i in range(0, nbr):
    spec[2*i] = fit_OI(names[i], extinctionname[i], redshift[i], sigma)
    spec[2*i+1] = fit_OIII(names[i], extinctionname[i], redshift[i], sigma)

# 2) Calculate the different flux ratios and errors

# Declare the matrix
fluxrat = np.zeros((nbr, 3))
fluxraterra = np.zeros((nbr, 3))
fluxratb = np.zeros((nbr, 1))
fluxraterrb = np.zeros((nbr, 1))
fluxlineha = np.zeros((nbr, 7))
fluxlinehb = np.zeros((nbr, 3))
broad = [False] * nbr

# Create a file "flux" with the flux of every emission lines
file_flux = open("flux", "w")
file_flux.write('Name redshift OI NIIa Halpha Halpha_broad '
                'NIIb SIIa SIIb Hbeta OIIIa OIIIb broad\n')
file_flux.close()
file_flux = open("flux", "a")

# Create a file "diagnostic" with the flux ratio and the errors
file_diagnostic = open("diagnostic", "w")
file_diagnostic.write('Name redshift OI/Halpha NIIb/Halpha SII/Halpha '
                      'OIIIb/Hbeta errOI/Halpha errNIIb/Halpha errSII/Halpha '
                      'errOIIb/Hbeta\n')
file_flux.close()
file_flux = open("flux", "a")

# Write the informations inside the files
for i in range(0, len(spec)/2):
    fluxrat[i, :], fluxraterra[i, :], fluxlineha[i, :], broad[i] = \
        grablinesha(spec[2*i])
    fluxratb[i, :], fluxraterrb[i, :], fluxlinehb[i, :] = \
        grablineshb(spec[2*i+1])
    file_flux.write(str(splitext(names[i])[0])+' '+str(redshift[i])+' ')
    for j in range(0, 7):  # write all the fluxes in the file flux
        file_flux.write(str(fluxlineha[i, j])+' ')
    for j in range(0, 3):
        file_flux.write(str(fluxlinehb[i, j])+' ')
    file_flux.write(str(broad[i]))
    file_flux.write('\n')
    file_diagnostic.write(str(splitext(names[i])[0])+' '+str(redshift[i])+' ')
    for j in range(0, 3):
        file_diagnostic.write(str(fluxrat[i, j])+' ')
    file_diagnostic.write(str(fluxratb[i, 0])+' ')
    for j in range(0, 3):
        file_diagnostic.write(str(fluxraterra[i, j])+' ')
    file_diagnostic.write(str(fluxraterrb[i, 0])+'\n')

# Close the files
file_flux.close()
file_diagnostic.close()

print 'Flux rate of Halpha :\n', fluxrat, '\n'
print 'Error flux rate of Halpha :\n', fluxraterra, '\n'
print '...................................\n'
print 'Flux rate of Hbeta :\n', fluxratb, '\n'
print 'Error flux rate of Hbeta :\n', fluxraterrb, '\n'

# 4) Make the diagnostic plots

# Declare the arrays
label = names
color = ['red']*nbr
for i in range(0, nbr):
    if broad[i] is True:
        color[i] = 'blue'
    else:
        color[i] = 'red'
fluxraterr = []
fluxraterr = [0]*nbr
print '.....................................\n'

# Plot the diagnostics
o1haplot(fluxrat[:, 0], fluxratb[:, 0], fluxraterra[:, 0],
         fluxraterrb[:, 0], label, color)
s2haplot(fluxrat[:, 2], fluxratb[:, 0], fluxraterra[:, 2],
         fluxraterrb[:, 0], label, color)
n2haplot(fluxrat[:, 1], fluxratb[:, 0], fluxraterra[:, 1],
         fluxraterrb[:, 0], label, color)
