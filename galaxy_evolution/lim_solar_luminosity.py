import numpy as np

# To convert Lim+2020 (ApJ 889 80) luminosity in units of log10(L*/L_\odot) to
# AB magnitudes we need to know the AB magnitude of L_\odot, which defines M=0
# in the log-luminosity scale. Lim+2020 quote L_IR, integrating from 8 to 1000
# microns.

# Data from http://svo2.cab.inta-csic.es/svo/theory/fps3/morefiles/sun.dat
# "Sun spectrum from calspec database: sun_reference_stis_002 composed with
# sun_mod_001 for larger wavelengths."
a = np.loadtxt('data/sun.dat')

min_l, max_l = 8, 1000  # microns
q = (a[:, 0] >= min_l * 1e4) & (a[:, 0] <= max_l * 1e4)
sun_avg_flux_dens = np.sum(a[q, 1]) / ((max_l - min_l) * 1e4)
# Most of this comes from http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php:
l_array = np.linspace(min_l, max_l, 1000000) * 1e4
f0v = 3631  # janksy
ab_avg_flux_dens = np.sum(2.9979246e-5 * l_array**-2 * f0v) / (l_array[-1] - l_array[0])
# With both flux densities, we can calculate the AB apparent magnitude of the Sun
# in the 8-1000 micron "bandpass".
sun_app_mag = -2.5 * np.log10(sun_avg_flux_dens / ab_avg_flux_dens)
# And convert to absolute magnitude, using 1AU/1pc = 4.848e-6.
sun_abs_mag = sun_app_mag - 5 * np.log10(4.848e-6) + 5
print(sun_avg_flux_dens, ab_avg_flux_dens, sun_app_mag, sun_abs_mag)
