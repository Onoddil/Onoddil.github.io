# Schechter Function Parameterisation for Galaxy Counts

Here I provide the data used to construct the galaxy count model described in Wilson (2022, RNAAS, ...). The model assumes galaxy luminosities are described by a Schechter function, where the density of galaxies as a function of absolute magnitude in a given bandpass is given by
$\phi(M) = 0.4 \ln(10) \phi^* [10^{-0.4 (M-M^*)}]^{\alpha+1} \exp(-10^{-0.4 (M - M*)})$, with $M^*(z) = M^*_0 - Qz$, $\phi^* = \phi^*_0 10^{0.4 P z}$ describing the redshift dependency of the characteristic absolute magnitude and normalising density of the luminosity function respectively and $\alpha$ the slope of the faint end of the distribution.

The observed galaxy count is then derived by a conversion from $\phi$ in units of $\mathrm{Mpc}^{-3}\,\mathrm{mag}^{-1}$ to sky number densities ($\mathrm{deg}^{-2}\,\mathrm{mag}^{-1}$) by consideration of the volume within a particular redshift range and sky region, and a sum over all redshifts. Final observed densities of galaxies in a bandpass are the sum of two Schechter functions, $\phi_\mathrm{tot} = \phi_b + \phi_r$, for "blue" (star-forming) and "red" (quiescent) galaxies, and hence we require a total of 10 parameters to describe the density of galaxies at a particular redshift.

To derive galaxy counts for an arbitrary bandpass, a parameterisation for these 10 parameters -- $M^*_0$, $\phi^*_0$, $P$, $Q$, and $\alpha$, for both blue and red galaxy distributions -- was found as a function of wavelength, with the exception of $Q$, which -- due to the degeneracy between it and $P$ -- was derived as a function of $P$. The functions used to describe these parameters are $y = m \times x + c$ for $Q(P)$, $\alpha(\log_{10}(\lambda))$, and $P(\log_{10}(\lambda))$; $y = A \exp(-m x) + c$ for $M^*_0(\log_{10}(\lambda))$ and blue galaxy $\phi^*_0(\log_{10}(\lambda))$; and $y = A \exp(-0.5 (x - u)^2 \times m) + c$ for red galaxy $\phi^*_0(\log_{10}(\lambda))$. The table below provides the values derived for each parameter fit for both galaxy types, with a machine-readable version of the data available [here](galaxy_count_parameter_table.csv).

Parameter | Galaxy type | c | m | a | u
| --- | --- | --- | --- | --- | --- |
| $M^*_0$ | b | -24.2865 | 1.142 | 2.6558 |  |
| $M^*_0$ | r | -23.1925 | 1.779 | 1.6683 |  |
| $\phi^*_0$ | b | 0.0015 | 2.919 | 0.0005 |  |
| $\phi^*_0$ | r | 0.0006 | 7.691 | 0.0033 | -0.066 |
| $\alpha$ | b | -1.2578 | 0.021 |  |  |
| $\alpha$ | r | -0.3091 | -0.067 |  |  |
| $P$ | b | -0.3020 | 0.034 |  |  |
| $P$ | r | -0.7131 | 0.233 |  |  |
| $Q$ | b | 1.2574 | -0.390 |  |  |
| $Q$ | r | 0.9125 | -0.549 |  |  |

Additionally, the final table provides the values used to calculate these functional forms, also available in machine-readable form [here]. The machine-readable .csv file also provides the original, un-corrected literature values of $M^*_0$ and $\phi^*_0$ -- in linear- and log-form where quoted solely in $\log_{10}(\phi^*)$ -- and $h$ and $z_0$ when conversion from non-zero redshift or values in arbitrary $H_0$ values are quoted.

| Citation | Band | Wavelength (nm) | Type | $M^*_0$ / AB mag | $\phi^*_0$ / Mpc$^{-3}$mag$^{-1}$ | $\alpha$ | Q / mag $z^{-1}$ | P / $z^{-1}$ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| MNRAS 420 1239 | u | 355.0 | b |  -18.495±0.117 | 4.493e-03±8.575e-04 | -1.430±0.070 | 5.500±0.600 | -7.100±1.500 |
| MNRAS 420 1239 | u | 355.0 | r |  -17.475±0.157 | 1.468e-02±4.939e-03 | -0.140±0.130 | 6.400±1.400 | -8.100±3.400 |
| MNRAS 420 1239 | g | 469.0 | b |  -20.045±0.092 | 2.504e-03±1.029e-04 | -1.400±0.030 | 3.100±0.700 | -1.200±1.500 |
| MNRAS 420 1239 | g | 469.0 | r |  -19.725±0.152 | 4.322e-03±5.831e-04 | -0.430±0.050 | 3.600±1.400 | -3.900±2.800 |
| MNRAS 420 1239 | r | 617.0 | b |  -20.975±0.067 | 1.303e-03±1.715e-04 | -1.490±0.030 | 0.800±0.300 | 2.900±0.600 |
| MNRAS 420 1239 | r | 617.0 | r |  -20.875±0.032 | 3.807e-03±5.145e-04 | -0.570±0.030 | 1.800±0.100 | -1.200±0.500 |
| MNRAS 420 1239 | i | 748.0 | b |  -21.285±0.072 | 1.441e-03±2.058e-04 | -1.450±0.020 | 1.700±0.400 | 1.200±0.900 |
| MNRAS 420 1239 | i | 748.0 | r |  -21.205±0.032 | 3.979e-03±5.145e-04 | -0.540±0.030 | 2.000±0.100 | -1.800±0.500 |
| MNRAS 420 1239 | z | 892.0 | b |  -21.675±0.045 | 1.166e-03±1.372e-04 | -1.450±0.030 | 0.900±0.200 | 3.600±0.500 |
| MNRAS 420 1239 | z | 892.0 | r |  -21.405±0.067 | 4.528e-03±5.488e-04 | -0.490±0.050 | 2.400±0.300 | -2.700±0.700 |
| MNRAS 451 1540 | r | 617.0 | b |  -21.265±0.051 | 1.842e-03±2.121e-04 | -1.380±0.060 | 1.090±0.100 | 1.300±0.250 |
| MNRAS 451 1540 | r | 617.0 | r |  -21.610±0.063 | 2.020e-03±2.325e-04 | -0.790±0.110 | 0.580±0.180 | 1.550±0.400 |
| MNRAS 439 1245 | u | 355.0 | b |  -18.530±1.250 | 9.640e-03±1.290e-02 | -0.910±7.170 | | |
| MNRAS 439 1245 | u | 355.0 | r |  -18.530±1.250 | 1.460e-03±9.390e-03 | 1.250±6.010 | | |
| MNRAS 439 1245 | g | 469.0 | b |  -20.280±0.260 | 3.510e-03±1.750e-03 | -1.290±0.140 | | |
| MNRAS 439 1245 | g | 469.0 | r |  -20.280±0.260 | 4.880e-03±1.220e-03 | 0.060±0.590 | | |
| MNRAS 439 1245 | r | 617.0 | b |  -20.900±0.260 | 4.510e-03±1.030e-03 | -1.130±0.070 | | |
| MNRAS 439 1245 | r | 617.0 | r |  -20.900±0.260 | 3.010e-03±8.300e-04 | 0.530±0.560 | | |
| MNRAS 439 1245 | i | 748.0 | b |  -21.450±0.200 | 2.200e-03±1.500e-03 | -1.350±0.210 | | |
| MNRAS 439 1245 | i | 748.0 | r |  -21.450±0.200 | 4.870e-03±1.340e-03 | -0.090±0.480 | | |
| MNRAS 439 1245 | z | 892.0 | b |  -21.780±0.250 | 1.400e-03±1.120e-03 | -1.460±0.210 | | |
| MNRAS 439 1245 | z | 892.0 | r |  -21.780±0.250 | 5.050e-03±8.800e-04 | -0.260±0.440 | | |
| MNRAS 439 1245 | Y | 1031.0 | b |  -21.760±0.240 | 1.440e-03±1.250e-03 | -1.450±0.230 | | |
| MNRAS 439 1245 | Y | 1031.0 | r |  -21.760±0.240 | 4.830e-03±8.400e-04 | -0.100±0.540 | | |
| MNRAS 439 1245 | J | 1248.0 | b |  -21.820±0.170 | 1.580e-03±7.600e-04 | -1.380±2.430 | | |
| MNRAS 439 1245 | J | 1248.0 | r |  -21.820±0.170 | 4.780e-03±7.200e-04 | 0.080±2.580 | | |
| MNRAS 439 1245 | H | 1631.0 | b |  -22.040±0.260 | 1.350e-03±6.330e-03 | -1.460±2.430 | | |
| MNRAS 439 1245 | H | 1631.0 | r |  -22.040±0.260 | 5.300e-03±6.690e-03 | 0.080±2.580 | | |
| MNRAS 439 1245 | K | 2201.0 | b |  -21.720±0.230 | 1.640e-03±3.130e-03 | -1.390±1.620 | | |
| MNRAS 439 1245 | K | 2201.0 | r |  -21.720±0.230 | 5.090e-03±3.180e-03 | 0.240±1.550 | | |
| ApJ 518 533L | U | 366.0 | b |  -19.942±0.248 | 3.259e-03±1.338e-03 | -1.140±0.130 | 0.510±0.660 | 2.670±0.920 |
| ApJ 518 533L | U | 366.0 | r |  -19.754±0.264 | 3.156e-03±8.918e-04 | -0.510±0.150 | 0.970±0.700 | 0.680±0.690 |
| ApJ 518 533L | B | 435.0 | b |  -19.981±0.266 | 2.470e-03±1.132e-03 | -1.230±0.120 | 0.180±0.710 | 3.080±0.990 |
| ApJ 518 533L | B | 435.0 | r |  -19.361±0.190 | 6.963e-03±1.235e-02 | 0.080±0.140 | 1.580±0.490 | -1.070±0.490 |
| ApJ 518 533L | Rc | 710.0 | b |  -20.852±0.286 | 1.921e-03±1.029e-03 | -1.340±0.120 | 0.110±0.740 | 3.170±1.030 |
| ApJ 518 533L | Rc | 710.0 | r |  -21.038±0.284 | 2.744e-03±7.889e-04 | -0.630±0.150 | 0.690±0.760 | 0.890±0.740 |
| ApJ 697 506 | [3.6] | 3544.0 | b |  -22.009±0.266 | 2.573e-03±6.860e-05 | -1.400±0.180 | 1.000±0.700 | 0.000±0.000 |
| ApJ 697 506 | [3.6] | 3544.0 | r |  -21.819±0.253 | 1.509e-03±1.715e-04 | -0.630±0.290 | 1.400±0.500 | 0.000±0.000 |
| ApJ 697 506 | [4.5] | 4487.0 | b |  -21.551±0.273 | 2.470e-03±6.860e-05 | -1.290±0.170 | 0.900±0.700 | 0.000±0.000 |
| ApJ 697 506 | [4.5] | 4487.0 | r |  -21.271±0.280 | 1.269e-03±1.372e-04 | -0.600±0.280 | 1.300±0.500 | 0.000±0.000 |
| ApJ 697 506 | [5.8] | 5710.0 | b |  -22.884±0.266 | 1.406e-03±3.430e-05 | -1.630±0.150 | 0.400±0.700 | 0.000±0.000 |
| ApJ 697 506 | [5.8] | 5710.0 | r |  -21.584±0.336 | 9.947e-04±1.029e-04 | -1.330±0.280 | 1.200±0.800 | 0.000±0.000 |
| ApJ 697 506 | [8.0] | 7841.0 | b |  -23.818±0.262 | 1.509e-03±3.430e-05 | -1.350±0.090 | 1.700±0.800 | 0.000±0.000 |
| ApJ 697 506 | [8.0] | 7841.0 | r |  -22.793±1.140 | 2.401e-04±3.430e-05 | -2.030±0.470 | 1.800±3.500 | 0.000±0.000 |
| ApJ 560 566 | Ks | 2159.0 | b |  -21.920±0.060 | 3.464e-03±4.459e-04 | -0.870±0.090 | | |
| ApJ 560 566 | Ks | 2159.0 | r |  -22.470±0.060 | 1.543e-03±2.058e-04 | -0.920±0.100 | | |
| MNRAS 427 3244 | FUV | 152.0 | b |  -17.905±0.045 | | -1.140±0.015 | | |
| MNRAS 427 3244 | FUV | 152.0 | r |  -16.975±0.200 | | -0.700±0.195 | | |
| MNRAS 427 3244 | NUV | 227.0 | b |  -18.335±0.040 | | -1.160±0.015 | | |
| MNRAS 427 3244 | NUV | 227.0 | r |  -17.355±0.140 | | -0.900±0.115 | | |
| MNRAS 427 3244 | u | 355.0 | b |  -19.465±0.035 | | -1.140±0.015 | | |
| MNRAS 427 3244 | u | 355.0 | r |  -18.635±0.060 | | -0.040±0.050 | | |
| MNRAS 427 3244 | g | 469.0 | b |  -20.805±0.035 | | -1.200±0.025 | | |
| MNRAS 427 3244 | g | 469.0 | r |  -20.275±0.060 | | -0.070±0.045 | | |
| MNRAS 427 3244 | r | 617.0 | b |  -21.455±0.030 | | -1.200±0.010 | | |
| MNRAS 427 3244 | r | 617.0 | r |  -21.045±0.065 | | -0.100±0.045 | | |
| MNRAS 427 3244 | i | 748.0 | b |  -21.935±0.040 | | -1.280±0.010 | | |
| MNRAS 427 3244 | i | 748.0 | r |  -21.485±0.060 | | -0.190±0.040 | | |
| MNRAS 427 3244 | z | 892.0 | b |  -22.145±0.035 | | -1.250±0.010 | | |
| MNRAS 427 3244 | z | 892.0 | r |  -21.715±0.060 | | -0.170±0.040 | | |
| MNRAS 427 3244 | Y | 1031.0 | b |  -22.315±0.040 | | -1.250±0.010 | | |
| MNRAS 427 3244 | Y | 1031.0 | r |  -21.875±0.060 | | -0.250±0.040 | | |
| MNRAS 427 3244 | J | 1248.0 | b |  -22.445±0.040 | | -1.240±0.015 | | |
| MNRAS 427 3244 | J | 1248.0 | r |  -22.015±0.055 | | -0.270±0.040 | | |
| MNRAS 427 3244 | H | 1631.0 | b |  -22.695±0.040 | | -1.210±0.015 | | |
| MNRAS 427 3244 | H | 1631.0 | r |  -22.315±0.055 | | -0.280±0.035 | | |
| MNRAS 427 3244 | K | 2201.0 | b |  -22.325±0.040 | | -1.170±0.015 | | |
| MNRAS 427 3244 | K | 2201.0 | r |  -22.035±0.065 | | -0.310±0.045 | | |
| MNRAS 465 672 | K | 2201.0 | b |  -23.590±0.510 | 4.270e-04±1.133e-03 | -1.516±0.230 | 0.249±0.340 | -0.454±0.880 |
| MNRAS 465 672 | K | 2201.0 | r |  -23.590±0.510 | 1.386e-03±1.228e-03 | -0.817±0.740 | 0.249±0.340 | -0.725±0.410 |
| A&A 476 137 | K | 2201.0 | b |  -23.379±0.055 | 1.195e-03±1.029e-03 | -1.300±0.053 | 0.493±0.044 | -0.543±0.073 |
| A&A 476 137 | K | 2201.0 | r |  -23.161±0.450 | 1.651e-03±1.102e-03 | -0.206±0.610 | 0.453±0.400 | -1.859±0.290 |
| MNRAS 380 585 | K | 2201.0 | b |  -23.167±0.047 | 6.637e-04±4.065e-04 | -1.232±0.032 | 0.677±0.039 | -0.500±0.074 |
| MNRAS 380 585 | K | 2201.0 | r |  -22.921±0.034 | 1.322e-03±6.613e-04 | -0.112±0.061 | 0.416±0.031 | -1.049±0.074 |
| IAUS 306 40 | r | 617.0 | b |  -21.305±0.080 | 1.334e-03±1.844e-04 | -1.470±0.060 | 0.580±0.050 | 2.740±0.250 |
| IAUS 306 40 | r | 617.0 | r |  -21.345±0.070 | 2.428e-03±2.796e-04 | -0.710±0.140 | 0.790±0.100 | 1.140±0.250 |
| A&A 508 1217 | B | 435.0 | b |  -19.961±0.600 | 4.527e-03±3.120e-03 | -1.143±0.360 | | |
| A&A 508 1217 | B | 435.0 | r |  -20.500±0.560 | 3.634e-03±1.323e-03 | -0.618±0.410 | | |
| ApJ 647 853 | B | 435.0 | b |  -20.113±0.015 | 3.905e-03±3.820e-05 | -1.300±0.500 | 1.232±0.017 | -0.456±0.012 |
| ApJ 647 853 | B | 435.0 | r |  -20.799±0.023 | 2.220e-03±4.020e-05 | -0.500±0.500 | 0.459±0.027 | -0.900±0.026 |
| ApJ 748 10 | r | 617.0 | b |  -20.919±0.120 | 3.958e-03±6.339e-04 | -1.110±0.500 | 1.660±0.090 | -0.380±0.210 |
| ApJ 748 10 | r | 617.0 | r |  -21.092±0.070 | 3.184e-03±4.777e-04 | -0.550±0.500 | 1.730±0.070 | -0.950±0.100 |
| ApJ 873 78 | K | 2201.0 | b |  -22.572±0.023 | 2.339e-03±2.090e-05 | -1.200±0.500 | 0.601±0.029 | -0.652±0.014 |
| ApJ 873 78 | K | 2201.0 | r |  -22.792±0.019 | 2.906e-03±5.670e-05 | -0.500±0.500 | 0.192±0.022 | -0.459±0.029 |
| ApJ 889 80 | IR | 504000.0 | b |  -20.774±0.304 | 9.528e-04±3.281e-04 | -0.540±0.241 | 0.336±0.125 | 0.046±0.156 |
| ApJ 889 80 | IR | 504000.0 | r |  -20.774±0.304 | 6.352e-04±2.187e-04 | -0.540±0.241 | 0.336±0.125 | 0.046±0.156 |
| MNRAS 494 1894 | FUV | 154.6 | b |  -18.329±0.120 | 5.861e-03±1.376e-03 | -1.391±0.080 | 0.886±0.008 | -0.405±0.010 |
| MNRAS 494 1894 | NUV | 234.5 | b |  -18.475±0.090 | 6.067e-03±1.305e-03 | -1.367±0.060 | 0.975±0.006 | -0.419±0.008 |


### Notes
Where applicable, relevant information from each literature reference used to calculate the quoted values in the above table is given below.

- In all cases, quoted uncertainties are statistical. These are inflated by adding 5% relative uncertainty and (0, 0.001, 0.1, 0.1, 0.01) systematic uncertainty to ($M^*_0$, $\phi^*_0$, $Q$, $P$, $\alpha$) respectively during the parameterisation fitting process.

Loveday et al. (2012, MNRAS, 420, 1239)  

- Table 5
- Conversion from units of $h^3\,\mathrm{Mpc}^{-3}$ and $M^* - 5\log_{10}(h)$ using $h=0.7$
- $M^*_0$ translated from $z_0 = 0.1$ to $z=0$

Loveday et al. (2015, MNRAS, 451, 1540)

- Tables 3 and 4, Mean Probability Petrosian
- Conversion from units of $h^3\,\mathrm{Mpc}^{-3}$ and $M^* - 5\log_{10}(h)$ using $h=0.7$
- $M^*_0$ translated from $z_0 = 0.1$ to $z=0$

Kelvin et al. (2014, MNRAS, 439, 1245)

- Table 9
- $H_0 = 70\,\mathrm{km}\,\mathrm{s}^{-1}\,\mathrm{kpc}$ gives $h=1$ so no conversion of values necessary

Lin et al. (1999, ApJ, 518, 533)

- Tables 1 and 2, $Q_0 = 0.5$
- Conversion from units of $h^3\,\mathrm{Mpc}^{-3}$ and $M^* - 5\log_{10}(h)$ using $h=0.7$
- $M^*_0$ translated from $z_0 = 0.3$ to $z=0$

Dai et al. (2009, ApJ, 697, 506)

- Table 2, using "early" and "late" type
- Conversion from units of $h^3\,\mathrm{Mpc}^{-3}$ and $M^* - 5\log_{10}(h)$ using $h=0.7$
- $M^*_0$ translated from $z_0 = 0.25$ to $z=0$
- $P$ and $Q$ not used in fitting $P(\log_{10}(\lambda))$ or $Q(P)$, but has its $Q$ values represented in Figure 1 of Wilson (2022, RNAAS, ...) at $P = 0$, the value fixed during the derivation for these parameters
- $M^*$ quoted in the Vega magnitude system, so for self-consistency had its values converted to the AB system using the [MIST zero-point table](https://waps.cfa.harvard.edu/MIST/BC_tables/zeropoints.txt) values

Kochanek et al. (2001, ApJ 560 566)

- Table 3
- Conversion from units of $h^3\,\mathrm{Mpc}^{-3}$ and $M^* - 5\log_{10}(h)$ using $h=0.7$
- $M^*$ quoted in the Vega magnitude system, so for self-consistency had its values converted to the AB system using the [MIST zero-point table](https://waps.cfa.harvard.edu/MIST/BC_tables/zeropoints.txt) values

Driver et al. (2012, MNRAS, 427, 3244)

- Tables 5 and 6
- Conversion from units of $M^* - 5\log_{10}(h)$ using $h=0.7$
- Due to an apparent bias towards assigning objects as "blue" galaxies, do not use $\phi^*$ values; if the values of $\phi^*$ were used, units are $(0.5\,\mathrm{mag})^{-1}$ and hence would need correcting by a factor two

Mortlock et al. (2017, MNRAS, 465, 672)

- Table 4
- Conversion from units of $h^3\,\mathrm{Mpc}^{-3}$ and $M^* - 5\log_{10}(h)$ using $h=0.7$
- Parameters quoted as $M^*$, $\phi^*$, $\alpha$ for various $z$ values; linear slopes fit for using the original bin values
- Both red and blue parameters have significant correlation and hence overall uncertainties underestimated, statistical uncertainty padded from least-squared fit quoted values

Arnouts et al. (2007, A&A, 476, 137)

- Table 1
- Conversion from units of $h^3\,\mathrm{Mpc}^{-3}$ and $M^* - 5\log_{10}(h)$ using $h=0.7$
- Parameters quoted as $M^*$, $\phi^*$, $\alpha$ for various $z$ values; linear slopes fit for using the original bin values
- Red parameters have significant correlation and hence overall uncertainties underestimated, statistical uncertainty padded from least-squared fit quoted values

Cirasuolo et al. (2007, MNRAS, 380, 585)

- Table 2
- Conversion from units of $h^3\,\mathrm{Mpc}^{-3}$ and $M^* - 5\log_{10}(h)$ using $h=0.7$
- Parameters quoted as $M^*$, $\phi^*$, $\alpha$ for various $z$ values; linear slopes fit for using the original bin values

Loveday (2014, IAUS, 306, 40)

- Table 1
- Conversion from units of $h^3\,\mathrm{Mpc}^{-3}$ and $M^* - 5\log_{10}(h)$ using $h=0.7$
- Uncertainties in $P$ and $Q$ estimated by eye from two-dimensional $\chi^2$ contours

Zucca et al. (2009, A&A, 508, 1217)

- Table 1, "early" and "spiral" morphology types
- $H_0 = 70\,\mathrm{km}\,\mathrm{s}^{-1}\,\mathrm{kpc}$ gives $h_{70}=1$ so no conversion of values necessary
- Parameters quoted as $M^*$, $\phi^*$, $\alpha$ for various $z$ values; linear slopes fit for using the original bin values
- Both red and blue parameters have significant correlation and hence overall uncertainties underestimated, statistical uncertainty padded from least-squared fit quoted values

Willmer et al. (2006, ApJ, 647, 853)

- Tables 4 and 5, "minimal" weights
- $H_0 = 70\,\mathrm{km}\,\mathrm{s}^{-1}\,\mathrm{kpc}$ gives $h=1$ so no conversion of values necessary
- Uncertainties on $\alpha$ set sufficiently large due to fixed value
- Parameters quoted as $M^*$, $\phi^*$, $\alpha$ for various $z$ values; linear slopes fit for using the original bin values

Cool et al. (2012, ApJ, 748, 10)

- Tables 4 and 5
- Conversion from units of $h^3\,\mathrm{Mpc}^{-3}$ and $M^* - 5\log_{10}(h)$ using $h=0.7$
- $M^*_0$ translated from $z_0 = 0.1$ to $z=0$
- Uncertainties on $\alpha$ set sufficiently large due to fixed value

Beare et al. (2019, ApJ, 873, 78)

- Table 4
- $H_0 = 70\,\mathrm{km}\,\mathrm{s}^{-1}\,\mathrm{kpc}$ gives $h_{70}=1$ so no conversion of values necessary
- Uncertainties on $\alpha$ set sufficiently large due to fixed value
- Parameters quoted as $M^*$, $\phi^*$, $\alpha$ for various $z$ values; linear slopes fit for using the original bin values

Lim et al. (2020, ApJ, 889, 80)

- Table 9, non-fixed faint-end slope
- Single Schechter slope, but split 60%-40% blue-red galaxies in $\phi^*$
- Parameters quoted as $\log_{10}(L^*/L_\odot)$, $\phi^*$, $\alpha$ for various $z$ values; linear slopes fit for using the original bin values
- Conversion from "IR" luminosity $\log_{10}(L*/L_\odot)$ in 8-1000$\mu\,\mathrm{m}$ "bandpass" to absolute AB magnitudes using the Solar spectrum from $\texttt{calspec}$; value not used in derivation of $M^*_0(\log_{10}(\lambda))$ but presented for reference in the sub-mm wavelength range

Moutard et al. (2020, MNRAS, 494, 1894)

- Table 2
- $H_0 = 70\,\mathrm{km}\,\mathrm{s}^{-1}\,\mathrm{kpc}$ gives $h=1$ so no conversion of values necessary
- Assume all galaxies are blue in the UV bands where single Schechter functions are fit
- Parameters quoted as $M^*$, $\phi^*$, $\alpha$ for various $z$ values; linear slopes fit for using the original bin values
- Only take the UV bands due to abnormal red galaxy $\phi$ normalisation values in double Schechter U-band parameterisation