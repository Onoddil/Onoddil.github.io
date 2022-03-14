import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy.optimize import minimize
import sys


def gridcreate(name, y, x, ratio, z, **kwargs):
    plt.figure(name, figsize=(z*x, z*ratio*y))
    gs = gridspec.GridSpec(y, x, **kwargs)
    return gs


def flat_fit(p, y, o):
    return np.sum((y - p[0])**2 / o**2), np.array([np.sum(-2 * (y - p[0]) / o**2)])


def flat_hess(p, y, o):
    return np.array([[np.sum(2 / o**2)]])


def lin_fit(p, x, y, o):
    f = p[0] + x * p[1]
    dfdc = 1
    dfdm = x
    return np.sum((y - f)**2 / o**2), np.array([np.sum(-2 * (y - f) * d / o**2) for d in
                                                [dfdc, dfdm]])


def lin_hess(p, x, y, o):
    hess = np.empty((2, 2), float)
    f = p[0] + x * p[1]
    dfdc = 1
    dfdm = x
    d2fdc2 = 0
    d2fdcdm = 0
    d2fdm2 = 0
    hess[0, 0] = np.sum(2 * dfdc**2 / o**2 - 2 * (y - f) * d2fdc2 / o**2)
    hess[0, 1] = hess[1, 0] = np.sum(2 * dfdc*dfdm / o**2 - 2 * (y - f) * d2fdcdm / o**2)
    hess[1, 1] = np.sum(2 * dfdm**2 / o**2 - 2 * (y - f) * d2fdm2 / o**2)
    return hess


def best_fit_log2_schechter_function(phi_mod, lstar_mod, alpha_mod, sigma_mod, xmin, xmax):
    # Fit the modified Schechter function, phi* (L/L*)^(1-alpha) \times
    # exp(-1/(2sigma^2) log10(1 + L/L*)**2) with the non-modified Schechter
    # function phi* (L/L*)^alpha exp(-L/L*).
    pass


def schechter_logl(p, x, y):
    # Schechter function in log10L space is ln(10) phi* (L/L*)^(alpha+1) exp(-L/L*)
    log_phi, log_lstar, alpha = p
    # f = np.log10(np.log(10) * 10**log_phi * (10**(x-log_lstar))**(alpha+1) * np.exp(-10**(x-log_lstar)))
    log_f = np.log10(np.log(10)) + log_phi + (alpha+1) * (x - log_lstar) - 10**(x - log_lstar) / np.log(10)

    return np.sum((y - log_f)**2)


def best_fit_double_power_law_schechter_function(
        log_phi_knee, logl_knee, alpha_low_pl, alpha_high_pl, xmin, xmax):
    # Fit a double power law phi_knee (L/Lknee)^alpha
    # (alpha = alpha_high > logl_knee, else alpha_low),
    # phi2 = phi_knee (L/Lo)^alpha / (L/Lo)^alpha2, with a Schecter function,
    # .
    log10_l_div_l0 = np.linspace(xmin, xmax, 10000)
    log10_phi = np.empty_like(log10_l_div_l0)
    q = logl_knee <= log10_l_div_l0
    log10_phi[q] = log_phi_knee + alpha_low_pl * (log10_l_div_l0[q] - logl_knee)
    log10_phi[~q] = log_phi_knee + alpha_high_pl * (log10_l_div_l0[~q] - logl_knee)

    res = minimize(schechter_logl, x0=[0.1*log_phi_knee, -5+logl_knee, -0.5],
                   args=(log10_l_div_l0, log10_phi))
    log_phi, log_lstar, alpha = res.x
    x = log10_l_div_l0
    plt.plot(log10_l_div_l0,
             np.log10(np.log(10) * 10**log_phi * (10**x / 10**log_lstar)**(alpha+1) * np.exp(-10**x / 10**log_lstar)), 'r-')
    print(log_phi_knee, logl_knee, -0.5, res.x)

    plt.plot(log10_l_div_l0, log10_phi, 'k-')
    plt.show()


if sys.argv[1] == 'mortlock':
    # From Mortlock+2017, MNRAS, 465, 672, table 4.
    z = np.array([0.5, 1, 1.5, 2, 2.5, 3.25])
    log10_phi1 = np.array([-2.59, -2.65, -2.82, -2.88, -3.49, -3.91])
    log10_phi1_err = np.array([0.06, 0.06, 0.07, 0.07, 0.18, 0.36])
    m_kstar = np.array([-22.77, -23.16, -23.36, -23.17, -23.18, -23.92])
    m_kstar_err = np.array([0.16, 0.15, 0.13, 0.16, 0.3, 0.66])
    alpha1 = np.array([-0.28, -0.72, -1.11, -0.77, 0.26, -0.59])
    alpha1_err = np.array([0.29, 0.23, 0.16, 0.23, 0.69, 1.46])
    log10_phi2 = np.array([-2.95, -3.44, -4.67, -3.97, -3.17, -3.81])
    log10_phi2_err = np.array([0.1, 0.29, 1.38, 0.33, 0.15, 0.45])
    alpha2 = np.array([-1.44, -1.59, -2.11, -1.97, -1.53, -1.91])
    alpha2_err = np.array([0.04, 0.11, 0.58, 0.15, 0.08, 0.15])
    name = 'mortlock'
if sys.argv[1] == 'cirasuolo':
    # From Cirasuolo+2007, MNRAS, 380, 585, table 2.
    z = np.array([0.5, 0.875, 1.125, 1.375, 1.625, 2])
    log10_phi1 = np.log10(np.array([1.5, 1.3, 1.3, 1, 0.9, 0.25]) * 1e-3)
    # Uncertainty on log10(x) is ox/x/ln(10)
    log10_phi1_err = np.array([0.2, 0.2, 0.3, 0.4, 0.3, 0.3]) * 1e-3 / 10**log10_phi1 / np.log(10)
    m_kstar = np.array([-22.71, -23.05, -23.17, -23.4, -23.45, -23.75])
    m_kstar_err = np.array([0.11, 0.18, 0.21, 0.25, 0.15, 0.2])
    alpha1 = np.array([-1.19, -1.44, -1.38, -1.43, -1.35, -1.35])
    # Last two bins are fixed for alpha, so set to have very large errorbars
    alpha1_err = np.array([0.05, 0.15, 0.19, 0.24, 2, 2])
    log10_phi2 = np.log10(np.array([2.3, 2, 1, 1, 0.7, 0.2]) * 1e-3)
    log10_phi2_err = np.array([0.2, 0.3, 0.3, 0.4, 0.3, 0.4]) * 1e-3 / 10**log10_phi2 / np.log(10)
    # Table 2 gives second bin as -23.48, but that's a massive outlier and makes way more
    # sense as -22.48, so we just move it to in line with the trend... (could also just remove
    # it and get the same answer)
    m_kstar2 = np.array([-22.36, -22.48, -22.65, -22.65, -22.81, -23.04])
    m_kstar2_err = np.array([0.07, 0.13, 0.17, 0.24, 0.1, 0.21])
    alpha2 = np.array([-0.12, -0.09, -0.05, -0.14, -0.1, -0.1])
    alpha2_err = np.array([0.1, 0.2, 0.4, 0.5, 2, 2])
    name = 'cirasuolo'
if sys.argv[1] == 'arnouts':
    # From Arnouts+2007, A&A, 476, 137, table 1.
    # First bin has asymmetric errorbars for phi* and M*, but just take the largest here:
    z = np.array([0.3, 0.5, 0.7, 0.9, 1.1, 1.35, 1.75])
    log10_phi1 = np.log10(np.array([2.38, 2.22, 2.21, 2.82, 1.99, 2.08, 1.27]) * 1e-3)
    log10_phi1_err = np.array([1.6, 1, 0.9, 1, 0.6, 0.5, 0.3]) * 1e-3 / 10**log10_phi1 / np.log(10)
    m_kstar = np.array([-22.83, -22.82, -22.95, -23.06, -23.19, -23.19, -23.49])
    m_kstar_err = np.array([0.56, 0.3, 0.2, 0.2, 0.2, 0.18, 0.15])
    alpha1 = np.array([-1.3, -1.3, -1.3, -1.3, -1.3, -1.3, -1.3])
    alpha1_err = np.array([0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2])
    log10_phi2 = np.log10(np.array([1.78, 1.47, 1.29, 1.58, 0.9, 0.41, 0.05]) * 1e-3)
    log10_phi2_err = np.array([0.5, 0.17, 0.13, 0.14, 0.08, 0.03, 0.02]) * \
        1e-3 / 10**log10_phi2 / np.log(10)
    m_kstar2 = np.array([-22.91, -22.55, -22.73, -22.83, -22.86, -22.84, -23.27])
    m_kstar2_err = np.array([0.49, 0.19, 0.15, 0.14, 0.17, 0.16, 0.15])
    alpha2 = np.array([-0.6, -0.3, -0.3, -0.3, 0, 0.3, 0.6])
    alpha2_err = np.array([0.2, 0.2, 0.2, 0.2, 0.3, 0.3, 0.3])
    name = 'arnouts'
if sys.argv[1] == 'willmer':
    # Willmer+2006, ApJ 647, 853
    # Take the largest of asymmetric errorbars
    z1 = np.array([0.3, 0.5, 0.7, 0.9, 1.1, 1.3])
    m_kstar = np.array([-20.36, -20.72, -21.15, -21.21, -21.38, -21.86])
    m_kstar_err = np.array([0.13, 0.07, 0.07, 0.03, 0.05, 0.08])
    log10_phi1 = np.log10(np.array([31.78, 33.4, 24.67, 27.27, 20.84, 13.44]) * 1e-4)
    log10_phi1_err = np.array([2.15, 1.77, 1.58, 0.42, 1.58,
                               2.71]) * 1e-4 / 10**log10_phi1 / np.log(10)
    alpha1 = np.array([-1.3, -1.3, -1.3, -1.3, -1.3, -1.3])
    alpha1_err = np.array([0.5, 0.5, 0.5, 0.5, 0.5, 0.5])
    z2 = np.array([0.3, 0.5, 0.7, 0.9, 1.1])
    log10_phi2 = np.log10(np.array([17.06, 14.15, 13.66, 10.72, 5.24]) * 1e-4)
    log10_phi2_err = np.array([1.65, 0.7, 1, 0.38, 0.95]) * \
        1e-4 / 10**log10_phi2 / np.log(10)
    m_kstar2 = np.array([-21.02, -20.97, -21.19, -21.11, -21.44])
    m_kstar2_err = np.array([0.18, 0.14, 0.06, 0.05, 0.08])
    alpha2 = np.array([-0.5, -0.5, -0.5, -0.5, -0.5])
    alpha2_err = np.array([0.5, 0.5, 0.5, 0.5, 0.5])
    name = 'willmer'
if sys.argv[1] == 'lim':
    # From Lim+2020, ApJ 889 80, table 9, taking largest of asymmetric errorbars
    z = np.array([0.6, 1.9, 3.75])
    alpha1 = np.array([-0.64, -0.39, 0.86])
    alpha1_err = np.array([0.44, 0.54, 1.16])
    m_kstar = -2.5 * np.array([11.91, 12.26, 12.36])
    m_kstar_err = -2.5 * np.array([0.32, 0.27, 0.31])
    # Correct for units being log(phi*) in units if Mpc^-3 dex^-1, with M*-M0=-2.5log(L*/L0)
    log10_phi1 = np.array([-3.24, -3.07, -3.17]) + np.log10(2.5)
    log10_phi1_err = np.array([0.36, 0.36, 0.38])
    name = 'lim'
if sys.argv[1] == 'beare':
    # From Beare+2019, ApJ 873 78, table 4.
    z = np.array([0.3, 0.5, 0.7, 0.9, 1.1])[1:]
    log10_phi1 = np.log10(np.array([2.34, 1.52, 1.54, 1.43, 0.99]) * 1e-3)[1:]
    log10_phi1_err = np.array([0.06, 0.16, 0.02, 0.11, 0.14])[1:] * 1e-3 / 10**log10_phi1 / np.log(10)
    m_kstar = np.array([-22.51, -22.88, -22.99, -23.11, -23.24])[1:]
    m_kstar_err = np.array([0.05, 0.09, 0.05, 0.06, 0.09])[1:]
    alpha1 = np.array([-1.2, -1.2, -1.2, -1.2, -1.2])[1:]
    alpha1_err = np.array([0.2, 0.2, 0.2, 0.2, 0.2])[1:]
    log10_phi2 = np.log10(np.array([3.14, 2.29, 2.25, 2.01, 1.75]) * 1e-3)[1:]
    log10_phi2_err = np.array([0.1, 0.1, 0.1, 0.2, 0.13])[1:] * \
        1e-3 / 10**log10_phi2 / np.log(10)
    m_kstar2 = np.array([-22.59, -22.88, -22.93, -22.97, -23])[1:]
    m_kstar2_err = np.array([0.04, 0.08, 0.05, 0.05, 0.05])[1:]
    alpha2 = np.array([-0.5, -0.5, -0.5, -0.5, -0.5])[1:]
    alpha2_err = np.array([0.2, 0.2, 0.2, 0.2, 0.2])[1:]
    name = 'beare'
if sys.argv[1] == 'zucca':
    # From Zucca+2009 A&A 508 1217, table 1 (early/spiral), take largest of asymmetric errobars
    z = np.array([0.225, 0.45, 0.65, 0.875])
    log10_phi1 = np.log10(np.array([3.59, 3.68, 2.92, 1.58]) * 1e-3)
    log10_phi1_err = np.array([0.53, 0.55, 0.62, 0.74]) * 1e-3 / 10**log10_phi1 / np.log(10)
    m_kstar = np.array([-20.46, -20.51, -21.08, -21.59])
    m_kstar_err = np.array([0.15, 0.14, 0.18, 0.29])
    alpha1 = np.array([-1.14, -0.97, -1.21, -1.78])
    alpha1_err = np.array([0.06, 0.13, 0.17, 0.27])
    log10_phi2 = np.log10(np.array([2.81, 2.54, 2.07, 1.16]) * 1e-3)
    log10_phi2_err = np.array([0.34, 0.22, 0.08, 0.37]) * \
        1e-3 / 10**log10_phi2 / np.log(10)
    m_kstar2 = np.array([-20.9, -20.63, -20.81, -21.69])
    m_kstar2_err = np.array([0.17, 0.14, 0.13, 0.27])
    alpha2 = np.array([-0.72, -0.38, -0.02, -1.29])
    alpha2_err = np.array([0.07, 0.14, 0.19, 0.29])
    name = 'zucca'
if sys.argv[1] == 'moutard_fuv':
    # Moutard+2020, MNRAS, 494, 1894, table 2 (for all three bands)
    z = np.array([0.175, 0.375, 0.525, 0.75, 1.1, 1.55, 2.15, 3])
    log10_phi1 = np.log10(np.array([4.85, 5.22, 4.13, 4.4, 4.97, 3.2, 2.82, 1.69]) * 1e-3)
    log10_phi1_err = np.array([0.35, 0.25, 0.5, 0.3, 0.57, 0.38, 0.11, 0.1]) * 1e-3 / \
        10**log10_phi1 / np.log(10)
    m_kstar = np.array([-18.269, -18.572, -18.797, -19.113, -19.554, -20.016, -20.261, -20.841])
    m_kstar_err = np.array([0.054, 0.038, 0.073, 0.041, 0.065, 0.074, 0.042, 0.046])
    alpha1 = np.array([-1.405, -1.369, -1.408, -1.402, -1.432, -1.446, -1.43, -1.43])
    alpha1_err = np.array([0.019, 0.017, 0.053, 0.038, 0.068, 0.074, 0.2, 0.2])
    name = 'moutard_fuv'
if sys.argv[1] == 'moutard_nuv':
    z = np.array([0.175, 0.375, 0.525, 0.75, 1.1, 1.55, 2.15, 3])
    log10_phi1 = np.log10(np.array([5.08, 6.01, 4.56, 4.24, 4.7, 3.13, 2.72, 1.71]) * 1e-3)
    log10_phi1_err = np.array([0.2, 0.26, 0.47, 0.4, 0.43, 0.23, 0.08, 0.11]) * 1e-3 / \
        10**log10_phi1 / np.log(10)
    m_kstar = np.array([-18.514, -18.798, -19.026, -19.416, -19.859, -20.367, -20.622, -21.152])
    m_kstar_err = np.array([0.025, 0.03, 0.062, 0.053, 0.052, 0.045, 0.034, 0.038])
    alpha1 = np.array([-1.399, -1.308, -1.364, -1.396, -1.385, -1.391, -1.4, -1.4])
    alpha1_err = np.array([0.011, 0.014, 0.043, 0.044, 0.046, 0.038, 0.2, 0.2])
    name = 'moutard_nuv'
if sys.argv[1] == 'moutard_u':
    z = np.array([0.175, 0.375, 0.525, 0.75, 1.1])
    log10_phi1 = np.log10(np.array([3.4, 2.41, 1.96, 1.46, 0.51]) * 1e-3)
    log10_phi1_err = np.array([0.36, 0.47, 0.08, 0.11, 0.63]) * 1e-3 / 10**log10_phi1 / np.log(10)
    m_kstar = np.array([-18.961, -19.5, -19.744, -20.244, -20.819])
    m_kstar_err = np.array([0.05, 0.04, 0.027, 0.023, 0.033])
    alpha1 = np.array([-1.568, -1.557, -1.56, -1.56, -1.56])
    alpha1_err = np.array([0.024, 0.051, 0.2, 0.2, 0.2])
    log10_phi2 = np.log10(np.array([7.16, 7.91, 5.59, 5.04, 4.28]) * 1e-3)
    log10_phi2_err = np.array([0.25, 0.3, 0.12, 0.09, 0.47]) * \
        1e-3 / 10**log10_phi2 / np.log(10)
    alpha2 = np.array([-0.213, -0.419, -0.506, -0.773, -1.187])
    alpha2_err = np.array([0.099, 0.092, 0.043, 0.042, 0.092])
    name = 'moutard_u'
# Missing data we could use from Ilbert et al. (2005, A&A, 439, 863), but not fitting for those
# bandpasses at the moment.

# best_fit_double_power_law_schechter_function(-2.52, 10.48, -2.2, -0.6, 10, 13)
# sys.exit()

gs = gridcreate('34342', 2, 3, 0.8, 6)
ax = plt.subplot(gs[0])
if name == 'willmer':
    z = z1
ax.errorbar(z, log10_phi1, yerr=log10_phi1_err, fmt='k.')
res = minimize(lin_fit, x0=[1, 0], args=(z, log10_phi1, log10_phi1_err), jac=True, hess=lin_hess,
               method='newton-cg')
sigs = 1 / np.sqrt(np.diag(lin_hess(res.x, z, log10_phi1, log10_phi1_err)))
ax.plot(z, res.x[0] + z * res.x[1], 'r-',
        label=r'0.4P$_1$ = {:.3f}$\pm${:.3f}, log$_{{10}}$($\phi*_1$) = {:.3f}$\pm${:.3f}'.format(
            res.x[1], sigs[1], res.x[0], sigs[0]))
ax.set_xlabel('z')
ax.set_ylabel(r'log${_10}$($\phi*_1$)')
ax.legend()
print(r'phi*1: 0.4P1 = {:.3f}$\pm${:.3f}, log10(phi*1) = {:.3f}$\pm${:.3f}'.format(
    res.x[1], sigs[1], res.x[0], sigs[0]))
print(r'P1 = {:.3f}$\pm${:.3f}, phi*1 = {:.7f}$\pm${:.7f}'.format(
    res.x[1] / 0.4, sigs[1]/0.4, 10**res.x[0], 10**res.x[0] * np.log(10) * sigs[0]))

ax = plt.subplot(gs[1])
ax.errorbar(z, m_kstar, yerr=m_kstar_err, fmt='k.')
res = minimize(lin_fit, x0=[1, 0], args=(z, m_kstar, m_kstar_err), jac=True, hess=lin_hess,
               method='newton-cg')
sigs = 1 / np.sqrt(np.diag(lin_hess(res.x, z, m_kstar, m_kstar_err)))
# Linear fit is y = m x + c but we have c - mx in M* = M*(0) - Qz
ax.plot(z, res.x[0] + z * res.x[1], 'r-',
        label=r'Q = {:.3f}$\pm${:.3f}, M* = {:.3f}$\pm${:.3f}'.format(-res.x[1], sigs[1],
                                                                      res.x[0], sigs[0]))
ax.set_xlabel('z')
ax.invert_yaxis()
ax.set_ylabel(r'M*')
ax.legend()
print(r'M*: Q = {:.3f}$\pm${:.3f}, M* = {:.3f}$\pm${:.3f}'.format(-res.x[1], sigs[1],
                                                                  res.x[0], sigs[0]))

ax = plt.subplot(gs[2])
if name == 'lim':
    y, dy, x = alpha1[:-1], alpha1_err[:-1], z[:-1]
else:
    y, dy, x = alpha1, alpha1_err, z
ax.errorbar(x, y, yerr=dy, fmt='k.')
res = minimize(flat_fit, x0=[0], args=(y, dy), jac=True, hess=flat_hess,
               method='newton-cg')
sigs = 1 / np.sqrt(np.diag(flat_hess(res.x, y, dy)))
ax.plot(x, res.x[0] * np.ones_like(x), 'r-',
        label=r'$\alpha_1$ = {:.3f}$\pm${:.3f}'.format(res.x[0], sigs[0]))
ax.set_xlabel('z')
ax.set_ylabel(r'$\alpha_1$')
ax.legend()
print(r'alpha1 = {:.3f}$\pm${:.3f}'.format(res.x[0], sigs[0]))

if name != 'lim' and name != 'moutard_fuv' and name != 'moutard_nuv':
    if name == 'willmer':
        z = z2
    ax = plt.subplot(gs[3])
    ax.errorbar(z, log10_phi2, yerr=log10_phi2_err, fmt='k.')
    res = minimize(lin_fit, x0=[1, 0], args=(z, log10_phi2, log10_phi2_err), jac=True,
                   hess=lin_hess, method='newton-cg')
    sigs = 1 / np.sqrt(np.diag(lin_hess(res.x, z, log10_phi2, log10_phi2_err)))
    ax.plot(z, res.x[0] + z * res.x[1], 'r-',
            label=r'0.4P$_2$ = {:.3f}$\pm${:.3f}, log$_{{10}}$($\phi*_2$) = {:.3f}$\pm${:.3f}'.
            format(res.x[1], sigs[1], res.x[0], sigs[0]))
    ax.set_xlabel('z')
    ax.set_ylabel(r'log${_10}$($\phi*_2$)')
    ax.legend()
    print(r'phi*2: 0.4P2 = {:.3f}$\pm${:.3f}, log10(phi*2) = {:.3f}$\pm${:.3f}'.format(
        res.x[1], sigs[1], res.x[0], sigs[0]))
    print(r'P2 = {:.3f}$\pm${:.3f}, phi*2 = {:.7f}$\pm${:.7f}'.format(
        res.x[1] / 0.4, sigs[1]/0.4, 10**res.x[0], 10**res.x[0] * np.log(10) * sigs[0]))

    if name == 'cirasuolo' or name == 'arnouts' or name == 'willmer' or name == 'beare' \
            or name == 'zucca':
        ax = plt.subplot(gs[4])
        ax.errorbar(z, m_kstar2, yerr=m_kstar2_err, fmt='k.')
        res = minimize(lin_fit, x0=[1, 0], args=(z, m_kstar2, m_kstar2_err), jac=True,
                       hess=lin_hess, method='newton-cg')
        sigs = 1 / np.sqrt(np.diag(lin_hess(res.x, z, m_kstar2, m_kstar2_err)))
        ax.plot(z, res.x[0] + z * res.x[1], 'r-',
                label=r'Q$_2$ = {:.3f}$\pm${:.3f}, M*$_2$ = {:.3f}$\pm${:.3f}'.format(
                    -res.x[1], sigs[1], res.x[0], sigs[0]))
        ax.set_xlabel('z')
        ax.invert_yaxis()
        ax.set_ylabel(r'M*')
        ax.legend()
        print(r'M*2: Q2 = {:.3f}$\pm${:.3f}, M*2 = {:.3f}$\pm${:.3f}'.format(-res.x[1], sigs[1],
                                                                             res.x[0], sigs[0]))

    ax = plt.subplot(gs[5])
    ax.errorbar(z, alpha2, yerr=alpha2_err, fmt='k.')
    res = minimize(flat_fit, x0=[0], args=(alpha2, alpha2_err), jac=True, hess=flat_hess,
                   method='newton-cg')
    sigs = 1 / np.sqrt(np.diag(flat_hess(res.x, alpha2, alpha2_err)))
    ax.plot(z, res.x[0] * np.ones_like(z), 'r-',
            label=r'$\alpha_2$= {:.3f}$\pm${:.3f}'.format(res.x[0], sigs[0]))
    ax.set_xlabel('z')
    ax.set_ylabel(r'$\alpha_2$')
    ax.legend()
    print(r'alpha2= {:.3f}$\pm${:.3f}'.format(res.x[0], sigs[0]))

plt.tight_layout()
plt.savefig('Fit Evolution/{}_evolution_parameters.pdf'.format(name))
