import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy.optimize import minimize
from astropy.modeling.models import Linear1D, Exponential1D
from astropy.cosmology import default_cosmology
import skypy.galaxies as skygal
import astropy.units as u
from speclite.filters import FilterResponse, load_filters


def gridcreate(name, y, x, ratio, z, **kwargs):
    plt.figure(name, figsize=(z*x, z*ratio*y))
    gs = gridspec.GridSpec(y, x, **kwargs)
    return gs


def combine_error(params, text):
    w = np.empty((2, len(params)), float)
    w[0, :] = params['{}p'.format(text)]
    w[1, :] = params['{}m'.format(text)]
    return w


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


def m_fit(p, x, y, o, factor1, factor2):
    c, m, a = p
    f = c/factor1 + a/factor2 * np.exp(-x * m)
    dfdc = 1/factor1
    dfdm = -x * a/factor2 * np.exp(-x * m)
    dfda = 1/factor2 * np.exp(-x * m)
    return np.sum((y - f)**2 / o**2), np.array([np.sum(-2 * (y - f) * d / o**2) for d in
                                                [dfdc, dfdm, dfda]])


def m_hess(p, x, y, o, factor1, factor2):
    c, m, a = p
    f = c/factor1 + a/factor2 * np.exp(-x * m)
    dfdc = 1/factor1
    dfdm = -x * a/factor2 * np.exp(-x * m)
    dfda = 1/factor2 * np.exp(-x * m)
    d2fdc2 = 0
    d2fdcdm = 0
    d2fdm2 = x**2 * a/factor2 * np.exp(-x * m)
    d2fdcda = 0
    d2fdmda = -x/factor2 * np.exp(-x * m)
    d2fda2 = 0

    hess = np.empty((3, 3), float)
    hess[0, 0] = np.sum(2 * dfdc**2 / o**2 - 2 * (y - f) * d2fdc2 / o**2)
    hess[0, 1] = hess[1, 0] = np.sum(2 * dfdc*dfdm / o**2 - 2 * (y - f) * d2fdcdm / o**2)
    hess[1, 1] = np.sum(2 * dfdm**2 / o**2 - 2 * (y - f) * d2fdm2 / o**2)

    hess[0, 2] = hess[2, 0] = np.sum(2 * dfdc*dfda / o**2 - 2 * (y - f) * d2fdcda / o**2)
    hess[1, 2] = hess[2, 1] = np.sum(2 * dfdm*dfda / o**2 - 2 * (y - f) * d2fdmda / o**2)
    hess[2, 2] = np.sum(2 * dfda**2 / o**2 - 2 * (y - f) * d2fda2 / o**2)
    return hess


def g_fit(p, x, y, o, factor1, factor2):
    c, m, a, mu = p
    f = c/factor1 + a/factor2 * np.exp(-0.5 * (x - mu)**2 * m)
    dfdc = 1/factor1
    dfdm = a/factor2 * np.exp(-0.5 * (x - mu)**2 * m) * (-0.5 * (x - mu)**2)
    dfda = 1/factor2 * np.exp(-0.5 * (x - mu)**2 * m)
    dfdu = a/factor2 * np.exp(-0.5 * (x - mu)**2 * m) * (x - mu) * m

    return np.sum((y - f)**2 / o**2), np.array([np.sum(-2 * (y - f) * d / o**2) for d in
                                                [dfdc, dfdm, dfda, dfdu]])


def g_hess(p, x, y, o, factor1, factor2):
    c, m, a, mu = p
    f = c/factor1 + a/factor2 * np.exp(-0.5 * (x - mu)**2 * m)
    dfdc = 1/factor1
    dfdm = a/factor2 * np.exp(-0.5 * (x - mu)**2 * m) * (-0.5 * (x - mu)**2)
    dfda = 1/factor2 * np.exp(-0.5 * (x - mu)**2 * m)
    dfdu = a/factor2 * np.exp(-0.5 * (x - mu)**2 * m) * (x - mu) * m

    d2fdc2 = 0
    d2fdcdm = 0
    d2fdm2 = a/factor2 * np.exp(-0.5 * (x - mu)**2 * m) * (-0.5 * (x - mu)**2)**2
    d2fdcda = 0
    d2fdmda = 1/factor2 * np.exp(-0.5 * (x - mu)**2 * m) * (-0.5 * (x - mu)**2)
    d2fda2 = 0
    d2fdcdu = 0
    d2fdmdu = a/factor2 * np.exp(-0.5 * (x - mu)**2 * m) * (x - mu) * (1 + m * -0.5 * (x - mu)**2)
    d2fduda = 1/factor2 * np.exp(-0.5 * (x - mu)**2 * m) * (x - mu) * m
    d2fdu2 = a/factor2 * np.exp(-0.5 * (x - mu)**2 * m) * m * ((x - mu)**2 * m - 1)

    hess = np.empty((4, 4), float)
    hess[0, 0] = np.sum(2 * dfdc**2 / o**2 - 2 * (y - f) * d2fdc2 / o**2)
    hess[0, 1] = hess[1, 0] = np.sum(2 * dfdc*dfdm / o**2 - 2 * (y - f) * d2fdcdm / o**2)
    hess[1, 1] = np.sum(2 * dfdm**2 / o**2 - 2 * (y - f) * d2fdm2 / o**2)

    hess[0, 2] = hess[2, 0] = np.sum(2 * dfdc*dfda / o**2 - 2 * (y - f) * d2fdcda / o**2)
    hess[1, 2] = hess[2, 1] = np.sum(2 * dfdm*dfda / o**2 - 2 * (y - f) * d2fdmda / o**2)
    hess[2, 2] = np.sum(2 * dfda**2 / o**2 - 2 * (y - f) * d2fda2 / o**2)

    hess[0, 3] = hess[3, 0] = np.sum(2 * dfdc*dfdu / o**2 - 2 * (y - f) * d2fdcdu / o**2)
    hess[1, 3] = hess[3, 1] = np.sum(2 * dfdm*dfdu / o**2 - 2 * (y - f) * d2fdmdu / o**2)
    hess[2, 3] = hess[3, 2] = np.sum(2 * dfdu*dfda / o**2 - 2 * (y - f) * d2fduda / o**2)
    hess[3, 3] = np.sum(2 * dfdu**2 / o**2 - 2 * (y - f) * d2fdu2 / o**2)

    return hess


def function_evaluation_lookup(cmau, ind1, ind2, x):
    c, m, a, u = cmau[ind1, ind2]
    if np.isnan(a) and np.isnan(u):
        return m * x + c
    elif np.isnan(u):
        return c + a * np.exp(-x * m)
    else:
        return c + a * np.exp(-0.5 * (x - u)**2 * m)


# Not currently used, but kept for completeness:
def galaxy_apparent_magnitudes(M_star, phi_star, magnitude_limit, sky_area, alpha, alpha0, alpha1,
                               weight, _filter, z_min, z_max, nz):
    cosmology = default_cosmology.get()
    z_range = np.linspace(z_min, z_max, nz)
    redshift, magnitude = skygal.schechter_lf(redshift=z_range, M_star=M_star, phi_star=phi_star,
                                              alpha=alpha, m_lim=magnitude_limit, sky_area=sky_area,
                                              cosmology=cosmology)
    spectral_coefficients = skygal.spectrum.dirichlet_coefficients(
        redshift=redshift, alpha0=alpha0, alpha1=alpha1, weight=weight)
    stellar_mass = skygal.spectrum.kcorrect.stellar_mass(coefficients=spectral_coefficients,
                                                         magnitudes=magnitude, filter=_filter)
    app_mags = skygal.spectrum.kcorrect.apparent_magnitudes(
        coefficients=spectral_coefficients, redshift=redshift, filters=_filter,
        cosmology=cosmology, stellar_mass=stellar_mass)

    return app_mags  # + 5 * np.log10(1 + redshift)


def galaxy_parameterisation_wrapper(wavs, m_array, c_array, mag_lim, sky_area, filters,
                                    z_min, z_max, nz, a_array, u_array,
                                    filter_set_name=None, filter_names=None, filter_wavs=None,
                                    filter_responses=None, filter_units=None):
    if filters == 'generate':
        if np.any([var is None for var in [
                filter_set_name, filter_names, filter_wavs, filter_responses, filter_units]]):
            raise UserWarning("Supply all of the parameters to generate a filter!")
        generate_speclite_filters(
            filter_set_name, filter_names, filter_wavs, filter_responses, filter_units)
        filters = ['{}-{}'.format(filter_set_name, q) for q in filter_names]
    else:
        # Convert a single string to an iterable object.
        filters = np.atleast_1d(filters)
    wavs = np.atleast_1d(wavs)
    alpha0_blue = [2.079, 3.524, 1.917, 1.992, 2.536]
    alpha1_blue = [2.265, 3.862, 1.921, 1.685, 2.480]
    weight_blue = [3.47e+09, 3.31e+06, 2.13e+09, 1.64e+10, 1.01e+09]
    alpha0_red = [2.461, 2.358, 2.568, 2.268, 2.402]
    alpha1_red = [2.410, 2.340, 2.200, 2.540, 2.464]
    weight_red = [3.84e+09, 1.57e+06, 3.91e+08, 4.66e+10, 3.03e+07]

    galaxies = []
    for wav, _filter in zip(wavs, filters):
        # wav should be in microns!
        log_wav = np.log10(wav)

        Mstar = function_evaluation_lookup(m_array, c_array, a_array, u_array, 0, 0, log_wav)
        alpha = function_evaluation_lookup(m_array, c_array, a_array, u_array, 2, 0, log_wav)
        phistar = function_evaluation_lookup(m_array, c_array, a_array, u_array, 1, 0, log_wav)
        P = function_evaluation_lookup(m_array, c_array, a_array, u_array, 3, 0, log_wav)
        Q = function_evaluation_lookup(m_array, c_array, a_array, u_array, 4, 0, P)
        # phi*(z) = phi* * 10**(0.4 P z) = phi* exp(0.4 * ln(10) P z)
        # and thus Exponential1D being of the form exp(x / tau), tau = 1 / (0.4 * ln(10) P)
        tau = 1 / (0.4 * np.log(10) * P)
        blue_gal_dict = {
            'M_star': Linear1D(-Q, Mstar), 'phi_star': Exponential1D(phistar, tau),
            'magnitude_limit': mag_lim, 'sky_area': sky_area * u.deg**2, 'alpha': alpha,
            'alpha0': alpha0_blue, 'alpha1': alpha1_blue, 'weight': weight_blue, '_filter': _filter}
        Mstar = function_evaluation_lookup(m_array, c_array, a_array, u_array, 0, 1, log_wav)
        alpha = function_evaluation_lookup(m_array, c_array, a_array, u_array, 2, 1, log_wav)
        phistar = function_evaluation_lookup(m_array, c_array, a_array, u_array, 1, 1, log_wav)
        P = function_evaluation_lookup(m_array, c_array, a_array, u_array, 3, 1, log_wav)
        Q = function_evaluation_lookup(m_array, c_array, a_array, u_array, 4, 1, P)
        tau = 1 / (0.4 * np.log(10) * P)
        red_gal_dict = {
            'M_star': Linear1D(-Q, Mstar), 'phi_star': Exponential1D(phistar, tau),
            'magnitude_limit': mag_lim, 'sky_area': sky_area * u.deg**2, 'alpha': alpha,
            'alpha0': alpha0_red, 'alpha1': alpha1_red, 'weight': weight_red, '_filter': _filter}

        blue_galaxies = galaxy_apparent_magnitudes(**blue_gal_dict, z_min=z_min, z_max=z_max, nz=nz)
        red_galaxies = galaxy_apparent_magnitudes(**red_gal_dict, z_min=z_min, z_max=z_max, nz=nz)
        galaxies.append(np.hstack([blue_galaxies, red_galaxies]))

    return galaxies


def generate_speclite_filters(group_name, filter_names, wavelength_list, response_list,
                              wavelength_unit):
    for filt_name, wavelength, response in zip(filter_names, wavelength_list, response_list):
        FilterResponse(wavelength=wavelength*wavelength_unit, response=response,
                       meta=dict(group_name=group_name, band_name=filt_name))

    # directory_name = '.'
    # f_name = whatever_filterresponse_gets_passed_to.save(directory_name)
    # f_name = os.path.join(directory_name, '{}-{}.ecsv'.format(group_name, filter_names[i]))
    # some_variable_name = speclite.filters.load_filter(f_name) -> this can just be passed through
    # galaxy_apparent_magnitude as a filters list, since speclite.filters.load_filters() takes
    # a list of .ecsv filepaths.


param_dtype = [('Citation', 'U100'), ('Ref', 'U100'), ('band', 'U5'), ('type', 'U1'),
               ('M', float), ('dMp', float), ('dMm', float), ('alpha', float), ('dalphap', float),
               ('dalpham', float), ('phi', float), ('dphip', float), ('dphim', float),
               ('Q', float), ('dQp', float), ('dQm', float), ('P', float), ('dPp', float),
               ('dPm', float), ('wavelength', float), ('h', float), ('z0', float),
               ('logphi', float), ('logphip', float), ('logphim', float), ('phi (unc)', float),
               ('dphip (unc)', float), ('dphim (unc)', float), ('M (unc)', float),
               ('dMp (unc)', float), ('dMm (unc)', float)]
params = np.genfromtxt('galaxy_count_parameters.csv', delimiter=',', dtype=param_dtype,
                       skip_header=1)


cmau_array = np.zeros((5, 2, 4), float) * np.nan

gs = gridcreate('12313', 3, 4, 0.8, 5)

for i, (varname, label, per_err, sys_err, w) in enumerate(zip(
    ['M', 'phi', 'alpha', 'P'], [r'$M^*_0$ / AB mag', r'$\phi^*_0$ / Mpc$^{{-3}}$ mag$^{{-1}}$',
                                 r'$\alpha$', 'P / $z^{{-1}}$'],
        [0.05, 0.05, 0.05, 0.05], [0, 0.001, 0.01, 0.1], [3, 4, 4, 3])):
    ax = plt.subplot(gs[i])
    for j, (typing, col) in enumerate(zip(['b', 'r'], ['b', 'r'])):
        q = (params['type'] == typing) & ~np.isnan(params[varname])
        x = np.log10(params['wavelength'][q] / 1000)
        y = params[varname][q]
        dy = np.sqrt(np.mean(combine_error(params[q], 'd{}'.format(varname)), axis=0)**2 +
                     np.abs(per_err * params[varname][q])**2 + sys_err**2)
        ax.errorbar(10**x, y, yerr=dy, fmt='{}.'.format(col), ls='None', alpha=0.6)

        if varname == 'P':
            _filter = y != 0
        elif varname == 'M':
            # Don't fit for 450um M* value since it has an unknown zero point, having
            # M_sun (in SCUBA-2 bandpass) as its defined M=0.
            _filter = 10**x < 400
        else:
            _filter = np.ones_like(y, bool)
        if i == 0 or (i == 1 and j == 0):
            run_flag = 1
        elif i == 1 and j == 1:
            run_flag = 2
        else:
            run_flag = 0
        if run_flag == 1:
            if i == 1 and j == 0:
                factor1, factor2 = 1000, 1e5
                x0 = [0.002 * factor1, 8, 0.004 * np.exp(-0.7*8) * factor2]
                args = (x[_filter], y[_filter], dy[_filter], factor1, factor2)
            else:
                factor1, factor2 = 1, 1
                x0 = [-22, 1.5, 5 * np.exp(-0.8 * 1.5)]
                args = (x[_filter], y[_filter], dy[_filter], factor1, factor2)
        elif run_flag == 2:
            factor1, factor2 = 1000, 1000
            x0 = [0.001 * factor1, 15, 0.004 * factor2, 0]
            args = (x[_filter], y[_filter], dy[_filter], factor1, factor2)
        else:
            args = (x[_filter], y[_filter], dy[_filter])
            x0 = [1, 0.1]
        f_fit = [lin_fit, m_fit, g_fit][run_flag]
        f_hess = [lin_hess, m_hess, g_hess][run_flag]
        res = minimize(f_fit, x0=x0, args=args, jac=True, hess=f_hess, method='newton-cg')
        sigs = 1 / np.sqrt(np.diag(f_hess(res.x, *args)))
        _x = np.linspace(np.amin(x), np.amax(x), 10000)
        if run_flag == 1:
            c, m, a = res.x
            c, a = c / factor1, a / factor2
            ax.plot(10**_x, c + a * np.exp(-_x * m), '{}-'.format(col),
                    label=r'a = {:.3e}$\pm${:.3e},' '\n' r'm = {:.{w}f}$\pm${:.{w}f},' '\n' r'c = {:.{w}f}$\pm${:.{w}f}'.format(
                        a, sigs[2]/factor2, m, sigs[1], c, sigs[0]/factor1, w=w))
        elif run_flag == 2:
            c, m, a, mu = res.x
            c, a = c / factor1, a / factor2
            ax.plot(10**_x, c + a * np.exp(-0.5 * (_x - mu)**2 * m), '{}-'.format(col),
                    label=r'u = {:.3e}$\pm${:.3e},' '\n' r'a = {:.3e}$\pm${:.3e},' '\n' r'm = {:.{w}f}$\pm${:.{w}f},' '\n' r'c = {:.{w}e}$\pm${:.{w}e}'.format(
                        mu, sigs[3], a, sigs[2]/factor2, m, sigs[1], c, sigs[0]/factor1, w=w))
        else:
            ax.plot(10**_x, res.x[0] + _x * res.x[1], '{}-'.format(col),
                    label=r'm = {:.{w}f}$\pm${:.{w}f},' '\n' r'c = {:.{w}f}$\pm${:.{w}f}'.format(
                        res.x[1], sigs[1], res.x[0], sigs[0], w=w))
        if run_flag == 1 or run_flag == 2:
            cmau_array[i, j, 0] = c
        else:
            cmau_array[i, j, 0] = res.x[0]
        cmau_array[i, j, 1] = res.x[1]
        if len(res.x) > 2:
            cmau_array[i, j, 2] = a
        if len(res.x) > 3:
            cmau_array[i, j, 3] = mu
    ax.set_xlabel(r'$\lambda$ / $\mu$m', fontsize=14)
    ax.set_ylabel(label, fontsize=14)
    ax.set_xscale('log')
    ax.legend(fontsize=10)

ax = plt.subplot(gs[4])
for typing, j in zip(['b', 'r'], [0, 1]):
    q = (params['type'] == typing) & ~np.isnan(params['P']) & ~np.isnan(params['Q'])
    x = params['P'][q]
    y = params['Q'][q]
    dx = np.sqrt(np.mean(combine_error(params[q], 'dP'), axis=0)**2 +
                 np.abs(0.05 * params['P'][q])**2 + 0.1**2)
    dy = np.sqrt(np.mean(combine_error(params[q], 'dQ'), axis=0)**2 +
                 np.abs(0.05 * params['Q'][q])**2 + 0.1**2)
    ax.errorbar(x, y, xerr=dx, yerr=dy, fmt='{}.'.format(typing), ls='None', zorder=0, alpha=0.6)

    rng = np.random.default_rng()
    N = 25000
    samples = np.empty((np.sum((y != 0)), N, 2), float)
    ii = 0
    for i in range(len(y)):
        if y[i] != 0:
            # Based on the 2-D ellipses of Loveday+2012, P-Q error circles are
            # tilted ~25degrees in the anti-correlation direction. Since we flip
            # P-Q to Q-P (i.e., Loveday+12 has Q on the x-axis, but we want Q on
            # the y-axis), we have to change the correlation angle to -(90-24)=-66.
            # Assume this angle unless that gives |rho| > 1, in which case reduce
            # to rho=0.9.
            deg = -66
            o22, o11 = max(dx[i], dy[i])**2, min(dx[i], dy[i])**2
            cos_t, sin_t = np.cos(np.radians(deg)), np.sin(np.radians(deg))
            o2p2 = (o11 - cos_t**2/sin_t**2 * o22) / (sin_t**2 - cos_t**4 / sin_t**2)
            o1p2 = (o11 - o2p2*sin_t**2) / cos_t**2
            corr = (o2p2 - o1p2) * cos_t * sin_t
            if np.abs(corr/dx[i]/dy[i]) >= 1:
                corr = np.sign(corr) * 0.9 * dx[i] * dy[i]
            samples[ii, :, [0, 1]] = rng.multivariate_normal(
                mean=[x[i], y[i]], cov=[[dx[i]**2, corr], [corr, dy[i]**2]], size=N).T
            ii += 1
    res = minimize(lin_fit, x0=[1, 5], args=(samples[:, :, 0].flatten(),
                   samples[:, :, 1].flatten(), 1), jac=True,
                   hess=lin_hess, method='newton-cg')
    cmau_array[4, j, 1] = res.x[1]
    cmau_array[4, j, 0] = res.x[0]
    sigs = 1 / np.sqrt(np.diag(lin_hess(res.x, samples[:, :, 0].flatten(),
                       samples[:, :, 1].flatten(), 1)))
    _x = np.linspace(np.amin(x), np.amax(x), 10000)
    xlims, ylims = ax.get_xlim(), ax.get_ylim()
    ax.plot(_x, res.x[0] + _x * res.x[1], '{}-'.format(typing),
            label=r'm = {:.4f}$\pm${:.4f}, c = {:.4f}$\pm${:.4f}'.format(
                res.x[1], sigs[1], res.x[0], sigs[0]))

    log_lam = np.linspace(np.log10(0.4), np.log10(450), 1000)
    _P = cmau_array[3, j, 1] * log_lam + cmau_array[3, j, 0]
    _Q = cmau_array[4, j, 1] * _P + cmau_array[4, j, 0]
    ax.plot(_P, _Q, c=typing, ls='-', lw=7)

    ax.set_xlim(*xlims)
    ax.set_ylim(*ylims)
ax.set_xlabel('P / $z^{{-1}}$', fontsize=14)
ax.set_ylabel('Q / mag $z^{{-1}}$', fontsize=14)
ax.legend(fontsize=8)


alpha0_blue = [2.079, 3.524, 1.917, 1.992, 2.536]
alpha1_blue = [2.265, 3.862, 1.921, 1.685, 2.480]
weight_blue = [3.47e+09, 3.31e+06, 2.13e+09, 1.64e+10, 1.01e+09]
alpha0_red = [2.461, 2.358, 2.568, 2.268, 2.402]
alpha1_red = [2.410, 2.340, 2.200, 2.540, 2.464]
weight_red = [3.84e+09, 1.57e+06, 3.91e+08, 4.66e+10, 3.03e+07]

for ind, wav, filters, file_name, offset, label, z_max in zip(
        [5, 6, 7, 8, 9, 10, 11], [0.435, 3.4, 2.2, 0.748, 0.1595, 3.6, 8, 3.4, 22.2],
        ['bessell-B', 'wise2010-w1', 'generate', 'sdss2010-i', 'generate', 'generate', 'generate',
         'wise2010-w1', 'wise2010-w4'],
        ['lin1999_b', 'jarrett_2017_w1', 'kochanek_2001_ks', 'metcalfe_2000_i',
         'gardnerbrownferguson_2000_fuv', 'ashby_2009_3.6', 'ashby_2009_8.0',
         'allwise_data_w1', 'allwise_data_w4'],
        [0, 2.699, 1.85, 0.43, 0, 2.785, 4.3916, 2.699, 6.6146],
        ['B', 'W1', 'K$_s$', 'i', 'FUV', '[3.6]', '[8.0]', 'W1', 'W4'],
        [1, 2, 5, 4, 3, 3, 3, 2, 2]):
    if filters == 'generate':
        if 'ks' in file_name:
            f_s_n = 'twomass'
            f_ns = ['k']
            data = np.loadtxt('filters/2MASS_2MASS.Ks.dat')
        elif '3.6' in file_name:
            f_s_n = 'spitzer'
            f_ns = ['I3_6']
            data = np.loadtxt('filters/Spitzer_IRAC.I1.dat')
        elif '8.0' in file_name:
            f_s_n = 'spitzer'
            f_ns = ['I8_0']
            data = np.loadtxt('filters/Spitzer_IRAC.I4.dat')
        else:
            f_s_n = 'hst'
            f_ns = ['fuv']
            data = np.loadtxt('filters/HST_STIS_FUV.F25QTZ.dat')
        filters_load = '{}-{}'.format(f_s_n, f_ns[0])
        f_ws = np.array([data[:, 0]]) / 10000
        f_rs = np.array([data[:, 1]])
        f_us = u.micron
        generate_speclite_filters(f_s_n, f_ns, f_ws, f_rs, f_us)
    else:
        filters_load = filters
    ax = plt.subplot(gs[ind])
    mag_lim, sky_area = 35, 30
    z_array_ = np.array([0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25,
                         1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5])
    z_array = z_array_[z_array_ <= z_max]
    col = np.array(['lightblue', 'b', 'aqua', 'g', 'y', 'goldenrod', 'r', 'm',
                    'violet', 'pink', 'lavenderblush', 'k'])

    app_bins = np.linspace(0, 30, 150)
    app_dens = np.zeros_like(app_bins)
    for i, (z_min, z_max) in enumerate(zip(z_array[:-1], z_array[1:])):
        nz = 11
        bins = np.linspace(-60, 50, 400)
        z_range = np.linspace(z_min, z_max, nz)
        cosmology = default_cosmology.get()
        dV_dz = (cosmology.differential_comoving_volume(z_range) * sky_area *
                 u.deg**2).to_value('Mpc3')
        dV = np.trapz(dV_dz, z_range)

        log_wav = np.log10(wav)
        Mstar = function_evaluation_lookup(cmau_array, 0, 0, log_wav)
        alpha1 = function_evaluation_lookup(cmau_array, 2, 0, log_wav)
        phistar = function_evaluation_lookup(cmau_array, 1, 0, log_wav)
        P = function_evaluation_lookup(cmau_array, 3, 0, log_wav)
        Q = function_evaluation_lookup(cmau_array, 4, 0, P)
        tau = 1 / (0.4 * np.log(10) * P)
        # Median-redshift Schechter function
        m_star1 = Linear1D(slope=-Q, intercept=Mstar)
        phi_star1 = Exponential1D(amplitude=phistar, tau=tau)
        L = 10 ** (0.4 * (m_star1(z_range) - bins[:, np.newaxis]))
        phi_model_z = 0.4 * np.log(10) * phi_star1(z_range) * L ** (alpha1+1) * np.exp(-L)
        phi_model1 = np.median(phi_model_z, axis=1) * dV / sky_area

        Mstar = function_evaluation_lookup(cmau_array, 0, 1, log_wav)
        alpha2 = function_evaluation_lookup(cmau_array, 2, 1, log_wav)
        phistar = function_evaluation_lookup(cmau_array, 1, 1, log_wav)
        P = function_evaluation_lookup(cmau_array, 3, 1, log_wav)
        Q = function_evaluation_lookup(cmau_array, 4, 1, P)
        tau = 1 / (0.4 * np.log(10) * P)
        m_star2 = Linear1D(slope=-Q, intercept=Mstar)
        phi_star2 = Exponential1D(amplitude=phistar, tau=tau)
        L = 10 ** (0.4 * (m_star2(z_range) - bins[:, np.newaxis]))
        phi_model_z = 0.4 * np.log(10) * phi_star2(z_range) * L ** (alpha2+1) * np.exp(-L)
        phi_model2 = np.median(phi_model_z, axis=1) * dV / sky_area

        w = skygal.spectrum.kcorrect.wavelength
        t = skygal.spectrum.kcorrect.templates
        rng = np.random.default_rng()

        redshift = rng.uniform(z_min, z_max, 100)
        spectral_coefficients = skygal.spectrum.dirichlet_coefficients(
            redshift=redshift, alpha0=alpha0_blue, alpha1=alpha1_blue, weight=weight_blue)

        kcorr = np.empty_like(redshift)
        for j in range(len(redshift)):
            z = redshift[j]
            f_list = load_filters(filters_load)
            f = f_list[0]
            fs = f.create_shifted(z)
            non_shift_ab_maggy, shift_ab_maggy = 0, 0
            for k in range(len(t)):
                non_shift_ab_maggy += spectral_coefficients[j, k] * f.get_ab_maggies(t[k], w)
                try:
                    shift_ab_maggy += spectral_coefficients[j, k] * fs.get_ab_maggies(t[k], w)
                except ValueError:
                    _t, _w = fs.pad_spectrum(t[k], w, method='edge')
                    shift_ab_maggy += spectral_coefficients[j, k] * fs.get_ab_maggies(_t, _w)
            # Backwards to Hogg+ astro-ph/0210394, our "shifted" bandpass is the rest-frame
            # as opposed to the observer frame.
            kcorr[j] = -2.5 * np.log10(1/(1+z) * shift_ab_maggy / non_shift_ab_maggy)
        # Loveday+2015 for absolute -> apparent magnitude conversion
        q = phi_model1 + phi_model2 > 0

        ax.plot(bins[q] + cosmology.distmod(np.mean(z_range)).value - offset +
                np.percentile(kcorr, 50), np.log10(phi_model1[q]+phi_model2[q]),
                c=col[i % len(col)], ls='-')
        app_dens += np.interp(app_bins, bins + cosmology.distmod(np.mean(z_range)).value -
                              offset + np.percentile(kcorr, 50), phi_model1+phi_model2)

    ax.plot(app_bins[app_dens > 0], np.log10(app_dens[app_dens > 0]), 'k-')

    if 'allwise' in file_name:
        f = np.loadtxt('data/wise_allwise.allwise_p3as_psd_19095.tbl')
        index = 2 if 'w1' in file_name else 4
        hist, bins = np.histogram(f[:, index], bins='auto')
        qqq = hist > 0
        # Polar cap >75 degrees galactic latitude:
        area = 2*np.pi*(np.sin(np.radians(90))-np.sin(np.radians(75))) * (180/np.pi)**2
        b, h = (bins[:-1]+np.diff(bins)/2)[qqq], np.log10(hist[qqq]/np.diff(bins)[qqq]/area)
        ax.plot(b, h, 'kx')
        ax.set_xlim(np.amin(b)-0.25, np.amax(b)+0.25)
        ax.set_ylim(np.amin(h)-0.1, np.amax(h)+0.1)
    else:
        f = np.loadtxt('data/{}_numdens.csv'.format(file_name), delimiter=',')
        if 'ashby' in file_name:
            f1 = np.loadtxt('data/{}_star_numdens.csv'.format(file_name), delimiter=',')
            f[:, 1] = np.log10(10**f[:, 1] - 10**f1[:, 1])
        ax.plot(f[:, 0], f[:, 1], ls='None', marker='x', c='k')
        ax.set_xlim(np.amin(f[:, 0])-0.25, np.amax(f[:, 0])+0.25)
        ax.set_ylim(np.amin(f[:, 1])-0.1, np.amax(f[:, 1])+0.1)

    ax.set_xlabel('{} / mag'.format(label), fontsize=14)
    ax.set_ylabel('log$_{10}$(sources / deg$^{-2}$ mag$^{-1}$)', fontsize=14)

plt.tight_layout()
plt.savefig('galaxy_count_parameters.pdf')

np.save('cmau_array.npy', cmau_array)
