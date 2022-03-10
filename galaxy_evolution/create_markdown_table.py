import numpy as np


def combine_error(params, text):
    w = np.empty((2, len(params)), float)
    w[0, :] = params['{}p'.format(text)]
    w[1, :] = params['{}m'.format(text)]
    return w


param_dtype = [('Citation', 'U100'), ('Ref', 'U100'), ('band', 'U5'), ('type', 'U1'),
               ('M', float), ('dMp', float), ('dMm', float), ('alpha', float), ('dalphap', float),
               ('dalpham', float), ('phi', float), ('dphip', float), ('dphim', float),
               ('Q', float), ('dQp', float), ('dQm', float), ('P', float), ('dPp', float),
               ('dPm', float), ('wavelength', float),
               ('h', float), ('z0', float), ('logphi', float), ('logphip', float),
               ('logphim', float), ('phi (unc)', float), ('dphip (unc)', float),
               ('dphim (unc)', float), ('M (unc)', float), ('dMp (unc)', float),
               ('dMm (unc)', float)]
a = np.genfromtxt('../../Postdocs/2019- Exeter/Work/Galaxy Evolution/'
                  'galaxy_count_parameters.csv', delimiter=',', dtype=param_dtype, skip_header=1)

f = open('schechter_markdown_table.md', 'w+')
f.write(r'| Citation | Band | Wavelength / nm | Type | $$M^*_0$$ / AB mag | $$\phi^*_0$$ / Mpc$$^{-3}$$mag$$^{-1}$$ | $$\alpha$$ | Q / mag $$z^{-1}$$ | P / $$z^{-1}$$|' '\n')
f.write(r'| --- | --- | --- | --- | --- | --- | --- | --- | --- |' '\n')
for i in range(len(a)):
    text = r'| {} | {} | {} | {} | '.format(a[i]['Citation'], a[i]['band'],
                                            a[i]['wavelength'], a[i]['type'])
    for x in ['M', 'phi', 'alpha', 'Q', 'P']:
        if np.isnan(a[i][x]):
            text += r' |'
        else:
            if x == 'phi':
                text += r' {:.3e}±{:.3e} |'.format(a[i][x],
                                                   np.mean(combine_error(a[i], 'd{}'.format(x))))
            else:
                text += r' {:.3f}±{:.3f} |'.format(a[i][x],
                                                   np.mean(combine_error(a[i], 'd{}'.format(x))))
    f.write(text + '\n')

f.close()

cmau_array = np.load('../../Postdocs/2019- Exeter/Work/Galaxy Evolution/cmau_array.npy')
c_array, m_array = cmau_array[:, :, 0], cmau_array[:, :, 1]
a_array, u_array = cmau_array[:, :, 2], cmau_array[:, :, 3]

big_table = np.empty(10, dtype=[('param', 'U10'), ('type', 'U1'), ('c', float),
                                ('m', float), ('a', float), ('u', float)])

f = open('parameter_markdown_table.md', 'w+')
f.write(r'| Parameter | Galaxy type | c | m | a | u' + '\n')
f.write(r'| --- | --- | --- | --- | --- | --- |' + '\n')
for i, param in enumerate([r'$$M^*_0$$', r'$$\phi^*_0$$', r'$$\alpha$$', r'$$P$$', r'$$Q$$']):
    for j, _type in enumerate(['b', 'r']):
        if np.isnan(a_array[i, j]) and np.isnan(u_array[i, j]):
            f.write(r'| {} | {} | {:.4f} | {:.3f} |  |  |'.format(
                param, _type, c_array[i, j], m_array[i, j]) + '\n')
            c, m, a, u = c_array[i, j], m_array[i, j], np.nan, np.nan
        elif np.isnan(u_array[i, j]):
            f.write(r'| {} | {} | {:.4f} | {:.3f} | {:.4f} |  |'.format(
                param, _type, c_array[i, j], m_array[i, j], a_array[i, j]) + '\n')
            c, m, a, u = c_array[i, j], m_array[i, j], a_array[i, j], np.nan
        else:
            f.write(r'| {} | {} | {:.4f} | {:.3f} | {:.4f} | {:.3f} |'.format(
                param, _type, c_array[i, j], m_array[i, j], a_array[i, j], u_array[i, j]) + '\n')
            c, m, a, u = c_array[i, j], m_array[i, j], a_array[i, j], u_array[i, j]
        big_table['param'][2*i+j] = param
        big_table['type'][2*i+j] = _type
        big_table['c'][2*i+j] = c
        big_table['m'][2*i+j] = m
        big_table['a'][2*i+j] = a
        big_table['u'][2*i+j] = u

np.savetxt('galaxy_count_parameter_table.csv', big_table, delimiter=',', fmt='%s, %s, %f, %f, %f, %f', header='Parameter, type, c, m, a, u')
