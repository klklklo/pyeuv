import numpy as np
import xarray as xr
import pyeuv._misc as _m


class Euv91:
    '''
    EUV91 model class
    '''
    def __init__(self):
        self.dataset = _m.get_euv91_dataset_full()
        self.bands_dataset, self.lines_dataset = _m.get_euv91_dataset()

        self.coeffs = np.vstack([self.dataset['a0'],
                                 self.dataset['a1'],
                                 self.dataset['a2'],
                                 self.dataset['a3'],
                                 self.dataset['a4'],
                                 self.dataset['a5'],
                                 self.dataset['a6'],
                                 self.dataset['a7'],
                                 self.dataset['a8'],
                                 self.dataset['a9']]).transpose()

        self.bands_coeffs = np.vstack([self.bands_dataset['a0'],
                                 self.bands_dataset['a1'],
                                 self.bands_dataset['a2'],
                                 self.bands_dataset['a3'],
                                 self.bands_dataset['a4'],
                                 self.bands_dataset['a5'],
                                 self.bands_dataset['a6'],
                                 self.bands_dataset['a7'],
                                 self.bands_dataset['a8'],
                                 self.bands_dataset['a9']]).transpose()

        self.lines_coeffs = np.vstack([self.lines_dataset['a0'],
                                 self.lines_dataset['a1'],
                                 self.lines_dataset['a2'],
                                 self.lines_dataset['a3'],
                                 self.lines_dataset['a4'],
                                 self.lines_dataset['a5'],
                                 self.lines_dataset['a6'],
                                 self.lines_dataset['a7'],
                                 self.lines_dataset['a8'],
                                 self.lines_dataset['a9']]).transpose()

    def _get_f(self, *proxies):
        l = proxies[0].size
        if not all(x.size == l for x in proxies):
            raise Exception(f'The number of proxy values does not match. '
                            f'lya contained {proxies[0].size} elements, '
                            f'hei contained {proxies[1].size} elements, '
                            f'f107 contained {proxies[2].size} elements, '
                            f'f107avg contained {proxies[3].size} elements.')

        return np.vstack([np.hstack([np.append(1., f), np.append(1., f)]) for f in zip(*proxies)],
                         dtype=np.float64)

    def _check_types(self, *proxies):
        if not (isinstance(proxies[0], (float, int, list, np.ndarray, type(None))) and
                isinstance(proxies[1], (float, int, list, np.ndarray, type(None))) and
                isinstance(proxies[2], (float, int, list, np.ndarray)) and
                isinstance(proxies[3], (float, int, list, np.ndarray))):

            raise TypeError(f'Only float, int, list and np.ndarray. '
                            f'lya was {type(proxies[0]).__name__}, '
                            f'hei was {type(proxies[1]).__name__}, '
                            f'f107 was {type(proxies[2]).__name__}, '
                            f'f107avg was {type(proxies[3]).__name__}')
        return True

    def prepare_data(self, lya, hei, f107, f107avg):
        if not (isinstance(lya, (float, int, list, np.ndarray, type(None))) and
                isinstance(hei, (float, int, list, np.ndarray, type(None))) and
                isinstance(f107, (float, int, list, np.ndarray)) and
                isinstance(f107avg, (float, int, list, np.ndarray))):

            raise TypeError(f'Only float, int, list and np.ndarray. '
                            f'lya was {type(lya).__name__}, '
                            f'hei was {type(hei).__name__}, '
                            f'f107 was {type(f107).__name__}, '
                            f'f107avg was {type(f107avg).__name__}')

        lya = np.array([lya]) if isinstance(lya, (type(None), int, float)) else np.array(lya)
        hei = np.array([hei]) if isinstance(hei, (type(None), int, float)) else np.array(hei)
        f107 = np.array([f107]) if isinstance(f107, (int, float)) else np.array(f107)
        f107avg = np.array([f107avg]) if isinstance(f107avg, (int, float)) else np.array(f107avg)

        if lya[0] is None and hei[0] is not None:
            lya = np.array([h * 3.7784687e9 + 8.4031723e10 for h in hei])

        elif hei[0] is None and lya[0] is not None:
            hei = np.array([(l - 8.4031723e10) / 3.7784687e9 for l in lya])

        elif lya[0] is None and hei[0] is None:
            raise Exception('lya and hei cannot be of type None at the same time.')

        return lya, hei, f107, f107avg

    def get_spectral_full(self, *, lya, hei, f107, f107avg):

        lya, hei, f107, f107avg = self. prepare_data(lya, hei, f107, f107avg)

        f = self._get_f(lya, hei, f107, f107avg)
        res = np.dot(self.coeffs, f.T)
        res_ener = res * 12400. * 1.602192e-12 / (self.dataset['center'].to_numpy().reshape(39, 1))/10
        return res, res_ener

    def get_spectral_bands(self, *, lya, hei, f107, f107avg):

        lya, hei, f107, f107avg = self. prepare_data(lya, hei, f107, f107avg)

        f = self._get_f(lya, hei, f107, f107avg)
        pflux = np.dot(self.bands_coeffs, f.T)
        eflux = pflux * 12400. * 1.602192e-12 / (self.bands_dataset['center'].to_numpy().reshape(23, 1))/10

        spectra = np.zeros((eflux.shape[1], eflux.shape[1], eflux.shape[1], eflux.shape[1], eflux.shape[0]))

        for i in range(eflux.shape[1]):
            spectra[i, i, i, i, :] = eflux[:, i]

        return xr.Dataset(data_vars={'euv_flux_spectra': (('Lya', 'HeI', 'F10.7', 'F10.7AVG', 'band_center'), spectra),
                                     'lband': ('band_number', self.bands_dataset['lband'].values),
                                     'uband': ('band_number', self.bands_dataset['uband'].values),
                                     'center': ('band_number', self.bands_dataset['center'].values)},
                          coords={'Lya': lya,
                                  'HeI': hei,
                                  'F10.7': f107,
                                  'F10.7AVG':  f107avg,
                                  'band_number': np.arange(23),
                                  'band_center': self.bands_dataset['center'].values})

    def get_spectral_lines(self, *, lya, hei, f107, f107avg):

        lya, hei, f107, f107avg = self.prepare_data(lya, hei, f107, f107avg)

        f = self._get_f(lya, hei, f107, f107avg)
        pflux = np.dot(self.lines_coeffs, f.T)
        eflux = pflux * 12400. * 1.602192e-12 / (self.lines_dataset['lambda'].to_numpy().reshape(16, 1))/10

        spectra = np.zeros((eflux.shape[1], eflux.shape[1], eflux.shape[1], eflux.shape[1], eflux.shape[0]))

        for i in range(eflux.shape[1]):
            spectra[i, i, i, i, :] = eflux[:, i]

        return xr.Dataset(data_vars={'euv_flux_spectra': (('Lya', 'HeI', 'F10.7', 'F10.7AVG', 'band_center'), spectra),
                       'wavelength': ('line_number', self.lines_dataset['lambda'].values)},
            coords={'Lya': lya,
                    'HeI': hei,
                    'F10.7': f107,
                    'F10.7AVG': f107avg,
                    'line_wavelength': self.lines_dataset['lambda'].values,
                    'line_number': np.arange(16)})

    def get_spectra(self, *, lya, hei, f107, f107avg):
        return self.get_spectral_bands(lya=lya, hei=hei, f107=f107, f107avg=f107avg), self.get_spectral_lines(f107=f107, lya=lya)