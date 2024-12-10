import xarray as xr
import pandas as pd
import numpy as np
from pyeuv91._pyeuv91 import Euv91

# print(xr.open_dataset('F:\MainProjects\InProgress\pyeuv91\src\pyeuv91\_coeffs\euv91_lines_dataset.nc'))
# print(xr.open_dataset('F:\MainProjects\InProgress\pyeuv91\src\pyeuv91\_coeffs\euv91_bands_dataset.nc'))
a = [1,2,3,4]
b = [5,6,7,8]
c = [9,1,2,3]
d = [4,5,6,7]

a1 = 1
b1 = 2
c1 = 3
d1 = 4

# mas = np.array([a,b])
# mas1 = np.array([1,2]).reshape(2,1)
# print(mas.shape, mas1.shape)
# print(mas)
# print(mas1)
# print(mas / mas1)

# пример из Тобишки
# date = 78269 day from
# lya = 5.16468e+11
# hei = 3.18297e+11
# f107 = 149
# f107avg = 150

# #78269
# lya = 5.16468e11
# hei = 3.18297e11
# f107 = 149
# f107avg = 150
#
# #78270
# lya = 5.19000e+11
# hei = 3.29632e+11
# f107 = 146
# f107avg = 151
#
# #78271
# lya = 5.20651e+11
# hei = 3.40968e+11
# f107 = 148
# f107avg = 152
#
# #78272
# lya = 5.22283e+11
# hei = 3.48525e+11
# f107 = 148
# f107avg = 152
#
# #78273
# lya = 5.25000e+11
# hei = 3.46377e+11
# f107 = 143
# f107avg = 153

# lya = [5.16468e11, 5.19000e+11, 5.20651e+11, 5.22283e+11, 5.25000e+11]
# hei = [3.18297e11, 3.29632e+11, 3.40968e+11, 3.48525e+11, 3.46377e+11]
# f107 = [149, 146, 148, 148, 143]
# f107avg = [150, 151, 152, 152, 153]

lya = [5.16468e11]
hei = [3.18297e11]
f107 = [149]
f107avg = [150]


e = Euv91()
res = e.get_spectral_bands(lya=lya, hei=hei, f107=f107, f107avg=f107avg)
print(res[0])

# xr.Dataset().from_dataframe(pd.read_csv('../src/pyeuv91/_coeffs/full_2_complete.csv')).to_netcdf('../src/pyeuv91/_coeffs/full_2_complete.nc')

# print(res['euv_flux_spectra'])
# print(e.get_spectral_bands(lya=a1, hei=b1, f107=c1, f107avg=d1))