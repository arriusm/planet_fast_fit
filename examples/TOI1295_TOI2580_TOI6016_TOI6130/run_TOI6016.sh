# 2024A+A...692A.220Ehrhardt__TOI1295_TOI2580_TOI6016_TOI6130.pdf
#
# TOI-6016:
#
# Table 8:
# Mass [Msun]     = 1.31 ± 0.08
# Radius [Rsun]   = 1.51 ± 0.03
#
# Table 9:
# P [d]           = 4.023687 ± 0.000003
# a [AU]          = 0.055 ± 0.002
# Mp [Mjup]       = 1.17 ± 0.09
# Rp [Rjup]       = 1.22 ± 0.03
# rho [g/cm3]     = 0.81 ± 0.08
# duration [min]  = 342 ± 13
#
# https://zenodo.org/records/13840492
# P [d]           =  4.0236869 +/- 0.0000030
# t_0 [d]         =  245 9877.79300 +/- 0.00030
# p               =  0.0820 +/- 0.0006
# b               =  0.44 +/- 0.04
# e               =  0.0 +- 0.0
# Ω [°]           =  90.0 +- 0.0
# rho_* [kg/m^3]  =  550 +/- 30
# K [m/s]         =  125 +/- 8
# M_{MaHPS}       =  -4 +/- 7
# q_{1,TESS}      =  0.38 +/- 0.07  u1 = 2*sqrt(q1)*q2      = 0.308220700
# q_{2,TESS}      =  0.25 +/- 0.08  u2 = sqrt(q1)*(1.-2*q2) = 0.308220700


# grep MaHPS toi6016.dat > TOI6016_RV_MaHPS.dat

python  ../../planet_fast_fit.py  -save TOI6016_planet_fast_fit.pdf  -tit "TOI 6016:" -P 4.0236869  -t0 9877.79300  -rp 0.0820  -b 0.44  -e 0.0  -w 90  -rhostar 0.550  -K 125  -off -4  -u1 0.30822   -u2 0.30822   -Ms 1.31   -zoom 3  -alpha 0.15  -rv TOI6016_RV_MaHPS.dat  TOI6016_TIC327369524_SEC58.fits

