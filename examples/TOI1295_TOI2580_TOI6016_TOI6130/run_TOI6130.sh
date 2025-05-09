# 2024A+A...692A.220Ehrhardt__TOI1295_TOI2580_TOI6016_TOI6130.pdf
#
# TOI-6130:
#
# Table 8:
# Mass [Msun]     = 1.16 ± 0.07
# Radius [Rsun]   = 1.16 ± 0.03
#
# Table 9:
# P [d]           = 2.392679 ± 0.000002
# a [AU]          = 0.036 ± 0.002
# Mp [Mjup]       = 1.05 ± 0.06
# Rp [Rjup]       = 1.28 ± 0.03
# rho [g/cm3]     = 0.64 ± 0.06
# duration [min]  = 145±5
#
# https://zenodo.org/records/13840492
# P [d]           =  2.3926789 +/- 0.0000020
# t_0 [d]         =  245 9849.63919 +/- 0.00016
# p               =  0.1111+/-0.0005
# b               =  0.788+/-0.007
# e               =  0.036+/-0.018
# Ω [°]           =  42+/-30
# rho_* [kg/m^3]  =  (9.9+/-0.6)e+02
# K [m/s]         =  145+/-6
# M_{MaHPS}       =  18+/-4
# q_{1,TESS}      =  0.21+/-0.06  u1 = 2*sqrt(q1)*q2      = 0.11915
# q_{2,TESS}      =  0.27+/-0.09  u2 = sqrt(q1)*(1.-2*q2) = 0.33911

# grep MaHPS toi6130.dat > TOI6130_RV_MaHPS.dat

python  ../../planet_fast_fit.py  -save TOI6130_planet_fast_fit.pdf  -tit "TOI 6130:" -P 2.3926789  -t0 9849.63919  -rp 0.1111  -b 0.788  -e 0.036  -w 42  -rhostar 0.99   -K 145.  -off 18.  -u1 0.11915   -u2 0.33911   -Ms 1.16   -zoom 3  -alpha 0.1  -rv TOI6130_RV_MaHPS.dat  TOI6130_TIC210083929_SEC56.fits

