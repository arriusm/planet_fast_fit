# 2024A+A...692A.220Ehrhardt__TOI1295_TOI2580_TOI6016_TOI6130.pdf
#
# TOI-1295:
#
# Table 8:
# Mass [Msun]     = 1.38 +/- 0.08
# Radius [Rsun]   = 1.70 +/- 0.03
#
# Table 9:
# P [d]           = 3.1968838 +/- 0.0000005
# a [AU]          = 0.047 +/- 0.002
# Mp [Mjup]       = 1.42 +/- 0.08
# Rp [Rjup]       = 1.40 +/- 0.08
# rho [g/cm3]     = 0.65 +/- 0.05
# duration [min]  = 370 +/- 15
#
# https://zenodo.org/records/13840492
# P [d]           = 3.1968840 +/- 0.0000005
# t_0 [d]         = 245 9913.37999 +/- 0.00020
# p               = 0.0840 +/- 0.0004
# b               = 0.555 +/- 0.020
# e               = 0.024 +/- 0.020
# Ω [°]           = 80 +/- 40
# rho_* [kg/m^3]  = 400 +/- 20
# K [m/s]         = 158 +/- 7
# M_{MaHPS}       = 53 +/- 5
# q_{1,TESS}      = 0.35 +/- 0.05   u1 = 2*sqrt(q1)*q2      = 0.15382
# q_{2,TESS}      = 0.13 +/- 0.07   u2 = sqrt(q1)*(1.-2*q2) = 0.43779

# grep MaHPS toi1295.dat > TOI1295_RV_MaHPS.dat

python ../../planet_fast_fit.py  -save TOI1295_planet_fast_fit.pdf  -tit "TOI 1295:" -P 3.1968840  -t0 9913.37999  -rp 0.0840  -b 0.555  -e 0.024  -w 80  -rhostar 0.400   -K 158.  -off 53.  -u1 0.15382  -u2 0.43779   -Ms 1.38   -zoom 3  -alpha 0.15  -rv TOI1295_RV_MaHPS.dat  TOI1295_TIC219852584_SEC59.fits
