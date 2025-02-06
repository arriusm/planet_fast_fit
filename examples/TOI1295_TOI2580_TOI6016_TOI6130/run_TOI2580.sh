# 2024A+A...692A.220Ehrhardt__TOI1295_TOI2580_TOI6016_TOI6130.pdf
#
# TOI-2580:
#
# Table 8:
# Mass [Msun]     = 1.33 ± 0.08
# Radius [Rsun]   = 1.81 ± 0.06
#
# Table 9:
# P [d]           = 3.397750 ± 0.000002
# a [AU]          = 0.048 ± 0.003
# Mp [Mjup]       = 0.63 ± 0.08
# Rp [Rjup]       = 1.55 ± 0.05
# rho [g/cm3]     = 0.22 ± 0.04
# duration [min]  = 503±25                       <===== richtig:  272.2 min
#
# https://zenodo.org/records/13840492
# P [d]           =  3.3977506 +/- 0.0000020
# t_0 [d]         =  245 8839.4534 +/- 0.0005
# p               =  0.0867 +/- 0.0006           <===== besser:   0.090
# b               =  0.16 +/- 0.11
# e               =  0.08 +/- 0.04
# Ω [°]           =  114 +/- 30
# rho_* [kg/m^3]  =  110 +/- 30                  <====  richtig:  316.2 = 1.33*Msun/(4/3*np.pi*(1.81*Rsun)**3)
# K [m/s]         =  70 +/- 8
# M_{MaHPS}       =  31 +/- 8
# q_{1,TESS}      =  0.31 +/- 0.09     u1 = 2*sqrt(q1)*q2       = 0.1002197585309404
# q_{2,TESS}      =  0.09 +/- 0.08     u2 = sqrt(q1)*(1.-2*q2)  = 0.4565566777520619

# grep MaHPS toi2580.dat > TOI2580_RV_MaHPS.dat

# python  ../../planet_fast_fit.py  -save TOI2580_planet_fast_fit0.pdf  -tit "TOI 2580:"  -P 3.3977506  -t0 8839.4534  -rp 0.0867  -b 0.16  -e 0.08  -w 114  -rhostar 0.3162  -K 70.  -off 31.  -u1  0.10022   -u2 0.45655   -Ms 1.33   -zoom 3  -alpha 0.1  -rv TOI2580_RV_MaHPS.dat  TOI2580_TIC102713734_SEC59.fits

# python  ../../planet_fast_fit.py  -save TOI2580_planet_fast_fit.pdf  -tit "TOI 2580:"  -P 3.3977506   -t0 8839.4534  -rp 0.0867  -b 0.16  -e 0.08  -w 114  -rhostar 0.110  -K 70.  -off 31.  -u1  0.10022   -u2 0.45655   -Ms 1.33   -zoom 3  -alpha 0.1  -rv TOI2580_RV_MaHPS.dat  TOI2580_TIC102713734_SEC59.fits

# python  ../../planet_fast_fit.py  -save TOI2580_planet_fast_fit.pdf  -tit "TOI 2580:"  -P 6.7955012  -t0 8839.4534  -rp 0.0867  -b 0.16  -e 0.08  -w 114  -rhostar 0.3162  -K 70.  -off 31.  -u1  0.10022   -u2 0.45655   -Ms 1.33   -zoom 3  -alpha 0.1  -rv TOI2580_RV_MaHPS.dat  TOI2580_TIC102713734_SEC59.fits

python  ../../planet_fast_fit.py  -save TOI2580_planet_fast_fit.pdf  -tit "TOI 2580:"  -P 3.3977506  -t0 8839.4534  -rp 0.090  -b 0.16  -e 0.08  -w 114  -rhostar 0.3162  -K 70.  -off 31.  -u1  0.10022   -u2 0.45655   -Ms 1.33   -zoom 3  -alpha 0.1  -rv TOI2580_RV_MaHPS.dat  TOI2580_TIC102713734_SEC59.fits
