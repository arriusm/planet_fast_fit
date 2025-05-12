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

# cd TOI1295_LC
# for i in *.fits; do python ../../../planet_fast_fit.py  -save ${i}.pdf  -tit "TOI 1295:" -P 3.1968840  -t0 9913.37999  -rp 0.0840  -b 0.555  -e 0.024  -w 80  -rhostar 0.400   -K 158.  -off 0.  -u1 0.15382  -u2 0.43779   -Ms 1.38   -zoom 3  -alpha 0.15 -y0 0.985 -y1 1.005  ${i}; done

python ../../planet_fast_fit.py  -save TOI1295_planet_fast_fit.pdf  -tit "TOI 1295:" -P 3.1968840  -t0 9913.37999  -rp 0.0840  -b 0.555  -e 0.024  -w 80  -rhostar 0.400   -K 158.  -off 53.  -u1 0.15382  -u2 0.43779   -Ms 1.38   -zoom 3  -alpha 0.15  -rv TOI1295_RV_MaHPS.dat  TOI1295_TIC0219852584_SEC59.fits


# AR_250511
#  output parameters:
#  x0 [d]           = 9911.79197493253
#  x1 [d]           = 9937.294080692209
#  y0               = 0.9810578376054764
#  y1               = 1.0176437944173813
#  y0_rv [km/s]     = -307.3
#  y1_rv [km/s]     = 529.2
#  t0 [d]           = 9926.17135
#  period [d]       = 3.199632
#  a/Rstar          = 6.380692182952806
#  a [AU]           = 0.03961224062084697
#  rp = Rplan/Rstar = 0.09095
#  ecc              = 0.0
#  w                = 90.0
#  u1               = 0.6
#  u2               = -0.3
#  q1 (juliet)      = 0.09
#  q2 (juliet)      = 1.0
#  inc [deg]        = 85.11819682188312
#  b                = 0.543
#  rhostar [g/cm3]  = 0.48
#  rhoplan [g/cm3]  = 0.7274145183768733
#  Mstar [Msun]     = 0.81
#  Mplan [Mjup]     = 0.9674199720488311
#  Mplan [Mearth]   = 307.47354854196334
#  Rstar [Rsun]     = 1.334950075084317
#  Rplan [Rjup]     = 1.1814960776048886
#  Rplan [Rearth]   = 13.243366767552828
#  bg               = 0.0
#  F0               = 1.0
#  K [m/s]          = 153.0
#  offset [m/s]     = 62.0
#  dur [h]          = 3.727860948402518
#  dur [min]        = 223.6716569041511
#  depth            = 8.85876492975035
#
# python ../../planet_fast_fit.py  -save TOI1295_planet_fast_fit_AR_250511.pdf  -tit "TOI 1295:" -P 3.199632  -t0 9926.17135  -rp 0.09095  -b 0.543  -e 0.  -w 90  -rhostar 0.48   -K 153.  -off 62.  -u1 0.6  -u2 -0.3   -Ms 0.81    -zoom 3  -alpha 0.15  -rv TOI1295_RV_MaHPS.dat  TOI1295_TIC0219852584_SEC59.fits

