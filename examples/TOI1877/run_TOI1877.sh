# https://exofop.ipac.caltech.edu/tess/target.php?id=233497719

# TOI-1877:
#
# Radius [Rsun]   = 1.00959 ± 0.0453467
# Mass [Msun]     = 1.13 ± 0.153259
# Density [g/cm3] = 1.548349 ± 0.341742

# t_0 [d]         = 2460534.067363 ± 0.00085197465
# P [d]           = 3.801591 ± 0.0000019
# Rp [Rearth]     = 14.4297 ± 0.657427
# rp              = 0.13
# Depth [ppm]     = 18950 ± 0.984908
# Duration [h]    = 2.752 ± 0.04

python ../../planet_fast_fit.py  -save TOI1877_nofit.pdf  -tit "TOI 1877: "  -zoom 10  -alpha 0.1  TOI1877_TIC0233497719_SEC49.fits

python ../../planet_fast_fit.py  -save TOI1877_planet_fast_fit.pdf  -tit "TOI 1877: " -P 3.801591  -t0 10534.06736   -b 0.4   -Rs 1.00959  -Ms 1.13 -Rp 14.4297 -rhop 1.0  -off 0.  -u1 0.5  -u2 0.  -e 0.  -w 90  -zoom 10  -alpha 0.1  TOI1877_TIC0233497719_SEC49.fits

python ../../planet_fast_fit.py  -save TOI1877_r_WST_43cm_nofit.pdf  -tit "TOI 1877: "  TOI1877_r_WST_43cm.tab

python ../../planet_fast_fit.py  -save TOI1877_r_WST_43cm_planet_fast_fit.pdf  -tit "TOI 1877: " -P 3.801591  -t0 10534.06736   -b 0.4   -Rs 1.00959  -Ms 1.13 -Rp 14.4297 -rhop 1.0  -off 0.  -u1 0.5  -u2 0.  -e 0.  -w 90  -zoom 10  -alpha 0.3   -prec 5001  -norm 0.999 TOI1877_r_WST_43cm.tab 

#python ../../planet_fast_fit.py  -save TOI1877_r_WST_43cm_planet_fast_fit.pdf  -tit "TOI 1877: " -P 3.801591  -t0 10534.06736   -b 0.4   -Rs 1.00959  -Ms 1.13 -Rp 14.4297 -rhop 1.0  -off 0.  -u1 0.5  -u2 0.  -e 0.  -w 90  -zoom 10  -alpha 0.3   -prec 100001 TOI1877_r_WST_43cm.tab  TOI1877_TIC0233497719_SEC49.fits




