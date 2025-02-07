
# rhostar = 0.344*rhosun = 0.48497892482916016
# ap = 0.055 *AU / (1.62*Rsun) = 7.2974325
# rp_paper = sqrt(0.01197) = 0.1094
# rp_2     = 1.72*Rjup/(1.36*Rsun) = 0.1270
# t0 = 6416.40138-7000. +0.30 = -583.29862
# http://exoplanet.eu/catalog/wasp-94_a_b/


# start:
# python  ../../planet_fast_fit.py  -save TOI0107_SEC27_start.png   TOI0107_SEC27_TIC0092352620.fits

# python  ../../planet_fast_fit.py  -rhostar 0.8   -per 3.9501907  -t0 2047.236  -b 0.0  -zoom 5 TOI0107_TIC0092352620_SEC27.fits

# python  ../../planet_fast_fit.py  -rhostar 0.49  -per 3.9502  -t0 2047.236  -b 0.0  -rp 0.104 -u1 0.5   -zoom 10   -alpha 0.1  TOI0107_TIC0092352620_SEC27.fits

# python  ../../planet_fast_fit.py  -rhostar 0.484978 -P 3.9501907 -t0 -583.3288 -i 88.7 -rp 0.1094  -u1 0.5 -zoom 10  -alpha 0.05 -tit "A&A 572, A49 (2014)"  TOI0107_TIC0092352620_SEC27.fits

# my bestfit:
python  ../../planet_fast_fit.py  -rhostar 0.49  -P 3.9502  -t0 9047.236  -b 0.0  -rp 0.104  -u1 0.5  -zoom 3  -alpha 0.05  -tit "TOI 107:"  -save TOI0107_planet_fast_fit.pdf  TOI0107_TIC092352620_SEC27.fits

# TIC:
# python  ../../planet_fast_fit.py  -rhostar 0.32  -P 3.9502  -t0 9047.236  -b 0.5  -rp 0.109  -u1 0.5  -zoom 3  -alpha 0.05   -tit "TOI 107:"  -save TOI0107_planet_fast_fit.pdf   TOI0107_TIC0092352620_SEC27.fits
