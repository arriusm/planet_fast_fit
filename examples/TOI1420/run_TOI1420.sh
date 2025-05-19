# https://exofop.ipac.caltech.edu/tess/target.php?id=321857016

# TOI-1420:
#
# Radius [Rsun]   = 0.940037 +/- 0.0535176
# Mass [Msun]     = 0.935    +/- 0.117864
# Density [g/cm3] = 1.587068 +/- 0.379394

# Exoplanet Archive
# t_0 [d]         = 2459517.43305 
# P [d]           = 6.9561063 +/- 0.0000017
# inc [deg]       = 88.58 +/- 0.13
# b               = 0.412 +/- 0.036
# rp              = 0.11816 +/- 0.00059
# ap =            = 16.53 +/- 0.47
# Rp [Rearth]     = 11.89 +/- 0.33
# Mp [Mearth]     = 25.1 +/- 3.8
# a [AU]          = 0.071 +/- 0.0012
# e               = 0.17
# q_{1,TESS}      = ...  u1 = 2*sqrt(q1)*q2      = 
# q_{2,TESS}      = ...  u2 = sqrt(q1)*(1.-2*q2) = 
# 
# Depth [ppm]     = 13960 +/- 140
# Duration [h]    = 3.372 +/- 0.15

python ../../planet_fast_fit.py  -save TOI1420_planet_nofit.pdf  -tit "TOI 1420 (Exoplanet_Archive):"  -zoom 15  -alpha 0.15  TOI1420_TIC0321857016_SEC56.fits
 
# Exoplanet Archive
# python ../../planet_fast_fit.py  -save TOI1420_Exoplanet_Archive.pdf  -tit "TOI 1420 (Exoplanet_Archive):" -P 6.9561063  -t0 9517.43305   -rp 0.11816  -b 0.412   -Rs 0.940037  -Ms 0.935  -Mp 25.1  -off 0.  -u1 0.5  -u2 0.  -e 0.  -w 90  -zoom 15  -alpha 0.15  TOI1420_TIC0321857016_SEC56.fits
# 
python ../../planet_fast_fit.py  -save TOI1420_planet_fast_fit.pdf  -tit "TOI 1420 (Exoplanet_Archive):" -P 6.9561063  -t0 9517.43305   -rp 0.11816  -b 0.5   -Rs 0.940037  -Ms 0.935  -rhop 1.0  -off 0.  -u1 0.5  -u2 0.  -e 0.  -w 90  -zoom 15  -alpha 0.15  TOI1420_TIC0321857016_SEC56.fits
 
# WST
python ../../planet_fast_fit.py  -save TOI1420_i_WST_43cm.tab_nofit.pdf  -tit "TOI 1420:"  -zoom 15 -alpha 0.5  TOI1420_i_WST_43cm.tab

python ../../planet_fast_fit.py  -save TOI1420_i_WST_43cm.tab_planet_fast_fit.pdf  -tit "TOI 1420:"  -P 6.9561063  -t0 9517.43305   -rp 0.11816  -b 0.5   -Rs 0.940037  -Ms 0.935  -rhop 1.0  -off 0.  -u1 0.5  -u2 0.  -e 0.  -w 90  -zoom 15 -alpha 0.5  -norm 0.997  TOI1420_i_WST_43cm.tab
