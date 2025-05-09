# https://de.wikipedia.org/wiki/Helvetios_(Stern)
#   https://ui.adsabs.harvard.edu/abs/2019A&A...623A..72K
#     Masse 	(1,158 ± 0,058) M☉
#     Radius 	(1,183 ± 0,059) R☉
# https://en.wikipedia.org/wiki/51_Pegasi
#   https://ui.adsabs.harvard.edu/abs/2024ApJ...960L...6M
#     Mass	1.09 ± 0.02 M☉
#     Radius	1.152 ± 0.009 R☉
# https://de.wikipedia.org/wiki/Dimidium
#   # The Extrasolar Planets Encyclopaedia: Planet 51 Peg b. Abgerufen am 9. Mai 2015.
#   Orbitdaten:
#     Große Halbachse 	0,052 AE
#     Exzentrizität 	0,0069 +0,0069−0,0066
#     Umlaufdauer 	4,2308 ± 0,00004 d [2]
#   Weitere Daten
#     Radius 	        1,9 ± 0,3 RJ
#     Mindestmasse 	0,46 +0,01−0,06 MJ
#     Entfernung 	        14.7 pc
#     Methode 	        Radialgeschwindigkeitsmethode
#     Bahnneigung 	80  +19−10 deg
# https://en.wikipedia.org/wiki/51_Pegasi_b
#   Orbital characteristics
#     Aphelion	        0.0534 AU (7,990,000 km)
#     Perihelion	0.0520 AU (7,780,000 km)
#     Semi-major axis   0.0527 ± 0.0030 AU (7,880,000 ± 450,000 km)
#     Eccentricity	0.013 ± 0.012
#     Orbital period	4.230785 ± 0.000036 d
#   Physical characteristics
#     Mean radius     	1.2±0.1 RJ
#     Mass	        0.46±0.02 MJ
#     Temperature	1,250 K


#  0.46 * Mjup/Mearth = 146.20106718467042
#  1.2  * Rjup/Rearth = 13.450776877126417

# python  ../../planet_fast_fit.py  -save 51Peg_planet_fast_fit.pdf  -tit "51Peg:"   -P 4.230785  -e 0.013  -w 58  -i 90  -t0 9521.16303  -Mplan 146.20  -K 0  -off 11  -Mstar 1.09  -Rstar 1.152  -rhostar 0.   -Rplan 13.45   -zoom 3  -alpha 0.2  -rv 51Peg_RV_MaHPS.dat

#  Mplan = 146.20 and inc = 90  =>  K = 54.586

python  ../../planet_fast_fit.py  -save 51Peg_planet_fast_fit.pdf  -tit "51Peg:"   -P 4.230785  -e 0.013  -w 58  -i 90  -t0 9521.16303           -K 54.586  -off 11  -Mstar 1.09  -Rstar 1.152  -rhostar 0.   -Rplan 13.45   -zoom 3  -alpha 0.2  -rv 51Peg_RV_MaHPS.dat

# python  ../../planet_fast_fit.py  -save 51Peg_planet_fast_fit.pdf  -tit "51Peg:"   -P 4.230785  -e 0.013  -w 58  -i 80  -t0 9521.16303  -rp 0          -K 54.586  -off 11  -Mstar 1.09  -Rstar 1.152  -rhostar 0.   -Rplan 13.45   -zoom 3  -alpha 0.2  -rv 51Peg_RV_MaHPS.dat

# python  ../../planet_fast_fit.py  -save 51Peg_planet_fast_fit.pdf  -tit "51Peg:"   -P 4.230785  -e 0.013  -w 58  -i 45  -t0 9521.16303  -rp 0          -K 54.586  -off 11  -Mstar 1.09  -Rstar 1.152  -rhostar 0.   -Rplan 13.45   -zoom 3  -alpha 0.2  -rv 51Peg_RV_MaHPS.dat
