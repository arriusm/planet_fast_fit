# https://www.aanda.org/articles/aa/full_html/2025/02/aa51676-24/aa51676-24.html

# Table 3. Stellar parameters of TOI 5108 and TOI 5786.
# 
# R∗ (R⊙)                                  1.36 +0.03 -0.05
# M∗ (M⊙)                                  1.23 +0.04 −0.03
# 
# Table 4. Priors and final parameter values from the juliet fit of TOI-5108.
# 
# Fitting parameters
# 
# P (d)                   U[12.7, 12.85]      12.779107 ± 0.000015
# T0 (BJDTDB − 2457000)   U[3140.5, 3140.7]   3140.6139 ± 0.0004
# Rp/R∗                   U[0.001, 1]         0.0576 +0.0009 -0.0008
# b                       U[0, 1]             0.33 +0.09 -0.14
# √e sinω                 fixed               0
# √e cosω                 fixed               0
# K (m s−1)               U[0, 40]            17.4 +2.1 -2.0
# ρ∗ (kg/m3)              N[690,80]           682 +90 -80
# 
# Instrumental parameters
# 
# q1,TESS                 N[0.32,0.1]         0.38 +0.07 −0.06    u1 = 2*sqrt(q1)*q2      = 0.3821936681840765
# q2,TESS                 N[0.36,0.1]         0.31 ± 0.09         u2 = sqrt(q1)*(1.-2*q2) = 0.2342477321128211
# μMaHPS (m s−1)          U[−100, 100]        1.1 ± 2.3
# μSOPHIE (m s−1)         U[−29000, −28000]   −28388 ± 2
# σMaHPS (m s−1)          L[0.001, 100]       14.12 ± 3 
# σSOPHIE (m s−1)         L[0.001, 100]       2.3 +4.0 −2.0
# 
# TOI-5786 b parameters
# 
# a/R∗                                        18.0 ± 0.8
# a (au)                                      0.114 ± 0.007	
# T_14 (h)                                    4.6 ± 0.3	
# i (deg)                                     88.9 ± 0.4	
# Rp (R⊕)                                     8.54 ± 0.13	
# Mp (M⊕)                                     73±9		
# rhoplan (kg/m3)                             640 ± 80	
# Teq (K)                                     1040 ± 40       


awk '($4=="SOPHIE" && $3<10.){print $1,$2+28388,$3,$4}' toi5786.dat > TOI5786_RV_SOPHIE.dat
# awk '($4=="MAHPS"  && $3<10.){print $1,$2      ,$3,$4}' toi5108.dat > TOI5786_RV_MaHPS.dat
awk '($4=="MAHPS"           ){print $1,$2      ,$3,$4}' toi5786.dat > TOI5786_RV_MaHPS.dat

#../../planet_fast_fit.py -save TOI5786_TIC040292751.pdf   -Mstar 1.23 -Rstar 1.36  -rhostar 0      -P 12.779107  -t0 10140.6139  -y0 0.993 -y1 1.003 -rp 0.0576  -u1 0.382 -u2 0.234 -b 0.33  -K 17.4  -off 0 -e 0.0  -F0 1.0004  -z 10 -alpha 0.5   TOI5786_TIC040292751_SEC40.fits  -rv TOI5786_RV_SOPHIE.dat

#../../planet_fast_fit.py -save TOI5786_TIC040292751.pdf   -Mstar 1.23 -Rstar 1.36  -rhostar 0      -P 12.779107  -t0 10140.6139  -y0 0.993 -y1 1.003 -rp 0.0576  -u1 0.382 -u2 0.234 -b 0.33  -K 17.4  -off 0 -e 0.0 -F0 1.0004   -z 10 -alpha 0.5   TOI5786_TIC040292751_SEC40.fits  -rv TOI5786_RV_MaHPS.dat
#../../planet_fast_fit.py -save TOI5786_TIC040292751.pdf   -Mstar 1.23   -a  18.0   -rhostar 0      -P 12.779107  -t0 10140.6139  -y0 0.993 -y1 1.003 -rp 0.0576  -u1 0.382 -u2 0.234 -b 0.33  -K 17.4  -off 0 -e 0.0  -F0 1.0004  -z 10 -alpha 0.5   TOI5786_TIC040292751_SEC40.fits  -rv TOI5786_RV_MaHPS.dat


## ../../planet_fast_fit.py -save TOI5786_TIC040292751.pdf   -Mstar 1.23              -rhostar 0.682  -P 12.779107  -t0 10140.6139  -y0 0.993 -y1 1.003 -rp 0.0576  -u1 0.382 -u2 0.234 -b 0.33  -K 17.4  -off 0 -e 0.0  -F0 1.0004 -z 10 -alpha 0.5   TOI5786_TIC040292751_SEC40.fits  -rv TOI5786_RV_MaHPS.dat


# >>> 8688.5 - 9446.5
# -758.0
# >>> 10323.4-8688.
# 1635.3999999999996
# >>> 10323.4-8688.5
# 1634.8999999999996
# >>> 10323.4-9446.5
# 876.8999999999996
# >>> (10323.4-8688.5)/2
# 817.4499999999998
# >>> (10323.4-8688.5)/3
# 544.9666666666666
# >>> (10323.4-8688.5)/4
# 408.7249999999999
# >>> (10323.4-8688.5)/5
# 326.9799999999999
# >>> (10323.4-8688.5)/6
# 272.4833333333333
# >>> (10323.4-8688.5)/7
# 233.5571428571428
# >>> (10323.4-8688.5)/8
# 204.36249999999995
# >>> (10323.4-8688.5)/9
# 181.6555555555555
# >>> 10527.8-10323.4
# 204.39999999999964


#for i in TOI5786_LC/*.fits
for i in TOI5786_LC/hlsp_tess-spoc_tess_phot_0000000040292751-s0014_tess_v1_lc.fits  TOI5786_LC/tess2024003055635-s0074-0000000040292751-0269-s_lc.fits  TOI5786_LC/tess2024196212429-s0081-0000000040292751-0276-s_lc.fits 
do
  #../../planet_fast_fit.py    -Mstar 1.23              -rhostar 0.682  -P 6.998405   -t0 8685.288    -y0 0.99 -y1 1.01  -rp 0.0258  -bg=-0.001    -u1 0.382 -u2 0.234 -b 0.38    -K 17.4  -off 0 -e 0.0  -F0 1.0004 -z 2 -alpha 0.03   $i   -rv TOI5786_RV_MaHPS.dat
  ../../planet_fast_fit.py    -Mstar 1.23              -rhostar 0.682  -P 204.347   -t0 10323.40   -y0 0.99 -y1 1.01  -rp 0.04  -bg=-0.001    -u1 0.382 -u2 0.234 -b 0.9    -K 17.4  -off 0 -e 0.0  -F0 1.0004 -z 50 -alpha 0.2   $i   -rv TOI5786_RV_MaHPS.dat
done




