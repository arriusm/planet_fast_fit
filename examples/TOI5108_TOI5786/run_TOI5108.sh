# https://www.aanda.org/articles/aa/full_html/2025/02/aa51676-24/aa51676-24.html

# Table 3. Stellar parameters of TOI 5108 and TOI 5786.
# 
# R∗ (R⊙)                                     1.29 +0.04 −0.04
# M∗ (M⊙)                                     1.10 +0.03 −0.03
# 
# Table 4. Priors and final parameter values from the juliet fit of TOI-5108.
# 
# Fitting parameters
# 
# P (d)                   U[6.7, 6.8]         6.753581  +- 0.000007
# T0 (BJDTDB − 2457000)   U[2569.1, 2569.8]   2569.4778 +- 0.0005
# Rp/R∗                   U[0.001, 1]         0.0472 +- 0.0007 −0.0007
# b                       U[0, 1]             0.871 +0.009 −0.010
# √e sinω                 fixed 0
# √e cosω                 fixed 0
# K (m s−1)               U[0, 40]            9.9 +1.4 −1.3
# ρ∗ (kg/m3)              N[720,90]           745 +80 −70
# 
# Instrumental parameters
# 
# q1,TESS                 N[0.28,0.1]         0.26 +0.06 −0.06         u1 = 2*sqrt(q1)*q2      = 0.40792156108742283 
# q2,TESS                 N[0.33,0.1]         0.40 +0.10 −0.1          u2 = sqrt(q1)*(1.-2*q2) = 0.10198039027185568
# μMaHPS (m s−1)          U[−100, 100]        3.2 +1.3 −1.3                
# μSOPHIE (m s−1)         U[−35000, −34000]   −34692.7 +1.3 −1.3         
# σMaHPS (m s−1)          L[0.001, 100]       11.5 +1.2 −1.1                
# σSOPHIE (m s−1)         L[0.001, 100]       6.2 +1.2 −1.0              
# 
# TOI-5108 b parameters
# 
# a/R∗                                        12.2 +0.4 −0.4
# a (au)                                      0.073 ± 0.004
# T_14 (h)                                    2.2 ± 0.1
# i (deg)                                     85.91 ± 0.14
# Rp (R⊕)                                     6.6 ± 0.1
# Mp (M⊕)                                     32 ± 5
# rhoplan (kg/m3)                             600 ± 90
# Teq (K)                                     1180 ± 40


# awk '($4=="SOPHIE" && $3<10.){print $1,$2+34691,$3,$4}' toi5108.dat > TOI5108_RV_SOPHIE.dat
# awk '($4=="MAHPS"  && $3<10.){print $1,$2      ,$3,$4}' toi5108.dat > TOI5108_RV_MaHPS.dat
# awk '($4=="MAHPS"  && $3<10.){print $1,$2      ,$3,$4}' toi5108.dat > TOI5108_RV_MaHPS.dat

# ../../planet_fast_fit.py -save TOI5108_TIC0350348197.pdf   -Mstar 1.10 -Rstar 1.29  -rhostar 0      -P 6.753581  -t0 9569.4778  -y0 0.993 -y1 1.003 -rp 0.0472  -u1 0.408 -u2 0.102 -b 0.871  -K 9.9  -off 0 -e 0.0  -z 10 -alpha 0.5   TOI5108_TIC0350348197_SEC45.fits  -rv TOI5108_RV_SOPHIE.dat

  ../../planet_fast_fit.py -save TOI0186_planet_fast_fit.pdf   -Mstar 1.10 -Rstar 1.29  -rhostar 0      -P 6.753581  -t0 9569.4778  -y0 0.993 -y1 1.003 -rp 0.0472  -u1 0.408 -u2 0.102 -b 0.871  -K 9.9  -off 0 -e 0.0  -z 10 -alpha 0.5   TOI5108_TIC0350348197_SEC45.fits  -rv TOI5108_RV_MaHPS.dat
# ../../planet_fast_fit.py -save TOI5108_TIC0350348197.pdf   -Mstar 1.10   -a 12.2    -rhostar 0      -P 6.753581  -t0 9569.4778  -y0 0.993 -y1 1.003 -rp 0.0472  -u1 0.408 -u2 0.102 -b 0.871  -K 9.9  -off 0 -e 0.0  -z 10 -alpha 0.5   TOI5108_TIC0350348197_SEC45.fits  -rv TOI5108_RV_MaHPS.dat
# ../../planet_fast_fit.py -save TOI5108_TIC0350348197.pdf   -Mstar 1.10              -rhostar 0.745  -P 6.753581  -t0 9569.4778  -y0 0.993 -y1 1.003 -rp 0.0472  -u1 0.408 -u2 0.102 -b 0.871  -K 9.9  -off 0 -e 0.0  -z 10 -alpha 0.5   TOI5108_TIC0350348197_SEC45.fits  -rv TOI5108_RV_MaHPS.dat

# vary rhostar <-> b and u1,u2 and rp
# ../../planet_fast_fit.py -save TOI5108_TIC0350348197.pdf   -Mstar 1.10              -rhostar 1.0    -P 6.753581  -t0 9569.4778  -y0 0.993 -y1 1.003 -rp 0.047   -u1 0.4   -u2 0.1   -b 0.85   -K 9.9  -off 0 -e 0.0  -z 10 -alpha 0.5   TOI5108_TIC0350348197_SEC45.fits  -rv TOI5108_RV_MaHPS.dat

