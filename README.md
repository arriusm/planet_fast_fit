usage: planet_fast_fit.py [-h] [-lc LC_FILENAME] [-rv RV_FILENAME] [-save SAVE] [-tex] 
                          [-tit TIT] [-Pmax PMAX] [-amax AMAX] [-t0 T0] [-P P] 
                          [-Mstar MSTAR] [-Mplan MPLAN] [-rp RP] [-Rstar RSTAR]
                          [-Rplan RPLAN] [-REB REB] [-a A] [-rhostar RHOSTAR] 
                          [-rhoplan RHOPLAN] [-b B] [-i I] [-u1 U1] [-u2 U2] 
                          [-e E] [-w W] [-norm NORM] [-x0 X0] [-x1 X1] [-y0 Y0] [-y1 Y1]
                          [-zoom ZOOM] [-alpha ALPHA] [-F0 F0] [-bg BG] [-K K] [-off OFFSET]
                          [files ...]

Vers. 2025-02-05 Arno Riffeser (arri@usm.lmu.de) 

python planet_fast_fit.py

positional arguments:
  files             files

options:
  -h, --help        show this help message and exit
  -lc LC_FILENAME   [] lc_filename
  -rv RV_FILENAME   [] rv_filename
  -save SAVE        [] save with filename (pdf/png)
  -tex              [False] tex Titel
  -tit TIT          [] Titel
  -Pmax PMAX        [0] Pmax
  -amax AMAX        [0] amax
  -t0 T0            [999] t0
  -P P              [2.] Periode P
  -Mstar MSTAR      [0] Masse Stern / Mstar
  -Mplan MPLAN      [0] Masse Planet / Mplan
  -rp RP            [0.1] Radiusverhaeltnis rp
  -Rstar RSTAR      [0] Radius Stern Rstar
  -Rplan RPLAN      [0] Radius Planet Rplan
  -REB REB          [0] Radius EB 2nd star
  -a A              [0.] Grosse Halbachse a
  -rhostar RHOSTAR  [1.0] mittl. Sterndichte rhostar
  -rhoplan RHOPLAN  [0.0] mittl. Planetendichte rhoplan
  -b B              [0.] Impaktparameter b
  -i I              [0.] Inklination i
  -u1 U1            [0.] Limb-Darkening u1
  -u2 U2            [0.] Limb-Darkening u2
  -e E              [0.] Excentrizitaet e
  -w W              [90.] Argument Periastron w
  -norm NORM        [0.] norm
  -x0 X0            [0] x0
  -x1 X1            [0] x1
  -y0 Y0            [0] y0
  -y1 Y1            [0] y1
  -zoom ZOOM        [1.0] zoom
  -alpha ALPHA      [0.05] alpha
  -F0 F0            [1.0] F0
  -bg BG            [0.0] bg
  -K K              [100.0] K
  -off OFFSET       [0.0] offset
