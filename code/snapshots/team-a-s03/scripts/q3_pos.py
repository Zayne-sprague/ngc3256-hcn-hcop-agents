# positions (offsets in arcsec from RA0,DE0) incl. Michiyama+2018 Table 1 coordinates
from q3_lib import *
def off(ra_s,de_s):   # sexagesimal h m s / d m s -> offsets from (RA0,DE0) in arcsec
    h,m,s=ra_s; d,mm,ss=de_s; ra=15*(h+m/60+s/3600); de=-(d+mm/60+ss/3600)
    return ((ra-RA0)*np.cos(np.deg2rad(DE0))*3600,(de-DE0)*3600)
# Michiyama+2018 Table 1 (J2000)
POS={'N':(0,2.55),'S':(-0.15,-2.55),'Bridge':(0.3,0.15),'Soutflow':(-0.21,-4.48),
 'M18 South-OF-Blue':off((10,27,51.17),(43,54,21.43)),'M18 OF-Red':off((10,27,51.20),(43,54,15.64)),'M18 North-OF-Blue':off((10,27,51.02),(43,54,13.12)),
 'M18 North-Nucleus':off((10,27,51.23),(43,54,14.00)),'M18 South-Nucleus':off((10,27,51.22),(43,54,19.20)),
 'mirror of Soutflow N of S':(-0.15+0.06,-2.55+1.93)}
