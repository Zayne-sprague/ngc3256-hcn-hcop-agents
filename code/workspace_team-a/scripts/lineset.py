# Target line list (rest GHz from Splatalogue CDMS/JPL; checked with astroquery at fetch time) and region list.
# name, rest GHz, band-OUS tuning, flags, half-window km/s for integration, baseline order, notes
VSYS=2780.0   # km/s, radio convention, LSRK (moment-1 of HCO+ 3-2 at both nuclei, Xb4f: 2780.2 / 2780.1)
RA0,DE0=156.96346,-43.90461   # phase centre (= all pointings)
# regions: offsets (dRA*,dDec) arcsec from phase centre, defined from Xb4d HCO+ 3-2 / HCN 3-2 moment-0 and wing peaks (regions.py)
REGIONS={'N':(0.00,2.55),'S':(-0.15,-2.55),'Bridge':(0.30,0.15),'NEclump':(6.75,3.90),'Eclump':(5.25,-2.85),
         'Soutflow':(-0.15,-4.50),'Sext':(-1.05,-5.40),'SWarm':(-3.75,-4.50)}
# (name, rest_GHz, W, blends)
TARGETS={
 'B5':[('HCN 2-1',177.261111,350,''),('HCO+ 2-1',178.375056,350,''),('HOC+ 2-1',178.972051,250,'near NH2CN 9(1,9)-8(1,8) 178.749 (-190 km/s)'),
       ('CH3CHO 9(1,8)-8(1,7) E',177.543458,200,'weak, Eu=45K'),('SO 4_5-3_4',178.605403,200,'Eu=24K; 4 ch from HOC+ not blended (-190 vs +..)'),
       ('CCS 13-12',166.662354,200,''),('CH3CN 9(0)-8(0)',165.569092,200,'K=0..3 ladder within 85 MHz (+-150 km/s blends)'),
       ('unid 163.54 GHz (rest~165.07)',165.0707,250,'unidentified')],
 'B7':[('HCN 4-3',354.505478,350,'13CH3OH 4(1,3)-3(0,3) at -50 km/s'),('HCO+ 4-3',356.734223,350,'SO2 10(4,6)-10(3,7) +18 km/s; CH3OCH3 8(4,x)-7(3,x) -90..-40 km/s'),
       ('CS 7-6',342.883000,300,'H2CS 10(0,10)-9(0,9) 342.9444 at +54 km/s (Eu=91K, logA -3.2) blends'),('H2CS 10(0,10)-9(0,9)',342.944369,150,'blended with CS 7-6 wing'),
       ('H13CN 4-3',345.339760,300,'SO2 13(2,12)-12(1,11) +0.5 km/s?; NH2CHO 16(1,15)-15(1,14) -10 km/s; CO 3-2 wing slope (edge of spw20)'),
       ('HC15N 4-3',344.200320,250,''),('SO 8_8-7_7',344.310612,250,'Eu=87K'),
       ('HCN v2=1f 4-3',356.255550,250,'vib. excited; Eu~1050K'),('H2C17O 5(2,4)-4(2,3)',354.908411,150,''),
       ('CO+ 3-2 J=7/2-5/2',354.014250,150,'tentative'),],
 'B6':[('HCO+ 3-2',267.557626,350,'SO2 13(3,11)-13(2,12) -22 km/s Eu=106K; OCS 22-21 -?'),('HCN 3-2',265.886434,350,''),('HOC+ 3-2',268.451094,250,'near CH3CHO 14(2,13)-13(2,12) 268.30 (-170 km/s)'),
       ('HCN v2=1f 3-2',267.199390,250,'vib. excited'),('HC3N 28-27',254.699500,250,'Eu=177K'),('SO 6_5-5_4',251.825770,250,'proposal SO(6_5-5_4); CH3OH 5(3)-5(2) A/E cluster within -+30 km/s (Eu 85K, logA -4.2)'),
       ('34SO 6(6)-5(5)',253.208020,200,''),('NS 11/2-9/2 f',253.968393,200,'hyperfine pair 253.9684/253.9706'),('NS 11/2-9/2 e',253.570476,200,''),('SiS 14-13',254.103211,200,''),('H2COH+ 4(0,4)-3(0,3)',252.870339,200,''),
       ('CH2NH 4(1,3)-3(1,2)',266.270020,200,''),('H2C18O 4(1,4)-3(1,3)',268.745789,200,''),('NH2CHO 13(0,13)-12(0,12)',267.063393,200,''),
       ('CH3OH 5(2)-4(1) E',266.838148,200,''),('CH3OH 3(3)-3(2) A cluster',251.905729,200,'CH3OH J_3-J_2 A-type Q-branch cluster 251.36-252.8 (blends w/ SO 5(6)-4(5) 251.826)'),
       ('l-C3H 23/2-21/2 f',251.434000,200,''),('CH3CHO 13(3,10)-12(3,9) E',251.363729,200,'')],
}
# which cubes (tuning->spw) contain which target (by sky frequency) is decided automatically from the cube frequency axis.
