# this file contains the geometrical description of the sailing boat
# notations follows the ITTC Symbols 2002, when possible
# all data must be provided at zero speed

###### HULL #######

DIVCAN  1.534 					# [m^3]         Displaced volume of canoe body
LWL     8.100                                   # [m] 		Design waterline length
BWL     1.958                                   # [m] 		Design waterline beam
B       2.920                                   # [m] 		Design maximum beam
AVGFREB 0.6                                     # [m]		Average freeboard
XFB     5.100                                   # [m] 		Longitudinal center of buoyancy LCB from fpp
XFF     5.100                                   # [m]		Longitudinal center of flotation LCF from fpp
CPL     0.5684                                  # [-]		Longitudinal prismatic coefficient
HULLFF  1                                       # [-]		Hull form factor
AW      10.640                                  # [m^2]		Design waterplane area
SC      11.386                                  # [m^2]		Wetted surface's area of canoe body
CMS     0.65   					# [-]		Midship section coefficient
T       1.750                                   # [m]		Total draft
TCAN    0.270                                   # [m]		Draft of canoe body
ALT     5.754                                   # [m]		Total lateral area of yacht
KG      1                                       # [m]		Center of gravity above moulded base or keel
KM      3.167                                   # [m]		Transverse metacentre above moulded base or keel

####### KEEL ########

DVK     0.242                                   # [m^3]		Displaced volume of keel
APK     1.1175                                  # [m^2]		Keel's planform area
ASK     2                                       # [-]		Keel's aspect ratio
SK      2.379                                   # [m^2]		Keel's wetted surface
ZCBK    0.860                                   # [m]		Keel's vertical center of buoyancy (above keel)
CHMEK   0.745                                   # [m]		Mean chord length
CHRTK   1.050                                   # [m]		Root chord length
CHTPK   0.440                                   # [m]		Tip chord length
KEELFF  1                                       # [-]		Keels form factor
DELTTK  0                                       # [-]		Mean thickness ratio of keel section
TAK     0.42                                    # [-]		Taper ratio of keel (CHRTK/CHTPK)

####### RUDDER #######

DVR     0                                       # [m^3]		Rudder's displaced volume
APR     0.3996                                  # [m^2]		Rudder's planform area
SR      0.800                                   # [m^2]		Rudder's wetted surface
CHMER   0.360                                   # [m]		Mean chord length
CHRTR   0.450                                   # [m]		Root chord length
CHTPR   0.270                                   # [m]		Tip chord length
DELTTR  0                                       # [m]		Mean thickness ratio of rudder section	
RUDDFF  1                                       # [-]		Rudder's form factor

####### SAILS ######## 

SAILSET 3  					# [-]		Sails used in calculation  1: main; 3: main+jib; 5: main + spinnaker 
P       11.250                                  # [m]		Mainsail height
E       4.150                                   # [m]		Mainsail base
MROACH  1.2                                     # [-]		Correction for mainsail roach
MFLB    1                                       # [0/1]		Full length battens in main: 0: no, 1: yes
BAD     1.050                                   # [m]		Boom height above deck
I       10.050                                  # [m]		Foretriangle height
J       3.120                                   # [m]		Foretriangle base
LPG     4.400                                   # [m]		Perpendicular of longest jib
SL      10.000                                  # [m]		Spinnaker length
EHM     12.300                                  # [m]		Mast's height above deck
EMDC    0.097                                   # [m]		Mast's average diameter

####### CREW ##########                        

MMVBLCRW 480                                    # [Kg]		Movable Crew Mass
