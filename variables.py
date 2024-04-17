import math
import numpy as np
#-------------------------------------- Initial Variables
# Constants
pi = np.pi  #GLOBAL
e = 1.602176634e-19 #elemental eletron
me = 9.1093837e-31 #mass electron
M = 131  # Xenon ion mass in atomic mass units (uma)
kB = 1.380649e-23 #boltz constant
sigmaI = 1e-18  # sigmaI meters² globalCossSection
k_cte = 0.0057  # W/K xenon gas thermal conductivity
epsilon0 = 8.854e-12 #F/m permittivity of vacuum
mi0 = 4e-7 * pi
lc = 1 #total length
c = 299792458 #m/s


# Initial conditions
n0 = 1 # Initial plasma density
ng0 = 1 # Initial neutral gas density
Te0 = 1 # Initial electron temperature
Tg0 = 1 # Initial neutral gas temperature 300K
vg0 = 1


# Initial Variables
R = 0.06  # meters Radius
L = 0.10  # meters Length
Rc = 0.07  # mmeters Radius Coil
T = 0.5


Q0 = 1.2e19  # s^-1 Initial gas injection
p = 2  # mTorr Initial pressure

betaI = 0.7  # ionsGridTransparency
betaG = 0.3  # neutralAtomsTrasnparency   0 < betaG < 1
Vgrid = 1000  # Volts voltageDiff
V = Vgrid
s = 0.001  # meters gapGrids
Jcl = (4/9)*epsilon0*((2*e/M)**0.5)*((Vgrid**1.5)/s**2) #Child-Langmuir limit current density

N = 5  # turns numberTurns
Lcoil = 4.84  # miH inductanceCoil ###(mi0*pi*Rc*Rc*N*N)/lc #coil in ducatance ### Perguntar
Rcoil = 2  # ohm resistanceCoil
Icoil = 26 # Amperes

lamb0 = (R/2.405) + (L/pi) # heat diffusion length

Eexc = 11.6  # V
Eiz = 12.127  # V
Kel = 10^-13 #m³s-1 Elastic scattering 



# Global variables
Vplasma = pi * R * R * L  # plasma volume
A = 2 * pi * R * (R + L)  # total Inside area
Ag = betaG * pi * R * R  # open area of thruster (NEUTRALS)
Ai = betaI * pi * R * R  # open area ions scape where, AI > AG


Q0 = 0.25*ng0*vg0*Ag #plasma turned off, gas injection rate 
#ng0 and vg0 are neutran gas dens and mean vel withouh plasma 

p0 = (4*kB*Tg0*Q0)/(vg0*Ag) #pressure in chamber 
vBeam = math.sqrt((2*e*Vgrid)/M)


