import math
from variables import *
from scipy.special import jn

# Function edge to egde center plasma density 
def calculate_hL_hR(ng):
    lambdaI = 1 / (ng * sigmaI)  # 3cm at 1mTorr
    hL = 0.86 * (3 + 0.5 * L / lambdaI)**(-0.5)
    hR = 0.8 * (4 + R / lambdaI)**(-0.5)
    return hL, hR

# Function Electron power balance 
def calculateZind(epsilon0, Lcoil, R, Rc, N, n):
    
    omega = 2*pi/T ## -> usar o valor da figura1 
    omegaPe = ((n * e**2) / (epsilon0 * me))**0.5
    vm = 1 # ve/(S+1) ## ALOOOOOOOOO QUE ISSO?
    epsilonP = 1 - ((omegaPe**2)/(omega*(omega-np.imag(vm)))) #real x imag VM?
    k0 = omega/c 
    k = k0*(epsilonP**0.5)
    
        ### MUDEI O PARAMETRO DE kR para kB, o que seria esse parametro??
    J0 = jn(0, kB)  #Bessel function of the first kind of order zero
    J1 = jn(1, kB)  #Bessel function of the first kind of order one
    Rind = ((2*pi*N*N)/(L*omega*epsilon0))*np.real((1j*R*J1)/(epsilonP*J0))
    Lind = Lcoil*(1-(R**2/Rc**2))+ ((2*pi*N**2)/(L*omega*omega*epsilon0))*np.imag((1j*R*J1)/(epsilonP*J0))
    return Rind, Lind

def effectAreaSurf(hR, hL):
    Aeff = (2*pi*R)*(hR*L + hL*R) #effective area
    Aeff1 = 2*hR*pi*R*L + (2-betaI)*hL*pi*R**2 #effective surface
    return Aeff, Aeff1

#def calculateFlux(n, hL,):
    
    #gamaI = hL*n*uB
    
    #return gamaI, gamaG

#Store particle balance equations in P, whereas P = f(n, ng)
def dPdt(t, P):
    
    n, ng, Tg, Te = P 
    #-------------------------------------n
    #edge to center plasma density
    hL, hR = calculate_hL_hR(P[1]) # indicar P[0]
    
    #effecive area and surface
    Aeff, Aeff1 = effectAreaSurf(hR, hL)
    
    #Bohm speed
    uB = (kB*Te/M)**0.5 #Bohm speed
    
    vg = ((8*kB*Tg)/pi*M)**0.5#mean velocity of the atoms
    
    #ion flux leaving thruster through grid holes
    gamaI = hL*n*uB 
    #Ji = 
    #if 
    #thermal flux of NEUTRAL gas across Ag
    gamaG = 0.25*ng*vg
    
    #Constants
    Kiz1 = 6.73*(10^-15)*((kB*Te/e)**0.5)*(3.97+0.643*(kB*Te/e)-0.0368*(kB*Te/e)**2)*math.exp((-e*Eiz)/(kB*Te)) #m³s-1
    Kiz2 = 6.73*(10^-15)*((kB*Te/e)**0.5)*(-0.0001031*(kB*Te/e)**2 + 6.386*np.exp((-e*Eiz)/(kB*Te))) #m³s-1
    #Electron-impact ionization
    Kiz = (Kiz1+Kiz2)/Kiz2 
    
    #-------------------------------------ng
    #m³s-1 Electron-impact excitation
    Kex = 1.2921*(10^-13)*math.exp(-e*Eexc/kB*Te) 
    
    
    #-------------------------------------Tg
    #mean velocity of the atoms
    
    
    # rate for ion-neutral elastic collisions
    Kin = sigmaI*vg  
    #heat diffusion length 
    lambda0 = (R/2.405) + (L/pi)
    
    #-------------------------------------Te
    Rind, Lind  = calculateZind( epsilon0, Lcoil, R, Rc, N, n)
    
    #Ploss and Pabs
    Ploss = Eiz*n*ng*Kiz + Eexc*n*ng*Kex + 3*(me/M)*kB*(Te-Tg)*n*ng*Kel + (7*kB*Te)*n*uB*(Aeff/V)
    Pabs = 0.5*Rind*Icoil**2/V
    
    return [n*ng*Kiz - n*uB*(Aeff/V), # n
            (Q0/V) + n*uB*(Aeff1/V) - n*ng*Kiz - gamaG*(Ag/V), #ng 
            3*(me/M)*kB*(Te-Tg)*n*ng*Kel + 0.25*M*n*ng*Kin*uB**2 - k_cte*((Tg-Tg0)/lamb0)*(A/V), #Tg
            Pabs - Ploss] #Te
   


