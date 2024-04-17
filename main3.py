from variables import *
from equations2 import *
import numpy as np 
from scipy.integrate import odeint


P0 = (n0, ng0, Tg0, Te0)

t = np.linspace(0,1,100)

sol = odeint(dPdt, y0=P0,t=t, tfirst=True)



Ji = e*gamaI 

## parametros de indicaçao de steady sate?

#Thrust and efficiency quantities ions
Ti =gamaI*M*vBeam*Ai

#ion thruhtst power

Pi = 0.5*M*gamaI*Ai*vBeam**2

##current density 

##where Ji < Jcl

#Neutrosn Thrust and power
Tn = gamaG*M*vg*Ag
Pn = 0.25&M*gamaG*Ag*vg**2
