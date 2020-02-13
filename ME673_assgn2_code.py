#!/usr/bin/python
''' 
 Python code for Numerical assignment 2, ME 678
 Written by Gowtham Kuntumalla, 140100091
 Fanno flow analysis.
 Both subsonic and supersonic inlet conditions are considered.
 M, P, T are the required properties
 
'''
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.optimize import fsolve


print("hey Code is running !")
# INLET CONDITIONS
gamma=1.4
R=287
T0=500 # K
P0=101000 # Pa
rho0=P0/(R*T0) #kg/m3 
D = 0.05 # Pipe diameter
M0 = 0.2 # inlet mach
h = 0.1 # x1=x0+h
fr = 0.005  # Friction factor


# Fanno flow DE #
def dlnMdx(x,M):
	return gamma * (M**2)*(2+(gamma-1)*M**2)/(1-M**2)*(fr/D)
def dlnPdx(x,M):	
	return -2*gamma * (M**2)*(1+(gamma-1)*M**2)/(1-M**2)*(fr/D)
def dlnTdx(x,M):
	return -2*gamma*(gamma-1)*(M**4)/(1-M**2)*(fr/D)	
	

# RK METHOD Functions #
def rk2(f,x,y):
	k1 = h*f(x,y)
	k2 = h*f(x+h,y+k1)
	return 0.5*(k1+k2)
def rk4(f,x,y):
	k1 = h*f(x,y)
	k2 = h*f(x+h/2,y+k1/2)
	k3 = h*f(x+h/2,y+k2/2)
	k4 = h*f(x+h,y+k3)
	return 1/6*(k1+2*k2+2*k3+k4)	
		 
if __name__ == "__main__":
	rk_method = rk2 # Use either rk2 or rk4
	x=[]; i=0
	M=[]; P=[]; T=[]; # individual lists of properties
	x.append(0);M.append(M0); P.append(P0); T.append(T0);
	while M[len(M)-1] < 0.8: # required exit mach. edit it
		M.append(M[i]*math.exp(rk_method(dlnMdx,x[i],M[i])))
		P.append(P[i]*math.exp(rk_method(dlnPdx,x[i],M[i])))
		T.append(T[i]*math.exp(rk_method(dlnTdx,x[i],M[i])))
		x.append(x[i]+h)
		i+=1
	print (M,P,T)
	# plotting
	f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
	ax1.plot(x, M)
	ax1.set_ylabel('Mach Number')
	ax1.set_title('Variation of properties along length')
	ax2.plot(x, P, color='g')
	ax2.set_ylabel('Pressure (Pa)')
	ax3.plot(x, T, color='r')
	ax3.set_ylabel('Temperature (K)')
	ax3.set_xlabel('Distance along the pipe (m)')	
	plt.show()
	
	print("Analysis done. Look at the plots")

#END OF PROGRAM