#import the necessary functions from scipy and python's base libraries
from scipy.optimize import fsolve
import math

#CONSTANTS
###################################################################################################
g = 9.81 #m/s^2
rho = 997 #kg/m^3
mu = 0.001 #Pa*s
eps = 0.0024 #pipeflow.com
d = 0.00794 #m
dT = 0.005 #delta time (s)
Atank = 0.32 * 0.26 #Area of tank
###################################################################################################

#Get the Reynold's number of the flow given a velocity
def getRe(v):
	return (rho*v*d)/mu

#This is the function that will be used to converge on an f value using fsolve. 
#This function is taken from the Moody diagram
#Inputs:
#x = 1/sqrt(f)
#data = a tuple containing the current Reynold's number of the flow
def fEq(x, *data):
	Re = data[0]
	return x + 2 * math.log((eps/d/3.7) + (2.51/Re) * x)

#Our derived expression for finding velocity, given pipe length, friction factor and current water level
#L = pipe length
#f = friction factor
#wL = water level
def findV(L, f, wL):
	return math.sqrt((2*g*(wL + L / 150)) / (0.5 + L * f / d))

#Function to find friction factor given a turbulent flow
#fguess = a starting guess for friction factor
#waterLevel = current height of water in the tank
#L = pipe length
def findFTurbo(fguess, waterLevel, L):
	f = fguess
	fold = 98465904856

	while (abs(f - fold) > 0.0001):
		fold = f
		v = findV(L, f, waterLevel)
		Re = getRe(v)

		data = tuple([Re])
		x = fsolve(fEq, 0, args=data)[0]
		f = (1/x)**2

	return f, v, Re

#Function to find friction factor given a laminar flow
#fguess = a starting guess for friction factor
#waterLevel = current height of water in the tank
#L = pipe length
def findFLammy(fguess, waterLevel, L):
	f = fguess
	fold = 98465904856

	while (abs(f - fold) > 0.0001):
		fold = f
		v = findV(L, f, waterLevel)
		Re = getRe(v)
		f = 64/Re

	return f, v, Re

#main function
#this function runs a loop that tests different pipe lengths and reports their time-to-drain to the output console
def main():
	L = 0.05 #m

	#test pipe lengths between 5cm and 1m inclusive, in increments of 5cm
	while(L <= 1.001):
		WL = 0.1 #initial water height of 10cm
		fguess = 0.019 #initial guess for friction factor
		t = 0
		ReOld = 57609456

		#continue running the simulation until the water level dips below 2cm above the pipe entrance
		while(WL > 0.02):
			#check the previous Reynold's number to see whether our flow is transitioning from laminar to turbulent
			#and use the appropriate findF function
			if (ReOld > 2300):
				f, v, Re = findFTurbo(fguess, WL, L)

			else:
				f, v, Re = findFLammy(fguess, WL, L)

			q = v * math.pi * (d/2)**2 #volumetric flow rate
			WL = WL - (q * dT) / Atank #calculate new water level
			t = t + dT #increment time
			ReOld = Re #update Re

		#after the simulation is complete for a given pipe length, output the time-to-drain and increment pipe length
		print("Length: ", L, "Time: ", t)
		L = L + 0.05

#start the script
main()