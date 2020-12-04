#import the necessary libraries
from scipy.optimize import fsolve
import math

#CONSTANTS
###################################################################################################
g = 9.81 #m/s^2
rho = 997 #kg/m^3
mu = 0.001 #Pa*s
eps = 0.0024 #pipeflow.com
dP = 0.00794 #m
dTJ = 0.0111125 #m
dT = 0.005 #delta time (s)
Atank = 0.32 * 0.26
###################################################################################################

#Get the Reynold's number of the flow given a velocity
def getRe(v, d):
	return (rho*v*d)/mu

#This is the function that will be used to converge on an f value using fsolve. 
#This function is taken from the Moody diagram
#Inputs:
#x = 1/sqrt(f)
#data = a tuple containing the current Reynold's number of the flow, as well as the pipe diameter
def fEq(x, *data):
	Re, d = data
	return x + 2 * math.log((eps/d/3.7) + (2.51/Re) * x)

#Our derived expression for finding velocity with the T-joint attached, given pipe length, friction factor and current water level
#L = pipe length
#fp = friction factor of the pipe
#ft = friction factor of the T-joint
#wL = water level
def findVTPipe(L, fp, ft, wL):
	return math.sqrt((2 * g * (wL + L / 150))/(125.94 * L * fp + 0.234 * ft + 1.565))

#Function to find friction factor given a laminar flow
#fguessPipe = a starting guess for friction factor of the pipe
#fguessTJoint = a starting guess for the friction factor of the T-joint
#waterLevel = current height of water in the tank
#L = pipe length
def findFTPipeLammy(fguessPipe, fguessTJoint, waterlevel, L):
	#at first, assign the friction factors of the pipe and T-joint to our initial guesses
	fPipe = fguessPipe
	fTJ = fguessTJoint

	#initially set our "old" f values to a random number since this is the first iteration
	fPipeOld = 8340983406
	fTJOld = 9384093460

	#while the f values haven't yet converged
	while (abs(fPipe - fPipeOld) > 0.001 and abs(fTJ - fTJOld) > 0.001):
		#calculate velocity, friction factor and Reynold's number for both pipe and T-joint each iteration
		fPipeOld = fPipe
		vPipe = findVTPipe(L, fPipe, fTJ, waterlevel)
		RePipe = getRe(vPipe, dP)
		fTJOld = fTJ
		vTJ = 0.255 * vPipe
		ReTJ = getRe(vTJ, dTJ)

		#f = 64/Re for laminar flows
		fTJ = 64/ReTJ
		fPipe = 64/RePipe

	return fPipe, vPipe, RePipe, fTJ, vTJ, ReTJ

#Function to find friction factor given a turbulent flow
#fguessPipe = a starting guess for friction factor of the pipe
#fguessTJoint = a starting guess for the friction factor of the T-joint
#waterLevel = current height of water in the tank
#L = pipe length
def findFTPipeTurbo(fguessPipe, fguessTJoint, waterlevel, L):
	#at first, assign the friction factors of the pipe and T-joint to our initial guesses
	fPipe = fguessPipe
	fTJ = fguessTJoint

	#initially set our "old" f values to a random number since this is the first iteration
	fPipeOld = 8340983406
	fTJOld = 9384093460

	#while the friction factors have not converged
	while (abs(fPipe - fPipeOld) > 0.001 and abs(fTJ - fTJOld) > 0.001):
		
		#calculate friction factors, velocities and Reynold's numbers each iteration
		fPipeOld = fPipe
		fTJOld = fTJ
		vPipe = findVTPipe(L, fPipe, fTJ, waterlevel)
		vTJ = 0.255 * vPipe
		RePipe = getRe(vPipe, dP)
		ReTJ = getRe(vTJ, dTJ)

		dataP = tuple([RePipe, dP])
		xP = fsolve(fEq, 0, args=dataP)[0]
		fPipe = (1/xP)**2

		dataTJ = tuple([ReTJ, dTJ])
		xTJ = fsolve(fEq, 0, args=dataTJ)[0]
		fTJ = (1/xTJ)**2

	return fPipe, vPipe, RePipe, fTJ, vTJ, ReTJ

#main function
#this function runs a loop that tests different pipe lengths and reports their time-to-drain to the output console
def main():
	L = 0.05 #m

	#test pipe lengths between 5cm and 1m inclusive, in increments of 5cm
	while(L <= 1.001):
		WL = 0.1
		fgp = 0.03
		fgtj = 0.03
		t = 0
		ReOld = 57609456

		#simulate the current pipe length until the water level dips below 2cm above the pipe
		while(WL > 0.02):
			#check the previous Reynold's number to see whether our flow is transitioning from laminar to turbulent
			#and use the appropriate findF function
			if (ReOld > 2300):
				fP, vP, ReP, fTJ, vTJ, ReTJ = findFTPipeTurbo(fgp, fgtj, WL, L)

			else:
				fP, vP, ReP, fTJ, vTJ, ReTJ = findFTPipeLammy(fgp, fgtj, WL, L)

			q = vP * math.pi * (dP/2)**2 #volumetric flow rate
			WL = WL - (q * dT) / Atank #update water level for this iteration
			t = t + dT #increment time
			ReOld = ReP #update Re
			
			#update our guess for friction factor each iteration in order to converge faster
			fgp = fP 
			fgtj = fTJ

		#after the simulation is complete for a given pipe length, output the time-to-drain and increment pipe length
		print("Length: ", L, "Time: ", t)
		L = L + 0.05

#start the script
main()
