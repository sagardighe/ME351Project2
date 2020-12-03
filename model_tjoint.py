from scipy.optimize import fsolve
import math
import cmath

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

def getDWeighted(L):
	tubeWeight = L / (L+0.04)
	tJointWeight = 0.04 / (L+0.04)

	return tubeWeight * dP + tJointWeight * dTJ

def getRe(v, d):
	return (rho*v*d)/mu

def fEq(x, *data):
	Re, d = data
	return x + 2 * math.log((eps/d/3.7) + (2.51/Re) * x)

def findVTPipe(L, fp, wL):
	D = getDWeighted(L)
	return math.sqrt((2 * g * (wL + L / 150)) / ((fp * (L + 0.0102)**2) / (D * (L + 0.04)) + 1.565))

def findFTPipeLammy(fguessPipe, waterlevel, L):
	fW = fguessPipe
	fWOld = 8340983406

	while (abs(fW - fWOld) > 0.001):
		fWOld = fW
		vW = findVTPipe(L, fW, waterlevel)
		ReW = getRe(vW, dP)

		fW = 64/ReW

	return fW, vW, ReW

def findFTPipeTurbo(fguessPipe, waterlevel, L):
	fW = fguessPipe
	fWOld = 8345983459834

	while(abs(fW - fWOld) > 0.001):
		fWOld = fW
		vWeighted = findVTPipe(L, fW, waterlevel)
		dWeighted = getDWeighted(L)
		ReWeighted = getRe(vWeighted, dWeighted)

		dataWeighted = tuple([ReWeighted, dWeighted])
		xW = fsolve(fEq, 0, args=dataWeighted)[0]
		fW = (1/xW)**2

	return fW, vWeighted, ReWeighted

def main():
	L = 0.3 #m

	while(L <= 0.3):
		WL = 0.1
		fgp = 0.019
		t = 0
		ReOld = 57609456

		while(WL > 0.02):
			if (ReOld > 2300):
				fP, vP, ReP = findFTPipeTurbo(fgp, WL, L)

			else:
				fP, vP, ReP = findFTPipeLammy(fgp, WL, L)

			q = vP * math.pi * (dP/2)**2
			WL = WL - (q * dT) / Atank
			t = t + dT
			ReOld = ReP
			fgp = fP

			print(fP, vP, ReP, WL, q, t)

		print("Length: ", L, "Time: ", t)
		L = L + 0.05

main()
