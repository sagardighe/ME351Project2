from scipy.optimize import fsolve
import math

#CONSTANTS
###################################################################################################
g = 9.81 #m/s^2
rho = 997 #kg/m^3
mu = 0.0008891 #Pa*s
eps = 0.0015 #pipeflow.com
d = 0.00794 #m
dT = 0.005 #delta time (s)
Atank = 0.32 * 0.26
###################################################################################################

def getRe(v):
	return (rho*v*d)/mu

def fEq(x, *data):
	Re = data[0]
	return x + 2 * math.log((eps/d/3.7) + (2.51/Re) * x)

def findV(L, f, wL):
	return math.sqrt((2*g*(wL - 0.02 + L / 150)) / (0.5 + L * f / d))


def findFTPipe(fguessPipe, fguessTJoint, waterlevel, L):
	fPipe = fguessPipe
	fTJ = fguessTJoint

	while (abs(fPipe - fPipeOld) > 0.0001 and abs(fTJ - fTJOld) > 0.0001):
		fPipeOld = fPipe

def findVTPipe(L, f, wL):
	return math.sqrt((-2*g*())/())

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

def findFLammy(fguess, waterLevel, L):
	f = fguess
	fold = 98465904856

	while (abs(f - fold) > 0.0001):
		fold = f
		v = findV(L, f, waterLevel)
		Re = getRe(v)
		f = 64/Re

	return f, v, Re

def main():
	L = 0.3 #m

	while(L <= 0.3):
		WL = 0.1
		fguess = 0.019
		t = 0
		ReOld = 57609456

		while(WL > 0.02):
			if (ReOld > 2300):
				f, v, Re = findFTurbo(fguess, WL, L)

			else:
				f, v, Re = findFLammy(fguess, WL, L)

			q = v * math.pi * (d/2)**2
			WL = WL - (q * dT) / Atank
			t = t + dT
			ReOld = Re
			print(f, v, q, WL, Re)

		print("Length: ", L, "Time: ", t)
		L = L + 0.01


main()