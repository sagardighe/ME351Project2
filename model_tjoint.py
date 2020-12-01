from scipy.optimize import fsolve
import math

#CONSTANTS
###################################################################################################
g = 9.81 #m/s^2
rho = 997 #kg/m^3
mu = 0.0008891 #Pa*s
eps = 0.0015 #pipeflow.com
dP = 0.00794 #m
dTJ = 0.0111125 #m
dT = 0.005 #delta time (s)
Atank = 0.32 * 0.26
###################################################################################################

def getRe(v, d):
	return (rho*v*d)/mu

def fEq(x, *data):
	Re, d = data
	return x + 2 * math.log((eps/d/3.7) + (2.51/Re) * x)

def findVTPipe(L, fp, ft, wL):
	#return math.sqrt((-2*g*(wL - 0.02 + L / 150))/(0.805 - L/(0.03176*fp) - 5.85e-2*ft))
	#print(fp,ft)
	#return math.sqrt((2*g*(wL - 0.02 + L / 150))/(125.945 * L * fp + 0.2341 * ft - 0.305))
	return (118816*math.sqrt(g*(0.00666667*L+wL-0.02)))/(math.sqrt(889000000000*fp*L + 1652155200*ft - 2611351267))

def findFTPipeLammy(fguessPipe, fguessTJoint, waterlevel, L):
	fPipe = fguessPipe
	fTJ = fguessTJoint
	fPipeOld = 8340983406
	fTJOld = 9384093460

	while (abs(fPipe - fPipeOld) > 0.001 and abs(fTJ - fTJOld) > 0.001):
		fPipeOld = fPipe
		vPipe = findVTPipe(L, fPipe, fTJ, waterlevel)
		RePipe = getRe(vPipe, dP)
		fTJOld = fTJ
		vTJ = 0.255 * vPipe
		ReTJ = getRe(vTJ, dTJ)

		fTJ = 64/ReTJ
		fPipe = 64/RePipe

	# while (abs(fTJ - fTJOld) > 0.0001):
	# 	fTJOld = fTJ
	# 	vPipe = findVTPipe(L, fPipe, fTJ, waterlevel)
	# 	vTJ = 0.255 * vPipe
	# 	ReTJ = getRe(vTJ, dTJ)

	# 	fTJ = 64/ReTJ

	return fPipe, vPipe, RePipe, fTJ, vTJ, ReTJ

def findFTPipeTurbo(fguessPipe, fguessTJoint, waterlevel, L):
	fPipe = fguessPipe
	fTJ = fguessTJoint
	fPipeOld = 8340983406
	fTJOld = 9384093460

	while (abs(fPipe - fPipeOld) > 0.001 and abs(fTJ - fTJOld) > 0.001):
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

	# while (abs(fTJ - fTJOld) > 0.0001):
	# 	fTJOld = fTJ
	# 	vPipe = findVTPipe(L, fPipe, fTJ, waterlevel)
	# 	vTJ = 0.255 * vPipe
	# 	ReTJ = getRe(vTJ, dTJ)

	# 	dataTJ = tuple([ReTJ, dTJ])
	# 	xTJ = fsolve(fEq, 0, args=dataTJ)[0]
	# 	fTJ = (1/xTJ)**2

	return fPipe, vPipe, RePipe, fTJ, vTJ, ReTJ

def main():
	L = 0.2 #m

	while(L <= 1.01):
		WL = 0.1
		fgp = 0.03
		fgtj = 0.03
		t = 0
		ReOld = 57609456

		while(WL > 0.02):
			if (ReOld > 2300):
				fP, vP, ReP, fTJ, vTJ, ReTJ = findFTPipeTurbo(fgp, fgtj, WL, L)

			else:
				fP, vP, ReP, fTJ, vTJ, ReTJ = findFTPipeLammy(fgp, fgtj, WL, L)

			q = vP * math.pi * (dP/2)**2
			WL = WL - (q * dT) / Atank
			t = t + dT
			ReOld = ReP
			fgp = fP
			fgtj = fTJ
			#print("Pipe:", fP, "TJoint", fTJ)
			#print(fP, vP, ReP, fTJ, vTJ, ReTJ, WL, q, t)

		print("Length: ", L, "Time: ", t)
		L = L + 0.05

main()
