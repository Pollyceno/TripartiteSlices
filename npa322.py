import ncpol2sdpa as nc
from recurso import *
from NCP322 import *

def NPA322(n, p44, pd,  pw, level):
##################################################################################################################
# - Receives parameters' precision n, the three boxes that define the slice of interest p45, pd,  pw, 
#   and hierarchy level
# - Returns the parameters greatest gamma, for each epsilon, of the region pabcDxyz = gamma*p45 +epsilon*pd 
#   + (1-gamma-epsilon)*pw, inside of Q_2;
##################################################################################################################
		
	E=np.arange(0, 1+(1/n), 1/n, dtype='float')
	G= 1 - E

##################################################################################################################
##################### BISECTION PARAMETERS #######################################################################
	iterac = 0
	gammaU = 1 
	gammaB = 0
	gamma = (gammaU - gammaB)/2
	atol=10**-8
##################################################################################################################
##################### VARING PARAMETERS ##########################################################################

	ç = n+1
	for i in range(0,ç):
		epsilon=(1/n)*i
		E[i]=epsilon
		max_it = 15
		aux = 0
		gammaU = 1-epsilon
		gamma = gammaB
		while(aux < max_it ):
			iterac = iterac +1

##################################################################################################################
##################### SLICE ######################################################################################

			pabcDxyz = gamma*p44 + epsilon*pd + (1-gamma-epsilon)*pw

##################################################################################################################
##################### Verifying Q_2 ##############################################################################
			Q = quantumedge322(pabcDxyz, level)


			if(Q == 'dual_infeas_cer'):
				gammaU = gamma
				if(abs(gammaU - gammaB) > atol):#Up to the precision atol, devide the current gamma's range in half
					gamma = gammaU - (gammaU - gammaB)/2
				else:
					break
			else:
				gammaB = gamma
				if(abs(gammaU - gammaB) > atol):#Up to the precision atol, devide the current gamma's range in half
					gamma = gammaB + (gammaU - gammaB)/2
				else:
					break
			aux = aux +1
	
		aux = 0
	
		gammaB = 0.0 #reset
		G[i] = gamma # Save the best gamma

	return (E, G)