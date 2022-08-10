import math
import numpy as np

from marginalization import marginalization
from bool import distribution
from resource import Index322


def H(prob):
##################################################################################################################
# - Receives a joint probability distribution;
# - Returns Shannon entropy value;
##################################################################################################################
	S=np.size(prob)# Number of elements
	Pro=np.reshape(prob, S) #One-dimensional vector
	h=np.zeros(S, dtype='float')
	for i in range(0,S):
		if Pro[i]!=0:
			h[i]=-Pro[i]*math.log(Pro[i],2)
	Ent=h.sum() # Shannon entropy
	return Ent



def noisy_ic_edge(n, p44, pd,  pw):
##################################################################################################################
# - Receives parameters' precision n, the three boxes that define the slice of interest p45, pd,  pw
# - Returns the parameters greatest gamma, for each epsilon, of the region pabcDxyz = gamma*p45 +epsilon*pd 
#   + (1-gamma-epsilon)*pw, for which the tripartite IC inequality is not violated;
##################################################################################################################

	E=np.arange(0., (1 + (1/n)), 1/n)
	G=1 - E	


#------------ALL CONDITIONAL PROBABILITIES (CP)------------------------------------------------------------------
	px0x1y0y1=(1/16)*np.ones((2,2,2,2), dtype='float')
	pz = (1/2)*np.ones((2), dtype='float')
	
	p2 = distribution(2)

	p3 = distribution(3)

	#p(Mx|a,x0)
	pMxDax0=p2

	#p(My|b,y0)
	pMyDby0=p2

	#p(x|x0,x1)
	pxDx0x1=p2

	#p(y|y0,y1)
	pyDy0y1=p2

	#p(Zj0|c,Mx,My)
	pZ0Dcd0d1 = p3

	#p(Z1|c,Mx,My)
	pZ1Dcd0d1 = p3


#-----------BISECTION PARAMETERS--------------------------------------------------------------------------------
	iterac = 0
	gammaU = 1 
	gammaB = 0
	#gamma = (gammaU - gammaB)/2
	gamma = 1
	atol=0.001
	bound = 0.0

#-----------VARING PARAMETERS-----------------------------------------------------------------------------------
	ç = n+1
	for i in range(0,ç):
		epsilon=(1/n)*i
		E[i]=epsilon
		max_it = n
		aux = 0
		gammaU = 1-epsilon
		gamma = gammaB
		while(aux < max_it ):
			Viola = False
			iterac = iterac +1

	#-----------SLICE ------------------------------------------------------------------------------------------

			pabcDxyz = gamma*p44 + epsilon*pd + (1-gamma-epsilon)*pw

	#-------------COMPUTING INEQUALITY -------------------------------------------------------------------------

			Viola = noisycriteria_tripartite(pabcDxyz, px0x1y0y1, pz, pMxDax0, pMyDby0, pxDx0x1, pyDy0y1, pZ0Dcd0d1, pZ1Dcd0d1)
			
		
	#-------------TESTING INEQUALITY----------------------------------------------------------------------------
			if Viola:
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
		
	
def noisycriteria_tripartite(pabcDxyz, px0x1y0y1, pz, pMxDax0, pMyDby0, pxDx0x1, pyDy0y1, pZ0Dcd0d1, pZ1Dcd0d1):
##################################################################################################################
# - Receives all joint and conditional distributions
# - Returns the answer if the noisy IC inequality is violated
##################################################################################################################

	Px0x1y0y1MxMyabcd0d1Z0=np.zeros((2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2), dtype='float')

	Px0x1y0y1MxMyabcd0d1Z1=np.zeros((2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2), dtype='float')


	p1=0.5001
	p2=0.5001
	pd0Dmx = np.array([[p1, (1-p1)],[(1-p1), p1]]) # Channel
	pd1Dmy = np.array([[p2, (1-p2)],[(1-p2), p2]]) # Channel
	
#-----------JOINT PROBABILITY---------------------------------------------------------
	
	N = 9
	for k in range(0,(2**N)):
		k = np.binary_repr(k)
		I = np.zeros((N-len(k)), dtype = 'int')
		for l in range(0,len(k)):
			I = np.hstack((I,int(k[l])))
		x0=I[0]
		x1=I[1]
		y0=I[2]
		y1=I[3]
		a=I[4]
		b=I[5]
		c=I[6]
		d0=I[7]
		d1=I[8]
		
		Px0x1y0y1MxMyabcd0d1Z0[x0][x1][y0][y1][(a+x0)%2][(b+y0)%2][a][b][c][d0][d1][(c+d0+d1)%2] = px0x1y0y1[x0][x1][y0][y1]*pMxDax0[(a+x0)%2][a][x0]*pMyDby0[(b+y0)%2][b][y0]*pd0Dmx[d0][(a+x0)%2]*pd1Dmy[d1][(b+y0)%2]*pabcDxyz[Index322((a,b,c,(x0+x1)%2,(y0+y1)%2,0))]*pZ0Dcd0d1[(c+d0+d1)%2][c][d0][d1]
		Px0x1y0y1MxMyabcd0d1Z1[x0][x1][y0][y1][(a+x0)%2][(b+y0)%2][a][b][c][d0][d1][(c+d0+d1)%2] = px0x1y0y1[x0][x1][y0][y1]*pMxDax0[(a+x0)%2][a][x0]*pMyDby0[(b+y0)%2][b][y0]*pd0Dmx[d0][(a+x0)%2]*pd1Dmy[d1][(b+y0)%2]*pabcDxyz[Index322((a,b,c,(x0+x1)%2,(y0+y1)%2,1))]*pZ1Dcd0d1[(c+d0+d1)%2][c][d0][d1]
		
	
	#----------MARGINALIZING-------------------------------------------------------------------
	
	#Maginalizando o que nao faz parte do protocolo	
	Px0x1y0y1MxMyd0d1Z0 = marginalizacao(Px0x1y0y1MxMyabcd0d1Z0, (0,1,2,3,4,5,9,10,11))
	Px0x1y0y1MxMyd0d1Z1 = marginalizacao(Px0x1y0y1MxMyabcd0d1Z1, (0,1,2,3,4,5,9,10,11))
		
#-------------COMPUTING INFORMATIONAL TERMS---------------------------------------------------------------------
	
	IX0 = H(marginalizacao(Px0x1y0y1MxMyd0d1Z0, (0,))) + H(marginalizacao(Px0x1y0y1MxMyd0d1Z0, (2,8))) - H(marginalizacao(Px0x1y0y1MxMyd0d1Z0, (0,2,8)))
	IX1 = H(marginalizacao(Px0x1y0y1MxMyd0d1Z1, (1,))) + H(marginalizacao(Px0x1y0y1MxMyd0d1Z1, (3,8))) - H(marginalizacao(Px0x1y0y1MxMyd0d1Z1, (1,3,8)))
	IY0 = H(marginalizacao(Px0x1y0y1MxMyd0d1Z0, (2,))) + H(marginalizacao(Px0x1y0y1MxMyd0d1Z0, (0,8))) - H(marginalizacao(Px0x1y0y1MxMyd0d1Z0, (0,2,8)))
	IY1 = H(marginalizacao(Px0x1y0y1MxMyd0d1Z1, (3,))) + H(marginalizacao(Px0x1y0y1MxMyd0d1Z1, (1,8))) - H(marginalizacao(Px0x1y0y1MxMyd0d1Z1, (1,3,8)))
	Id0mx = H(marginalizacao(Px0x1y0y1MxMyd0d1Z0, (4,))) + H(marginalizacao(Px0x1y0y1MxMyd0d1Z0, (6,))) - H(marginalizacao(Px0x1y0y1MxMyd0d1Z0, (4,6)))
	Id1my = H(marginalizacao(Px0x1y0y1MxMyd0d1Z0, (5,))) + H(marginalizacao(Px0x1y0y1MxMyd0d1Z0, (7,))) - H(marginalizacao(Px0x1y0y1MxMyd0d1Z0, (5,7)))		
	Ix0x1 = H(marginalizacao(Px0x1y0y1MxMyd0d1Z0, (0,))) + H(marginalizacao(Px0x1y0y1MxMyd0d1Z0, (1,))) - H(marginalizacao(Px0x1y0y1MxMyd0d1Z0, (0,1)))
	Iy0y1 = H(marginalizacao(Px0x1y0y1MxMyd0d1Z0, (2,))) + H(marginalizacao(Px0x1y0y1MxMyd0d1Z0, (3,))) - H(marginalizacao(Px0x1y0y1MxMyd0d1Z0, (2,3)))
	
	desig1 = IX0 + IX1 + Ix0x1 - Id0mx 
	desig2 = IY0 + IY1 + Iy0y1 - Id1my 
	desig = desig1 + desig2

#-------------TESTING INEQUALITY----------------------------------------------
	#if((desig > 10**-9)):
	if((desig1 > 10**-9) or (desig2 > 10**-9)):
		return True
	
	return False # Do not violate
				