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


def IC_single_copy(n, p45, pd,  pw):
##################################################################################################################
# - Receives parameters' precision n, the three boxes that define the slice of interest p45, pd,  pw
# - Returns the parameters greatest gamma, for each epsilon, of the region pabcDxyz = gamma*p45 +epsilon*pd 
#   + (1-gamma-epsilon)*pw, for which the tripartite IC inequality is not violated;
##################################################################################################################
	
	E=np.arange(0., (1 + (1/n)), 1/n)
	G=1 - E	

	Px0x1y0y1MxMyabcxyZ0=np.zeros((2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2), dtype='float')
	Px0x1y0y1MxMyabcxyZ1=np.zeros((2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2), dtype='float')


#------------ALL CONDITIONAL PROBABILITIES (CP)------------------------------------------------------------------
	px0x1y0y1=(1/16)*np.ones((2,2,2,2), dtype='float')#Uniform dist
	pz = (1/2)*np.ones((2), dtype='float')#Uniform dist
	
	p2 = distribution(2)#CP for f(x,y) = x+y of two inputs variable

	p3 = distribution(3)#CP for f(x,y,z) = x+y+z of three inputs variable

	#p(Mx|a,x0)
	pMxDax0=p2

	#p(My|b,y0)
	pMyDby0=p2

	#p(x|x0,x1)
	pxDx0x1=p2

	#p(y|y0,y1)
	pyDy0y1=p2

	#p(Zj0|c,Mx,My)
	pZ0DcMxMy = p3

	#p(Z1|c,Mx,My)
	pZ1DcMxMy = p3


#-----------BISECTION PARAMETERS--------------------------------------------------------------------------------
	iterac = 0
	gammaU = 1 
	gammaB = 0
	gamma = (gammaU - gammaB)/2
	atol=0.001
	bound = 0.0

#-----------VARING PARAMETERS------------------------------------------------------------
	ç = n+1
	for i in range(0,ç):
		epsilon=(1/n)*i
		E[i]=epsilon
		max_it = 10
		aux = 0
		gammaU = 1-epsilon
		gamma = gammaB
		while(aux < max_it ):
			iterac = iterac +1

	#-----------SLICE --------------------------------------------------------------------

			pabcDxyz = gamma*p45 + epsilon*pd + (1-gamma-epsilon)*pw
			
		#-----------JOINT PROBABILITY-----------------------------------------------------
			
			N = 12
			for k in range(0,(2**N)):
				k = np.binary_repr(k)
				I = np.zeros((N-len(k)), dtype = 'int')
				for l in range(0,len(k)):
					I = np.hstack((I,int(k[l])))
				x0=I[0]
				x1=I[1]
				y0=I[2]
				y1=I[3]
				Mx=I[4]
				My=I[5]
				a=I[6]
				b=I[7]
				c=I[8]
				x=I[9]
				y=I[10]
				#z=I[10]
				Z0=I[11]
				Z1=I[11]
				
				Px0x1y0y1MxMyabcxyZ0[x0][x1][y0][y1][Mx][My][a][b][c][x][y][Z0] = px0x1y0y1[x0][x1][y0][y1]*pMxDax0[Mx][a][x0]*pMyDby0[My][b][y0]*pabcDxyz[Index322(a,b,c,x,y,0)]*pxDx0x1[x][x0][x1]*pyDy0y1[y][y0][y1]*pz[0]*pZ0DcMxMy[Z0][c][Mx][My]/pz[0]
				Px0x1y0y1MxMyabcxyZ1[x0][x1][y0][y1][Mx][My][a][b][c][x][y][Z1] = px0x1y0y1[x0][x1][y0][y1]*pMxDax0[Mx][a][x0]*pMyDby0[My][b][y0]*pabcDxyz[Index322(a,b,c,x,y,1)]*pxDx0x1[x][x0][x1]*pyDy0y1[y][y0][y1]*pz[1]*pZ1DcMxMy[Z1][c][Mx][My]/pz[1]
					
			#----------MARGINALIZING-------------------------------------------------------------------
			
			Px0x1y0y1MxMyZ0 = marginalizacao(Px0x1y0y1MxMyabcxyZ0, (0,1,2,3,4,5,11))
			Px0x1y0y1MxMyZ1 = marginalizacao(Px0x1y0y1MxMyabcxyZ1, (0,1,2,3,4,5,11))
		

	#-------------COMPUTING INFORMATIONAL TERMS---------------------------------------------------------------------

			IX0 = H(marginalizacao(Px0x1y0y1MxMyZ0, (0,))) + H(marginalizacao(Px0x1y0y1MxMyZ0, (2,6))) - H(marginalizacao(Px0x1y0y1MxMyZ0, (0,2,6)))
	
			IX1 = H(marginalizacao(Px0x1y0y1MxMyZ1, (1,))) + H(marginalizacao(Px0x1y0y1MxMyZ1, (3,6))) - H(marginalizacao(Px0x1y0y1MxMyZ1, (1,3,6)))
	
			IY0 = H(marginalizacao(Px0x1y0y1MxMyZ0, (2,))) + H(marginalizacao(Px0x1y0y1MxMyZ0, (0,6))) - H(marginalizacao(Px0x1y0y1MxMyZ0, (0,2,6)))
	
			IY1 = H(marginalizacao(Px0x1y0y1MxMyZ1, (3,))) + H(marginalizacao(Px0x1y0y1MxMyZ1, (1,6))) - H(marginalizacao(Px0x1y0y1MxMyZ1, (1,3,6)))
	
			HMxMy = H(marginalizacao(Px0x1y0y1MxMyZ1, (4,5)))
	
			Ix0x1 = H(marginalizacao(Px0x1y0y1MxMyZ0, (0,))) + H(marginalizacao(Px0x1y0y1MxMyZ0, (1,))) - H(marginalizacao(Px0x1y0y1MxMyZ0, (0,1)))
	
			Iy0y1 = H(marginalizacao(Px0x1y0y1MxMyZ0, (2,))) + H(marginalizacao(Px0x1y0y1MxMyZ0, (3,))) - H(marginalizacao(Px0x1y0y1MxMyZ0, (2,3)))
			

	#-------------TESTING INEQUALITY-------------------------------------------------------
	
			desig = IX0 + IX1 + IY0 + IY1 - HMxMy + Ix0x1 + Iy0y1
	
	
			if(desig>bound):
				print('U',epsilon, gamma, abs(desig - bound), desig)
				gammaU = gamma
				if(abs(desig - bound) > atol):#Up to the precision atol, devide the current gamma's range in half
					gamma = gammaU - (gammaU - gammaB)/2
				else:
					break
			else:
				print('B',epsilon, gamma, abs(desig - bound), desig)
				gammaB = gamma
				if(abs(desig - bound) > atol):#Up to the precision atol, devide the current gamma's range in half
					gamma = gammaB + (gammaU - gammaB)/2
				else:
					break
			aux = aux +1
	
		aux = 0
	
		gammaB = 0.0 #reset
		G[i] = gamma # Save the best gamma

	return (E, G)
		
	
	