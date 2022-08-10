import math
import numpy as np
import itertools

def Index322(variables):
#######################################################################################################
# Receives the set o in/out-puts;
# Returns the respective index;
#######################################################################################################
	cardinalities = 6*[2]

	Possibilities = list(itertools.product(*[range(i) for i in cardinalities]))

	return Possibilities.index(variables)

def box45():
#######################################################################################################
# It computes the complete tripartite behavior for which (a+b+c) == ((x*z)+(y*z))
#######################################################################################################
	p = np.zeros(2**6, dtype='float')
	for k in range(0,(2**6)):
		k = np.binary_repr(k)
		I = np.zeros((6-len(k)), dtype = 'int')
		for l in range(0,len(k)):
			I = np.hstack((I,int(k[l])))
		a = I[0]
		b = I[1]
		c = I[2]
		x = I[3]
		y = I[4]
		z = I[5]
		if((a+b+c)%2 == ((x*z)+(y*z))%2):
			p[Index322((a,b,c,x,y,z))] = 1/4
		else:
			p[Index322((a,b,c,x,y,z))] = 0
	return p

def deterministic():
#######################################################################################################
# It computes the complete tripartite deterministic behavior for which a = 0, b = 0, and c = 0
#######################################################################################################
	p = np.zeros(2**6, dtype='float')
	for k in range(0,(2**6)):
		k = np.binary_repr(k)
		I = np.zeros((6-len(k)), dtype = 'int')
		for l in range(0,len(k)):
			I = np.hstack((I,int(k[l])))
		a = I[0]
		b = I[1]
		c = I[2]
		x = I[3]
		y = I[4]
		z = I[5]
		if(a==0 and b==0 and c ==0):
			p[Index322((a,b,c,x,y,z))] = 1
		else:
			p[Index322((a,b,c,x,y,z))] = 0
	return p