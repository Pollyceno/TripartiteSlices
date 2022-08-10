import numpy as np

def marginalization(P, marg):
##################################################################################################################
# - Receives a joint probability P and a set o variables to be marginalized;
# - Returns marginal joint probability distribution;
##################################################################################################################
	n =P.ndim
	x=np.arange(n)
	
	elim = np.delete(x, marg)
	elim = np.sort(elim)[::-1].astype(int)
	
	for i in elim:
		P = P.sum(axis=i,dtype='float')
	
	return P