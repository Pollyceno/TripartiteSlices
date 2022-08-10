import numpy as np

def multiple(n, p45, pd,  pw):
##################################################################################################################
# - Receives parameters' precision n, the three boxes that define the slice of interest p45, pd,  pw
# - Returns the parameters greatest gamma, for each epsilon, of the region pabcDxyz = gamma*p45 +epsilon*pd 
#   + (1-gamma-epsilon)*pw, for which the multiple copies inequality (EI)**2 + (EII)**2 <= 1 is not violated;
##################################################################################################################

	EM=np.arange(0., (1+1/n), 1/n)
	GM= 1 - EM
	
#-----------BISECTION PARAMETERS---------------------------------------------------------------
	iterac = 0
	gammaU = 1 
	gammaB = 0
	gamma = (gammaU - gammaB)/2
	atol=10**-5
	bound = 1+10**-8
#-----------VARING PARAMETERS------------------------------------------------------------------
	ç = n+1
	for i in range(0,ç):
		epsilon=(1/n)*i
		EM[i]=epsilon
		max_it = n
		aux = 0
		gammaU = 1-epsilon
		gamma = gammaU
		while(aux < max_it ):
			iterac = iterac +1
		
		#-----------SLICE----------------------------------------------------------------------

			pabcDxyz = gamma*p45 +epsilon*pd + (1-gamma-epsilon)*pw

		#-------------TESTING INEQUALITY-------------------------------------------------------
			desig = Computing_Inequality(3, pabcDxyz)

			if(desig>bound):
				gammaU = gamma
				if(abs(desig - bound) > atol): #Up to the precision atol, devide the current gamma's range in half
					gamma = gammaU - (gammaU - gammaB)/2
				else:
					break
			else:
				gammaB = gamma
				if(abs(desig - bound) > atol): #Up to the precision atol, devide the current gamma's range in half
					gamma = gammaB + (gammaU - gammaB)/2
				else:
					break
			aux = aux +1

		aux = 0

		gammaB = 0.0 #reset
		GM[i] = gamma # Save the best gamma
	
	return (EM, GM)


def Computing_Inequality(N_parts, p):
##################################################################################################################
# - Receive number of parts N_parts for the scenario and the behavior p (out/in = 0,1);
# - Return the inequality parameters P_I and P_II computed;
##################################################################################################################


	N = (2*N_parts) # Behavior's size
	P_I = 0 
	P_II = 0

	for i in range(0,2**N):

#-------Computing the bit-string I with (out/in)puts for p(a...|x...), respective to 'i'
		k = np.binary_repr(i)
		I = str(0)*(N-len(k))
		I = I+k

#-------Sum module 2 of all outputs a_i
		String_a = I[:int(len(I)/2)]
		String_a = np.array(list(map(int, list(String_a))))
		Sum_String_a = sum(String_a)%2

#-------Sum module 2 of all inputs x_i less the last
#		It make sense because we are conditionating the last x_n.
#		So, its a simple sum of the inputs.
		String_x = I[int(len(I)/2):]
		String_x = np.array(list(map(int, list(String_x))))
		String_xm1 = np.delete(String_x, len(String_x)-1) 
		Sum_String_xm1 = sum(String_xm1)%2

		#print(I[len(I)-1], String_a, String_xm1)
		#if(int(I[len(I)-1]) == 0):
		#	print(Sum_String_a, Sum_String_xm1)
#-------Computing success probability for the first and second boolean functions respectively	
		if(Sum_String_a == 0 and int(I[len(I)-1]) == 0 ):

			P_I = P_I + p[i]

		if(Sum_String_a == Sum_String_xm1 and int(I[len(I)-1]) == 1 ):
			P_II = P_II + p[i]
	

#-------Re-normalizing
	P_I = (1/2**(N_parts-1))*P_I
	P_II = (1/2**(N_parts-1))*P_II

	EI = 2*P_I - 1
	EII = 2*P_II - 1
	Ineq = ((EI)**2 + (EII)**2)

	return Ineq

