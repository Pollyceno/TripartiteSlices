import numpy as np
import ncpol2sdpa as nc
from recurso import Index322
from grupo import permuta

def quantumedge322(p, level):

	#p = outros(3)
	
	#p = recurso(1)
	#p = permuta(p)
	#print(len(p))
	
##################################################################################################
##################### DEFININDO CORRELATORES ####################################################
	
	Ax = np.zeros((2))
	By = np.zeros((2))
	Cz = np.zeros((2))
	
	AxBy = np.zeros((2,2))
	ByCz = np.zeros((2,2))
	AxCz = np.zeros((2,2))
	
	AxByCz = np.zeros((2,2,2))
	
	for x in range(2):
		for y in range(2):
			for z in range(2):
				ax = 0
				by = 0
				cz = 0
				axby = 0
				bycz = 0
				axcz = 0
				axbycz = 0
				for a in range(2):
					for b in range(2):
						for c in range(2):
							ax = ax + ((-1)**(a)*p[Index322((a,b,c,x,y,z))])
							by = by + ((-1)**(b)*p[Index322((a,b,c,x,y,z))])
							cz = cz + ((-1)**(c)*p[Index322((a,b,c,x,y,z))])
							axby = axby + ((-1)**(a+b)*p[Index322((a,b,c,x,y,z))])
							bycz = bycz + ((-1)**(b+c)*p[Index322((a,b,c,x,y,z))])
							axcz = axcz + ((-1)**(a+c)*p[Index322((a,b,c,x,y,z))])
							axbycz = axbycz + ((-1)**(a+b+c)*p[Index322((a,b,c,x,y,z))])
	
				Ax[x] = ax
				By[y] = by
				Cz[z] = cz
				AxBy[x][y] = axby
				ByCz[y][z] = bycz
				AxCz[x][z] = axcz
				AxByCz[x][y][z] = axbycz
	
	#print(Ax)
	#print(By)
	#print(Cz)
	#print(AxBy)
	#print(ByCz)
	#print(AxCz)
	#print(AxByCz)
	
	
##################################################################################################
##################### DEFININDO OBSERVAVEIS ####################################################
	
	n_A = 2 # Number of dichotomic observables of party A
	n_B = 2 # Number of dichotomic observables of party B
	n_C = 2 # Number of dichotomic observables of party C
	#P = nc.Probability([2,2],[2,2],[2,2])
	
	#level = 2 # Level of relaxation
	
	A = nc.generate_operators('A', n_A, hermitian=True)
	B = nc.generate_operators('B', n_B, hermitian=True)
	C = nc.generate_operators('C', n_C, hermitian=True)
	
	subs = {A[i] ** 2 :1 for i in range(n_A)}
	subs.update({B[i] ** 2 :1 for i in range(n_B)})
	subs.update({C[i] ** 2 :1 for i in range(n_C)})
	
	subs.update({A[i]*B[j]:B[j]*A[i] for i in range(n_A) for j in range(n_B)})
	subs.update({A[i]*C[j]:C[j]*A[i] for i in range(n_A) for j in range(n_C)})
	subs.update({B[i]*C[j]:C[j]*B[i] for i in range(n_B) for j in range(n_C)})
	'''
	subs.update({A[i]*B[j]*C[k]:AxByCz[i][j][k] for i in range(n_B) for j in range(n_B) for k in range(n_C)})
	
	subs.update({A[i]*B[j]:AxBy[i][j] for i in range(n_A) for j in range(n_B)})
	subs.update({A[i]*C[j]:AxCz[i][j] for i in range(n_A) for j in range(n_C)})
	subs.update({B[i]*C[j]:ByCz[i][j] for i in range(n_B) for j in range(n_C)})
	
	subs.update({A[i]:Ax[i] for i in range(n_A)})
	subs.update({B[i]:By[i] for i in range(n_B)})
	subs.update({C[i]:Cz[i] for i in range(n_C)})
	'''
######################################################################################################
##################### RESTRICOES DO COMPORTAMENTO ####################################################
	cons =[
		A[0]-Ax[0],
		A[1]-Ax[1],
		B[0]-By[0],
		B[1]-By[1],
		C[0]-Cz[0],
		C[1]-Cz[1],
		A[0]*B[0]-AxBy[0][0],
		A[0]*B[1]-AxBy[0][1],
		A[1]*B[0]-AxBy[1][0],
		A[1]*B[1]-AxBy[1][1],
		A[0]*C[0]-AxCz[0][0],
		A[0]*C[1]-AxCz[0][1],
		A[1]*C[0]-AxCz[1][0],
		A[1]*C[1]-AxCz[1][1],
		B[0]*C[0]-ByCz[0][0],
		B[0]*C[1]-ByCz[0][1],
		B[1]*C[0]-ByCz[1][0],
		B[1]*C[1]-ByCz[1][1],
		A[0]*B[0]*C[0]-AxByCz[0][0][0],
		A[0]*B[0]*C[1]-AxByCz[0][0][1],
		A[0]*B[1]*C[0]-AxByCz[0][1][0],
		A[0]*B[1]*C[1]-AxByCz[0][1][1],
		A[1]*B[0]*C[0]-AxByCz[1][0][0],
		A[1]*B[0]*C[1]-AxByCz[1][0][1],
		A[1]*B[1]*C[0]-AxByCz[1][1][0],
		A[1]*B[1]*C[1]-AxByCz[1][1][1]]

######################################################################################################
##################### MONOMIOS EXTRAS ####################################################

	Extra = []
	Extra.append([
	    A[i]*B[j]*C[k]
	    for i in range(n_A) 
	    for j in range(n_B) 
	    for k in range(n_C)
	    ])
	ExtraMonomials = [item for sublist in Extra for item in sublist]
	
##################################################################################################
##################### RESOLVENDO ################################################################
	
	sdp = nc.SdpRelaxation(A+B+C)
	if(level == 1):
		sdp.get_relaxation(level, objective=-0, substitutions=subs, momentequalities=cons, extramonomials=ExtraMonomials )
	if(level == 2):
		sdp.get_relaxation(level, objective=-0, substitutions=subs, momentequalities=cons)
	sdp.solve(solver='mosek')
	
	#print(sdp.primal,sdp.dual,sdp.status)
	
	return sdp.status
