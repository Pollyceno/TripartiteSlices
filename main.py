import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import matplotlib.colors as mcolors

from tripartite_criterion import *
from noisy_tripartite_criterion import *
from resource import *
from npa322 import NPA32
from multiple_copies_criterion import *





n=10

###########################################################################################################################
			# TRIPARTITE BOXES
###########################################################################################################################

p45 = box45(Eclass)

pd = deterministic()

pw = (1/8)*np.ones(64)



###########################################################################################################################
			# TRIPARTITE EDGE
###########################################################################################################################
print('#### One copy IC 322 ############')
(E,G)= IC_single_copy(n, p45, pd,  pw)

###########################################################################################################################
			# MULTIPLE 322
###########################################################################################################################
print('#### Multiple copies 322 ############')
(EM3, GM3) = multiple(n, p45, pd,  pw)


###########################################################################################################################
			# NPA Q2
###########################################################################################################################
print('#### NPA 322 Q2 ############')
(Eq, Gq) = NPA322(n, p45, pd,  pw, 2)

###########################################################################################################################
			# NOISY TRIPARTITE EDGE
###########################################################################################################################
print('#### Noisy 322 IC ############')
(En,Gn)= noisy_ic_edge(n, p45, pd,  pw)


###########################################################################################################################
			# PLOT
###########################################################################################################################


plt.figure()
plt.style.use('classic') #coloca ticks "in"
b = np.arange(0., 1, 0.01)

plt.plot(b, 1-b , label='NS',linewidth=1.5,c = 'k',linestyle='--')
plt.plot(E,G,label='IC single copy',linewidth=2.5, c='tab:red')
plt.plot(EM3,GM3,label='IC 322 multiple copies',linewidth=2.5, c='lightskyblue')
plt.plot(En,Gn,label='IC noisy',linewidth=1.5, c='tab:orange',linestyle='--')
plt.plot(EU,GU,label='Uffink',linewidth=1.5, c='tab:green', linestyle='--')
plt.plot(Eq,Gq,label='Q$_2$',linewidth=2.5, c='k')

plt.xlabel(r'$\varepsilon$', fontsize='30')
plt.xticks(fontsize=20)
plt.ylabel(r'$\gamma$', fontsize='30')
plt.yticks(fontsize=20)
plt.legend(fontsize=20,fancybox=True, shadow=True).get_frame().set_edgecolor('lightgray')
plt.ylim(0.0,1.01)
plt.xlim(-0.01,1)
ax=plt.gca()
ax.xaxis.set_major_locator(MultipleLocator(0.2))


fig = plt.gcf()
fig.set_size_inches(10, 8)
plt.tight_layout()
nome = 'test'+str(caixa2)+'.pdf'
plt.savefig(nome)
plt.close(fig)
print('#-------END---------#')
