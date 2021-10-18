from scipy.special import binom
from numpy import array, concatenate, abs, mean, std, sqrt, arange, linspace, polyfit, inf, where, poly1d, RankWarning
from numpy import nan, nanmean, nanstd
from numpy.random import choice

import matplotlib.pyplot as plt

from struct import pack,unpack

import os






if __name__ == "__main__":

	#Dimension of Hilbert space
	L , N, L_b = 10, 5, 1
	dim = int(binom(L,N))

	#Open the file
	f = open("../Data/phase-g.dat","rb")
	
	#Read in the data            
	seed = unpack('i',f.read(4))[0]
	g_step, g_stop = unpack('dd',f.read(16))
	eps_start, eps_stop, eps_step = unpack('ddd',f.read(24))
		

	#Construct the g & epsilon values			
	g_vals = arange(0,g_stop + g_step/2,g_step)
	g_len = len(g_vals)
	eps_vals = arange(eps_start,eps_stop+eps_step/2,eps_step)
	eps_len = len(eps_vals)

			
	#Read in the g_cs
	g_cs = unpack('d'*eps_len,f.read(8*eps_len))

	#Read in spectra
	open_spectrum = array(unpack('d'*dim, f.read(8*dim)))
	ring_spectrum = array(unpack('d'*dim, f.read(8*dim)))

	#Read in amplitudes
	amplitudes = array(unpack('d'*(L-L_b),f.read(8*(L-L_b))))

	print("Now plotting")
		for j,eps in enumerate(eps_vals):
			print("On eps = "+str(eps))
			for W_c in W_c_guesses[U]:
				plt.figure(figsize=(10,10))
				plt.subplot(211)
				for k,L in enumerate(L_vals):
					plt.errorbar((W_vals-W_c)*L,means[k,:,j]/L,fmt='-^',yerr=stds[k,:,j]/L,label=r'$L = '+str(L)+'$')

				plt.xlabel(r"$(W-W_c)L^{1/\nu}$",fontsize=12)
				plt.ylabel(r"$\xi/L$",fontsize=12)
				plt.legend()
				
				plt.subplot(212)
				for k,L in enumerate(L_vals):
					plt.errorbar((W_vals-W_c)*L,means[k,:,j]/L,fmt='-^',yerr=stds[k,:,j]/L,label=r'$L = '+str(L)+'$')

				plt.xlim([-20,20])
				plt.ylim([0,2])
				plt.xlabel(r"$(W-W_c)L^{1/\nu}$",fontsize=12)
				plt.ylabel(r"$\xi/L$",fontsize=12)
				plt.legend()
				

				plt.savefig("../Figures/U="+str(int(U))+"/new/scaling-eps/eps={:.2f}/W_c=".format(eps)+str(W_c)+".png")
				plt.close()

	