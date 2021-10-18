from matplotlib import rc
import matplotlib.pyplot as plt

from scipy.special import binom
from numpy import array, abs, mean, nanmean, std, var, sqrt, arange, inf, nan, meshgrid, around
from struct import unpack

rc('text',usetex=True)



def my_mean(arr):
	sum = 0
	num = 0

	for val in arr:
		if val:
			sum += val
			num += 1

	if num:
		return sum/num
	else:
		return nan




if __name__ == "__main__":

	L_vals = [12]
	U_vals = [0,1,3,5]
	


	for L in L_vals:
		dim = int(binom(L,L/2))

		for U in U_vals:
			print('Now accessing L = '+str(L)+', U = '+str(U))

			f = open('../Data/U='+str(U)+'/g-counts.dat','rb')

			#Extract the values of W and epsilon
			W_num = unpack('i',f.read(4))[0]
			W_vals = array(unpack('d'*W_num,f.read(8*W_num)))
			eps_num = unpack('i',f.read(4))[0]
			eps_vals = array(unpack('d'*eps_num,f.read(8*eps_num)))

			plus = []
			minus = []

			for W in  W_vals:
				plus.append(list(unpack('d'*eps_num,f.read(8*eps_num))))
				minus.append(list(unpack('d'*eps_num,f.read(8*eps_num))))


			#### Make a density plot of the results
			plus = array(plus).T
			minus = array(minus).T

			print(plus[:,45])
			print(minus[:,45])


			fig_plus = plt.figure(figsize=(12,12))
			fig_minus = plt.figure(figsize=(12,12))
			ax_plus = fig_plus.gca()
			ax_minus = fig_minus.gca()

			W,eps = meshgrid(W_vals , eps_vals)
			#W_inds = [abs(W[0,:] - W_c).argmin() for W_c in W_cs]
			#eps_inds = [abs(eps[:,0] - epsilon).argmin() for epsilon in epsilons]

			im_plus = ax_plus.pcolormesh(W,eps,plus,shading='nearest',cmap='hot',vmin=0.0,vmax=1.0)
			cbar_plus = fig_plus.colorbar(im_plus,ax=ax_plus)
			cbar_plus.set_label(label=r'$x^+$',size=15)

			im_minus = ax_minus.pcolormesh(W,eps,minus,shading='nearest',cmap='hot',vmin=0.0,vmax=1.0)
			cbar_minus = fig_plus.colorbar(im_minus,ax=ax_minus)
			cbar_minus.set_label(label=r'$x^-$',size=15)

			#ax.scatter(W[0,W_inds],eps[eps_inds,0],color='b',marker=">")


			ax_plus.set_xlabel(r'Disorder Strength $W$',fontsize=15)
			ax_plus.set_ylabel(r'Energy density $\epsilon$',fontsize=15)
			ax_plus.tick_params(which='both',labelsize=12)
			fig_plus.savefig('../Figures/Phase-Portraits/U='+str(U)+'-plus.png')
		
			ax_minus.set_xlabel(r'Disorder Strength $W$',fontsize=15)
			ax_minus.set_ylabel(r'Energy density $\epsilon$',fontsize=15)
			ax_minus.tick_params(which='both',labelsize=12)
			fig_minus.savefig('../Figures/Phase-Portraits/U='+str(U)+'-minus.png')
			
			plt.close(fig_plus)
			plt.close(fig_minus)