from scipy.special import binom
from numpy import array, abs, mean, nanmean, std, var, count_nonzero, sqrt, log, arange, inf, nan, round
from numpy.linalg import norm

from numpy.random import choice

from struct import pack,unpack

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

def my_std(arr):
	return sqrt(my_mean(arr**2) - my_mean(arr)**2)



def bootstrap_error(arr,samples):

	data = choice(arr,size=(samples,len(arr)))
	means = nanmean(data,axis=1)

	return std(means)



if __name__ == "__main__":

	L_vals = [(12,3)]
	U_vals = [0,1,3,5]
	W_starts = [1.0,1.0,3.0,4.0]
	W_stops = {12:[7.0,7.0,9.0,10.0]}
	W_step = 0.1

	eps_start = .05
	eps_end = .95
	eps_step = .01
	eps_vals = arange(eps_start,eps_end,eps_step)
	eps_num = len(eps_vals)

	eps_inds = [1,1,1,1]

	beta = 1

	for L,L_b in L_vals:
		L_loc = L-L_b

		dim = int(binom(L,L/2))

		for U,W_start,W_stop,eps_ind  in zip(U_vals,W_starts,W_stops[L],eps_inds):
			print('Now accessing L = '+str(L)+', U = '+str(U))
			W_num = int(round((W_stop-W_start)/W_step)) + 1 
			W_vals = [W_start + W_step*i for i in range(W_num)]



			counts_plus = []
			counts_minus = []

			for W in W_vals:
				f = open('../Data/U='+str(int(U))+'/uni'+'{:.6f}'.format(W)+'-eps.dat','rb')

				trials = unpack('i',f.read(4))[0]

				g_cs = []
				energies = []
				max_amps = []	
				counts_plus_this_W = []
				counts_minus_this_W = []	

				for i in range(trials):
					g_vals = array(unpack('d'*dim,f.read(8*dim)))
					energies = array(unpack('d'*dim,f.read(8*dim)))
					amplitudes = array(unpack('d'*L_loc,f.read(8*L_loc)))


					E_min , E_max = min(energies) , max(energies)
					delta_E = E_max - E_min
					E_norm = (energies - E_min)/delta_E
					indices = [abs(energies - eps).argsort()[0:eps_ind] for eps in eps_vals]

					g_cs.append([my_mean(g_vals[inds]) for inds in indices]) 
					max_amps.append([max([(2*log(1/amplitudes[i]) - L*log(2) + log((eps*delta_E + E_min)*beta))/(2*(i+1)) for i in range(L_loc)]) for eps in eps_vals])
				
				g_cs = array(g_cs)
				max_amps = array(max_amps)

				for i in range(len(eps_vals)):
					counts_plus_this_W.append(count_nonzero(g_cs[:,i] >= max_amps[:,i])/trials)
					counts_minus_this_W.append(count_nonzero(-g_cs[:,i] >= max_amps[:,i])/trials)

				counts_plus.append(counts_plus_this_W)
				counts_minus.append(counts_minus_this_W)

				f.close()


			
			f = open('../Data/U='+str(U)+'/g-counts.dat','wb')

			f.write(pack('i',W_num))
			f.write(pack('d'*W_num,*W_vals))
			f.write(pack('i',eps_num))
			f.write(pack('d'*eps_num,*eps_vals))

			for i in range(W_num):
				f.write(pack('d'*eps_num,*counts_plus[i]))
				f.write(pack('d'*eps_num,*counts_minus[i]))
			
			f.close()
