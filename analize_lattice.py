import cPickle as pkl
import numpy as np
import argparse,glob
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
import ROOT as r



def e_fit(x,t_c,c_exp):
	return pow(abs((x-t_c)/t_c),c_exp)


def log_e(x,t_c):
	return np.log(abs((x-t_c)/t_c))

def lin_e(x,exp):
	return exp*x


def bound(T,f):													# this function defines the initialization of
	bound_g = (-2,-0.75)											# the fit paramenters, returning the range
	bound_n = (-1.5,-0.5)											# of the correct quantity to fit
	bound_b = ([2*i_dim-2,0],[2*i_dim-1,0.75])
	tc = np.asarray([2.27,4.5,6.68,8.78])
	if i_dim == 2 or i_dim == 3: idx = (np.absolute(np.asarray(T)-tc[i_dim-2])).argmin()+1
	else: idx = (np.absolute(np.asarray(T)-tc[i_dim-2])).argmin()+2
	pt = 0
	if f == 'b':
		return [pt,idx,bound_b,tc[i_dim-2]]
	elif f == 'g':
		return [pt,idx,bound_g,tc[i_dim-2]]
	elif f == 'n':
		return [idx, bound_n,tc[i_dim-2]]


def magnetization(plot=False):											
	th_exp = [0.125,0.327,0.5,0.5]
	if str(args.accumulate()) == 'magnetization(True)' and root:						# first condition propmts ROOT interface
		graph = r.TGraphErrors(len(T),np.asarray(T),np.asarray(M),np.zeros(len(T)),np.asarray(sigma))	# second runs the fit with scipy
		graph.SetMarkerStyle(24)
		graph.Draw('ap')
		raw_input()
	else:
		b = bound(T,'b')
		if easy: p_b,cov_b = curve_fit(lambda x,bt:e_fit(x,b[3],bt),T[b[0]:b[1]],M[b[0]:b[1]])		# fit carried out with one parameter or two
		else: p_b,cov_b = curve_fit(e_fit,T[b[0]:b[1]],M[b[0]:b[1]],bounds=b[2],p0=(b[3],th_exp[i_dim-2]))
		if plot:
			print p_b,cov_b
			plt.scatter(T,M)
			if easy: plt.plot(np.linspace(T[b[0]],b[3],1000),e_fit(np.linspace(T[b[0]],b[3],1000),b[3],p_b[0]))
			else: plt.plot(np.linspace(T[b[0]],p_b[0],1000),e_fit(np.linspace(T[b[0]],p_b[0],1000),p_b[0],p_b[1]))
			plt.plot(np.linspace(T[0],b[3],1000),e_fit(np.linspace(T[0],b[3],1000),b[3],th_exp[i_dim-2]))
			plt.show()
		return p_b[0]


def susceptibility(T,M,plot=False):										# susceptibility function returning the
	th_exp = [-1.75,-1.237,-1,-1]
	S = np.absolute(np.divide(np.gradient(M),np.gradient(T)))
	if root and plot:
		graph = r.TGraphErrors(len(T),np.asarray(T),np.asarray(S),np.zeros(len(T)),math.sqrt(2)*np.asarray(sigma))
		graph.SetMarkerStyle(24)
		graph.Draw('ap')
		raw_input()
	elif plot and not root:
		tc = magnetization()
		g = bound(T,'g')
		p_g,cov_g = curve_fit(lin_e,log_e(T[0:g[1]],tc),np.log(S[0:g[1]]),bounds=g[2],p0=th_exp[i_dim-2])
		if plot:
			print p_g
			plt.scatter(T,S)
			plt.plot(np.linspace(T[0],T[-1],10),e_fit(np.linspace(T[0],T[-1],10),g[3],p_g[0]))
			plt.plot(np.linspace(T[0],T[-1],10),e_fit(np.linspace(T[0],T[-1],10),g[3],th_exp[i_dim-2]))
			plt.show()
	return S


def cfunct(in_spin,plot=False):											# correlation function (see the paper for reference)
	def expo(x,l): return np.exp(-x/l)
	N = len(in_spin[0][1])
	length = pow(N,1./float(i_dim))
	rg_length = np.arange(length/2)
	dim =  np.full(i_dim,int(round(length)),dtype=int)
	glob_corr_funct = []
	l_vect = []
	for t,e in enumerate(T):
		corr_funct = []
		c_spin = np.asarray(in_spin[t][1])
		c_spin = np.reshape(c_spin,dim)
		for l in rg_length:
			sum_c_funct = 0
			for i in xrange(i_dim):
       		 		sum_c_funct +=  np.sum(np.roll(c_spin,int(l),axis=i)*c_spin+np.roll(c_spin,int(-l),axis=i)*c_spin)
			c_funct_val = abs(sum_c_funct/(2.*i_dim*N))-pow(np.mean(c_spin),2)
			corr_funct.append(c_funct_val)
		p_l,cov_l = curve_fit(expo,rg_length,corr_funct)
		glob_corr_funct.append([e,corr_funct])
		l_vect.append(float(p_l))

	if plot:												# plot the correlation function if the flag is True
		color = plt.cm.rainbow(np.linspace(0,1,len(glob_corr_funct)))
		for i,e in enumerate(glob_corr_funct):
			plt.plot(rg_length,e[1],c=color[i],label=str(e[0]))
		plt.legend(ncol=3,prop={'size':8})
		plt.xlabel('distance')
		plt.ylabel('correlation function')
		plt.show()
	return T,l_vect


def ext_cfunct(T,i_dim,in_spin,plot=False):									# utlity function (not called in this script)
	def expo(x,l): return np.exp(-x/l)
	N = len(in_spin[0][1])
	length = pow(N,1./float(i_dim))
	rg_length = np.arange(length/2)
	dim =  np.full(i_dim,int(round(length)),dtype=int)
	glob_corr_funct = []
	l_vect = []
	for t,e in enumerate(T):
		corr_funct = []
		c_spin = np.asarray(in_spin[t][1])
		c_spin = np.reshape(c_spin,dim)
		for l in rg_length:
			sum_c_funct = 0
			for i in xrange(i_dim):
       		 		sum_c_funct +=  np.sum(np.roll(c_spin,int(l),axis=i)*c_spin+np.roll(c_spin,int(-l),axis=i)*c_spin)
			c_funct_val = abs(sum_c_funct/(2.*i_dim*N))-pow(np.mean(c_spin),2)
			corr_funct.append(c_funct_val)
		p_l,cov_l = curve_fit(expo,rg_length,corr_funct)
		glob_corr_funct.append([e,corr_funct])
		l_vect.append(float(p_l))
	return l_vect


def c_len(in_spin, plot=False):											# plot and fit correlation length 
	T,L = cfunct(in_spin)
	th_exp = [-1.75,-1.237,-1,-1]
	if root and plot:
		graph = r.TGraphErrors(len(T),np.asarray(T),np.asarray(L),np.zeros(len(T)),np.asarray(sigma))
		graph.SetMarkerStyle(24)
		graph.Draw('ap')
		raw_input()
	elif plot and not root:
		tc = magnetization()
		n = bound(T,'n')
		if easy: p_n,cov_n = curve_fit(lin_e,log_e(T[0:n[0]],n[2]),np.log(L[0:n[0]]),bounds=n[1])
		else: p_n,cov_n = curve_fit(lin_e,log_e(T[0:n[0]],tc),np.log(L[0:n[0]]),bounds=n[1])
		print p_n,cov_n
		plt.scatter(T,L)
		if easy: plt.plot(np.linspace(T[0],T[n[0]],10),e_fit(np.linspace(T[0],T[n[0]],10),n[2],p_n[0]))
		else: plt.plot(np.linspace(T[0],T[n[0]],10),e_fit(np.linspace(T[0],T[n[0]],10),tc,p_n[0]))
		plt.show()

def main(args):
	eval(args.accumulate())											# the quantity analyzes is the one parsed from line of command

root = False

if __name__ == '__main__':
	parser = argparse.ArgumentParser()									# script run as python analyze_lattice.py -f n_dim
	parser.add_argument('n_dim',help='n-dim lattice')							# where -f is the quantity to analyze
	parser.add_argument('-e',help='returns only the critical exponent, leave empty to return the critical temperature too',action='store_const',dest='easy',const='a',default=None)
	parser.add_argument('-r',help='prompts ROOT fit panel',action='store_const',dest='root',const='b',default=None)
	parser.add_argument('-m',help='plots the magnetization and prints beta (and t_c if -e active)',dest='accumulate',action='store_const',const=lambda:'magnetization(True)')
	parser.add_argument('-s',help='plots the magnetic susceptibility and prints gamma',dest='accumulate',action='store_const',const=lambda:'susceptibility(T,M,True)')
	parser.add_argument('-cf',help='plots the correlation function',dest='accumulate',action='store_const',const=lambda:'cfunct(spin,True)')
	parser.add_argument('-l',help='plots the correlation length and prints nu',dest='accumulate',action='store_const',const=lambda:'c_len(spin,True)')
	args = parser.parse_args()
	i_dim = int(args.n_dim)

	with open('pickle/pickle'+str(i_dim)+'d_multip.p','rb') as rfp:						# the file pickle file produced in cluster_ising.py
		in_spin = list(pkl.load(rfp))									# is opened and initialized

	T = []
	M = []
	sigma = []
	spin = []
	for e in in_spin:
		T.append(e[0])
		M.append(e[1]+e[2])
		sigma.append(e[2])
		spin.append([e[0],e[3]])
	easy = bool(args.easy)
	root = bool(args.root)
	main(args)

