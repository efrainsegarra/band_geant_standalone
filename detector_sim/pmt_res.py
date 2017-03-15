from __future__ import division
import matplotlib.pyplot as plt
import numpy as np

# Relevant reference papers:
# https://arxiv.org/pdf/1610.05667.pdf
# http://citeseerx.ist.psu.edu/viewdoc/download;jsessionid=308EB1E49D97380D43DB8245785F4964?doi=10.1.1.659.2431&rep=rep1&type=pdf

# Fudge factor -- will have to extract from data. In sigma_t(x), the
# A(x) factors should really be N(x), however I don't have a conversion
# from amplitude to the number of photoelectrons in a single PMT based
# on position. They should be proportional, so this factor just scales
# it to the paper I'm referencing.
def N(x):
	return (1600.)*4
def other_fudge(x):
	return 0.08

# Amplitude of signal in PMT based on position
def A(x):
	A0 = 0.07959
	A1 = 0.1015
	A2 = 4.582e-6
	return A0*np.exp(-x/8.86) + A1*np.exp(-x/172.6) + A2*np.exp(x/36.96)
 
# Time resolution of single PMT based on position
def sigma_t(x):
	sigma_pl = 1.61 #ns
	sigma_l = 0.89 #ns
	sigma_el = 0.0 #ns

	return np.sqrt(  sigma_pl**2/A(x) + (sigma_l*x)**2/A(x) + sigma_el**2  )/N(x) + other_fudge(x)


# Effects on time resolution included above:

# (1) The amount of time for the photoelectron to traverse the PMT
# to the first dynode and for the electron shower to traverse
# the dynode chain can both vary. This is referred to as the 'PMT' jitter

# (2) Emission rise/decay time of scintillator

# --> (1) and (2) depend on # photoelectrons and are combined to sigma_pl

# (3) Light pulse diffusion due to light transmission through scintillator. 
# this is defined as sigma_l, and also depends on # of photoelectrons, as
# well as grows in x.

# (4) Electronic readout / setup noise, referred to as sigma_el


# Combine the 2 PMT resolutions by error propagation
def t_weight(x):
	#num = t1 / sigma_t(x)**2 + t2 / sigma_t(300-x)**2
	summ = 1. / sigma_t(x)**2 + 1. / sigma_t(300-x)**2
	
	return np.sqrt(1/summ)

def main():
	x = np.arange(0.,300.1,0.1)
	plt.plot(x,sigma_t(x),label='PMT1')
	plt.plot(x,sigma_t(300-x),label='PMT2')
	plt.plot(x,t_weight(x),label='Weighted Avg')
	plt.legend(loc='best',numpoints=1)
	plt.xlabel('Distance from PMT1 [cm]')
	plt.ylabel('Time resolution [ns]')
	plt.title('PMT Time Resolution by Position')
	plt.text(85,0.28,"$\\sigma_T=\sqrt{\\frac{\\sigma^2_{p}}{N(x)} + \\frac{(\\sigma_l x)^2}{N(x)} + \\sigma_{el}^2 } $",fontsize=20)
	plt.grid(True)
	plt.show()

main()

