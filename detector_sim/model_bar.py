from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import collections

# Assume 2D bar for now (no z-direction thickness)


# Generates random propagation vectors for the photons
# in the scintillator bar
def generate_propagation(dims, number):
    vecs = np.random.normal(size=(number,dims))
    mags = np.linalg.norm(vecs, axis=-1)

    return vecs / mags[..., np.newaxis]

def check_boundary(positions,length,width,thickness):
	# input is a list of position vectors for each photon

	# Wish to check when element 1 or 2 is >= width, thickness (respectively)
	check_x = np.absolute(positions[:,0])
	check_y = np.absolute(positions[:,1])
	check_z = np.absolute(positions[:,2])

	# returns list of bools (true/false) so that
	# if an index is true in check_y, that photon
	# must be flipped in y direction. Also checks
	# if the photon at end of bar
	return check_x >= length/2, check_y >= width/2, check_z >= thickness/2

def check_pmt(positions,length,diameter):
	# get x position of all photons
	check_x = np.absolute(positions[:,0])

	# get radial distance from x-axis
	check_yz = positions[:,1]**2 + positions[:,2]**2

	# check if radial position is inside aperture
	check_yz = check_yz <= diameter/2
	# check if also photon is at the end of the bar
	check_x =  check_x >= length/2

	return check_yz & check_x

def run_along_x(x_pos):
	# quicker method would be to use a look up table
	# but I think this is good enough for the qualitative
	# answer I want


	# note i only want a single PMT number so really
	# x just goes from 0 to 300 and I lose half of my guys

	# Number of initial photons created
	N0 = 1000.

	# Number of dimensions
	d = 3

	# Create randomized propagation direction vectors
	prop_vecs = generate_propagation(d, N0)
	image_vecs = prop_vecs[:,0:2]


	# Build bar
	l = 300. # length in cm (x)
	w = 5. # width in cm (y)
	t = 5. # thickness in cm (z)

	# Build PMT aperture
	diam = 2.5 # diameter in cm

	c = 29.9792 # in cm/ns
	n = 1.5 # index of refraction of bar

	step_size = 0.05
	# the range of how far this program goes
	# depends on the window size of the PMT accept
	# signal
	t_steps = np.arange(0.,100,step_size)
	
	init_pos = np.zeros((N0,2))
	init_pos = np.c_[np.zeros((N0,1))+x_pos ,init_pos]


	path_x=[]
	path_y=[]
	number_photons= []
	for step in t_steps:
		path_x.append(init_pos[:,0])
		path_y.append(init_pos[:,1])

		# Propagate position of photon
		pos = init_pos + c*step_size*prop_vecs/n
		init_pos = pos

		# Check if photon at edge in y/z to reflect
		# and check if photon has reached end of scintillator
		is_at_boundary = check_boundary(pos,l,w,t)
		flip_x, flip_y,flip_z = is_at_boundary[0], is_at_boundary[1], is_at_boundary[2]

		# if photon at end of scintillator in x, kill
		# propagation vector
		for ind in range(len(flip_x)):
			if flip_x[ind] == True:
				prop_vecs[ind] = np.zeros(d)

		# if photon is at a y/z wall, reflect the propagation
		# vector
		flip_y = [1 if x==0 else x for x in flip_y*-1] 
		flip_z = [1 if x==0 else x for x in flip_z*-1]
		reflection = np.column_stack((np.ones(N0),flip_y,flip_z))

		# reflect any propagation, or kill any photons, if needed
		prop_vecs = reflection*prop_vecs
		
		# Now check how many are inside PMT aperture
		N = np.sum(check_pmt(pos,l,diam))

		number_photons.append(N)


	#plt.figure(1)
	# plot # photons in PMT as function of time
	#plt.plot(t_steps,number_photons)
	
	# plot sample paths
	#plt.figure(2)
	#plt.plot(path_x,path_y)

	# plots images of unit vectors that were randomized
	# figure, axis = plt.subplots()
	# vectors = np.insert(image_vecs[:, np.newaxis], 0, 0, axis=1)
	# axis.add_collection(collections.LineCollection(vectors))
	# axis.axis((-1, 1, -1, 1))
	#plt.show()
	return number_photons[-1]

def main():
	x_samples = np.arange(-140,141,10)
	photons = []
	for x in x_samples:
		print x
		photons.append(run_along_x(x))
	plt.plot(x_samples,photons)
	plt.show()

main()

