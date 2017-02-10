from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import collections
import math


def get_time_to_wall(photon_start_pos,length_of_scint,width_of_scint):
	# We can estimate the scintillator as cylindrical
	# so we have a symmetry and don't have to worry about
	# 3 dimensions, only need to check x and y.

	# generate a bunch of unit vectors
	angle_step = 0.005
	theta = np.arange(0,np.pi+angle_step,angle_step)

	vecs = np.vstack([np.sin(theta),np.cos(theta)]).T
	mags = np.linalg.norm(vecs, axis=-1)
	vecs = vecs / mags[..., np.newaxis]

	# For each vector, find how long it will go a distance x
	start = photon_start_pos
	x = length_of_scint - start
	c = 29.9792 # in cm/ns
	n = 1.5 # index of refraction of bar

	# i.e. how long until this is x
	time_to_hit = np.zeros(len(vecs[:,0]))
	for time in np.arange(14,500,0.01):
		pos = vecs[:,0]*c/n*time
		for ind in range(len(pos)):
			if pos[ind]>=x:
				if time_to_hit[ind]==0:
					time_to_hit[ind] = time

	# based on how long it takes to get to the end
	# of the bar, find out where in the radial direction
	# that photon will be -- going to be used as acceptance
	y_proj = np.cos(np.asarray(theta))

	y_dist =  c/n*time_to_hit*y_proj
	x_dist = c/n*time_to_hit*np.sin(np.asarray(theta))
	number_of_reflections = np.floor(y_dist/(2*width_of_scint))

	y_dist = y_dist - number_of_reflections*2*width_of_scint
	y_dist -= width_of_scint

	return np.column_stack([theta,time_to_hit,y_dist])


def read_data_set(length_of_scint,width_of_scint):
	# Read the array from disk
	new_data = np.loadtxt('scintLength_'+str(length_of_scint)+\
		'_scintWidth_'+str(width_of_scint)+\
		'.txt')
	
	new_data = new_data.reshape((4,630,3))

	#.                | is for x slice.     | is for entry in the [theta,time,y] row
	plt.plot(new_data[0][:,0],new_data[0][:,2])
	plt.show()
	return

def compile_lookup(length_of_scint,width_of_scint):
	x_samples = np.arange(0,length_of_scint+1,1)
	data = []
	# open look up table to fill
	with file('scintLength_'+str(length_of_scint)+\
		'_scintWidth_'+str(width_of_scint)+\
		'.txt','w') as outfile:
		outfile.write('#  Angle | Time to Wall | Radial Pos on Wall |\n')

		# sample in start position
		for photo_start in x_samples:
			outfile.write('# New slice, x = '+str(photo_start)+'\n')
		
			# get data for that start position
			data = get_time_to_wall(photo_start,length_of_scint,width_of_scint)
			np.savetxt(outfile,data)

	return


#get_time_to_wall(0,300,5)
#test_data_set(0,300,5)
compile_lookup(300,5)