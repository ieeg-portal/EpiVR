import os
import sys
import re

import numpy as np
from numpy import linalg as LA
from matplotlib.mlab import PCA
import scipy
import scipy.spatial

def centeroidnp(arr):
    length = arr.shape[0]
    sum_x = np.sum(arr[:, 0])
    sum_y = np.sum(arr[:, 1])
    sum_z = np.sum(arr[:, 2])
    return sum_x/length, sum_y/length, sum_z/length



csv_fn = '/Users/lkini/Documents/LittLab/aim1/results/HUP065/aim1/HUP065_T1_19991230_soz_coord.csv'

# open CSV file
lines = open(csv_fn,'r').readlines()

# Get unque grids/strips
elex = []
for line in lines:
	label = line.split(',')[3]
	m = re.compile(r'''[A-Za-z]+''').match(label)
	if m.group() not in elex:
		elex.append(m.group())

elex_coord = {}
# For each grid strip
for ele in elex:
	elex_coord[ele] = []
	for line in lines:
		x = np.float(line.split(',')[0])
		y = np.float(line.split(',')[1])
		z = np.float(line.split(',')[2])
		label = line.split(',')[3]
		m = re.compile(r'''[A-Za-z]+''').match(label)
		if m.group() == ele:
			elex_coord[ele].append(np.array([x,y,z]))
for ele in elex:
	elex_coord[ele] = np.array(elex_coord[ele])

for ele in elex:
	hull = scipy.spatial.ConvexHull(elex_coord[ele])
	vx = elex_coord[ele][hull.vertices]
	c_vx = centeroidnp(vx)
	
	x_min = min(vx[:,0])
	x_max = max(vx[:,0])
	y_min = min(vx[:,1])
	y_max = max(vx[:,1])
	z_min = min(vx[:,2])
	z_max = max(vx[:,2])

	verts = [(x_max,y_max,z_min),(x_max,y_min,z_min),(x_min,y_max,z_min),(x_min,y_max,z_min),(x_max,y_max,z_max),(x_max,y_min,z_max),(x_min,y_min,z_max),(x_min,y_max,z_max)]
	
	faces = [(0, 1, 2, 3),(4, 7, 6, 5),(0, 4, 5, 1),(1, 5, 6, 2),(2, 6, 7, 3),(4, 0, 3, 7)]
	print verts
	blah


verts = [(59.565979029600001, 50.282608038799999, -22.208232879600001), (59.565979029600001, 43.3451080432, -22.208232879600001), (20.284729028400001, 50.282608038799999, -22.208232879600001), (20.284729028400001, 50.282608038799999, -22.208232879600001), (59.565979029600001, 50.282608038799999, -7.9082328796399999), (59.565979029600001, 43.3451080432, -7.9082328796399999), (20.284729028400001, 43.3451080432, -7.9082328796399999), (20.284729028400001, 50.282608038799999, -7.9082328796399999)]
faces = [(0, 1, 2, 3),(4, 7, 6, 5),(0, 4, 5, 1),(1, 5, 6, 2),(2, 6, 7, 3),(4, 0, 3, 7)]
import bpy
mesh_data = bpy.data.meshes.new("cube_mesh_data")
mesh_data.from_pydata(verts, [], faces)
mesh_data.update() # (calc_edges=True) not needed here

cube_object = bpy.data.objects.new("Cube_Object", mesh_data)

scene = bpy.context.scene  
scene.objects.link(cube_object)  
cube_object.select = True
