#!/usr/bin/env python
# Authors: Sevan Adourian and Dan Frost
# sevan.adourian@berkeley.edu

import json
import numpy as np
import os
import pickle
import pyspl
import scipy.sparse as sparse
import time
import sys
import yaml
import pdb
from numba import jit
from numba.typed import Dict

from scipy.sparse.linalg import lsqr
# from scipy.interpolate import SmoothSphereBivariateSpline
from scipy.interpolate import RectSphereBivariateSpline

from Model1D import Model1D, parameter_column_map

#######
def setup(conf):
    # Read in grid and spline data and populate to a structure
    # read in grid data ...
    print('reading spherical spline knot locations and levels')
    tmp = np.loadtxt(conf['sspl_knots_file'], skiprows=1)
    conf['sspl_lons'] = tmp[:,0]
    conf['sspl_lats'] = tmp[:,1]
    conf['sspl_levels'] = tmp[:,2].astype(int)

    print('reading b-spline knot radii ...')
    conf['bspl_radii'] = np.loadtxt(conf['bspl_knots_file'], skiprows=1)
    print('reading lateral sampling grid file: <lon> <lat>')
    tmp = np.loadtxt(conf['lateral_grid_file'])
    conf['grid_lons'] = np.deg2rad(tmp[:,0])
    conf['grid_lats'] = np.deg2rad(tmp[:,1])

    print('reading radial sampling grid file: <radius>')
    conf['grid_r'] = np.loadtxt(conf['radial_grid_file'])

    # read in target model samples
    print('reading model samples file: <sample>')
    conf['values'] = np.loadtxt(conf['model_aggregated_file']).ravel()

    return conf


######################################
def grid2spline(conf, rad_index):
    ## Interpolating the dv values from gridded to spline data
    nb_grid_lat_lon = len(conf['grid_lons'])
    nb_spln_lat_lon = len(conf['sspl_lons'])
    
    colats = (np.pi/2)-conf['grid_lats'][:]
    lons   = conf['grid_lons'][:]

    conf['colats_spln'] = (np.pi/2)-np.deg2rad(conf['sspl_lats'])
    conf['lons_spln']   = np.deg2rad(conf['sspl_lons'])

    values_slice = conf['values'][rad_index*nb_grid_lat_lon:(rad_index+1)*nb_grid_lat_lon]
    colats_uniq = np.unique(colats)
    lons_uniq   = np.unique(lons)
    colats_uniq[0]  = colats_uniq[0]  + 1e-6
    colats_uniq[-1] = colats_uniq[-1] - 1e-6
    lons_uniq[0]    = lons_uniq[0]    + 1e-6
    values_slice = np.reshape(values_slice, (len(colats_uniq), len(lons_uniq)))
    # pdb.set_trace()
    conf['interpolator'] = RectSphereBivariateSpline(colats_uniq, lons_uniq, values_slice)
  
    return conf
    

# const

# earth radius
rn = 6371.0

########
# config
####################
## model construction
## CREATION OF THE INTERPOLATOR ON THE SPHERE
for conf_file in sys.argv[1:]:
    print('using configuration loaded from: %s' % (conf_file))
    with open(conf_file) as f:
        conf = yaml.safe_load(f)
        conf = setup(conf)

# compute interp
sspl = pyspl.SphericalSplines(conf['sspl_lons'], conf['sspl_lats'], conf['sspl_levels'])
H = sspl.evaluate(conf['sspl_lons'], conf['sspl_lats'])
# kn_lon, kn_lat, lev = sspl.get_knots_and_levels()

print('Have %i x %i spherical spline interpolant' % H.shape)

# compute radial sampling and interpolant
sampling_knots = conf['bspl_radii'].copy()

bspl = pyspl.CubicBSplines(conf['bspl_radii'])
V = bspl.evaluate(sampling_knots)
print('Have %i x %i radial b-spline interpolant' % V.shape)

# compute sampling ...
t0 = time.time()
b_gc = 0.0
b_gs = 0.0

pert_gc = []
pert_gs = []

## EVALUATION OF THE GRIDDED MODEL
# Looping through the bsplines
for no,r in enumerate(sampling_knots):
    pert_r_gc = np.zeros(H.shape[0])
    pert_r_gs = np.zeros(H.shape[0])

    print("... computing interpolant for bspline knot %i/%i ..." % (no+1, len(sampling_knots)))
    # Compute the interpolator for this grid level
    conf = grid2spline(conf, no)
    # pdb.set_trace()
    colats_sspl = (np.pi/2) - np.deg2rad(conf['sspl_lats'])
    lons_sspl   = np.deg2rad(conf['sspl_lons'])
    colats_sspl[np.where(colats_sspl == 0)] = 1e-6
    colats_sspl[np.where(colats_sspl == np.pi)] = np.pi - 1e-6
    lons_sspl[np.where(lons_sspl == 0)] = 1e-6

    pert_r_gc = conf['interpolator'].ev(colats_sspl, lons_sspl)
    
    # save anomaly at this radius
    pert_gc.append(pert_r_gc)
    pert_gs.append(pert_r_gs)

b_gc += np.hstack(pert_gc)
b_gs += np.hstack(pert_gs)

t1 = time.time()
print('Built model sampling in %.1f s' % (t1 - t0))

# compute product for fitting
t0 = time.time()
A = sparse.kron(V, H)
t1 = time.time()
print('Computed Kronecker product in %.1f s' % (t1 - t0))

# solve
t0 = time.time()
x_gc = lsqr(A, b_gc, show=True, iter_lim=2000)
x_gs = lsqr(A, b_gs, show=True, iter_lim=2000)

# Extract Mean 1D from A3d Model
print('Extracting the 1D Arithmetic mean from the A3d Model')
weights = np.cos(np.radians(conf['sspl_lats']))
weights = weights.reshape((weights.size,1))

# load the reference model, derive reference parameters
ref = Model1D(conf['path_to_ref_model'])
ref.load_from_file()

vsv = ref.get_values(1000 * sampling_knots, parameter='vsv')
vsh = ref.get_values(1000 * sampling_knots, parameter='vsh')
vs0 = 0.001 * np.sqrt((2 * vsv ** 2 + vsh ** 2) / 3)

c_gc = x_gc[0].reshape((V.shape[1], H.shape[1]))
c_gs = x_gs[0].reshape((V.shape[1], H.shape[1]))

H = sspl.evaluate(conf['sspl_lons'], conf['sspl_lats'])

gc = H * (V * c_gc).T
gs = H * (V * c_gs).T

vsv2 = (1.0 + np.sqrt(gc**2 + gs**2)) * vsv**2
mean = np.sqrt((weights * vsv2).sum(axis = 0) / weights.sum())

# PRINTING A3d

t1 = time.time()
print('Computed model coefs in %.1f s' % (t1 - t0))

# write out
model_fname = conf['output_file']
# num_xi = 642
num_s  = H.shape[0]

anis = np.zeros([26,642])
with open(model_fname, 'w') as f:
    f.write('1\n') # Number of parameters
    f.write('%i %i S\n'  % (num_s,  V.shape[1])) # Size of the grid
    # f.write('%i %i X\n'  % (num_xi, V.shape[1]))

    f.write('0\n')
    for r in conf['bspl_radii']: # Write the bspline radii
        f.write('%.1f ' % (r))
    f.write('\n')
    np.savetxt(f, x_gc[0].reshape((V.shape[1], H.shape[1]))) # Write the perturbation model
    # np.savetxt(f, anis)
    f.close()
