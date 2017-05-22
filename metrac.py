#!/usr/bin/env python
"""
Meridional Energy Transport Calculator

Usage:
  metrac.py [options] <INPUT> <OUTPUT>

Options:
  -h --help     Show this screen.
  -c format Variable naming convention. Options 'erai'
"""

# cpT:  [J / kg K] * [K]     = [J / kg]
# Lvq:  [J / kg] * [kg / kg] = [J / kg]
# gz in [m2 / s2] = [ kg m2 / kg s2 ] = [J / kg]

# multiply by v: [J / kg] * [m / s] => [J m / kg s]
# sum over longitudes [J m / kg s] * [ m ] = [J m2 / kg s]

# integrate over pressure: dp: [Pa] = [N m-2] = [kg m2 s-2 m-2] = [kg s-2]
# [J m2 / kg s] * [Pa] = [J m2 / kg s] * [kg / s2] = [J m2 / s3]
# and factor 1/g: [J m2 / s3] * [s2 /m2] = [J / s] = [Wat]

from docopt import docopt
import sys

import netCDF4 as cdf
import numpy as np

import logging
logging.basicConfig(level=logging.DEBUG)

import matplotlib.pyplot as plt

constants = {
  'c_p': 1004.64, # [J / kg K]
  'g':  9.80616,  # [m / s2]
  'L': 2264670.0, # [J / kg]
  'R': 6371000.0  # [m]
}

def calculateDeltas(grid):
  # Calculate dp in [Pa] 
  p = grid['p']
  dp = np.zeros(len(p[:]))
  dp[:] = p[:]

  if p[0] < p[-1]:
    # level 0 is at top of atmosphere
    dp[1:] -= p[:-1]  
  elif p[0] > p[-1]:
    # level 0 is at surface
    dp[:-1] -= p[1:]

  grid['dp'] = dp

  # Calculate dlon
  longitude = grid['longitude']
  latitude = grid['latitude']

  dlon = np.zeros(len(longitude[:]))
  dlon[:-1 ] = longitude[1:] - longitude[:-1]
  dlon[-1] = longitude[0] - longitude[-1]

  # fix periodicity issues
  dlon = np.where(dlon < 0, dlon + 360, dlon)

  # convert to [degrees east] => [m]
  deg_to_km_at_eq = 2.0 * np.pi * constants['R'] / 360.0
  dlon = np.outer(deg_to_km_at_eq * np.cos(latitude * np.pi / 180.0), dlon)

  # fix rounding errors
  dlon = np.where(dlon < 0, 0, dlon)

  grid['dlon'] = dlon

def zonalAndVerticalSum (x, grid, variables):
  """Perform zonal and vertical integral of meridional transport

  Arguments:

  x[level, lat, lon] -- variable to sum
  grid -- Constant variables in dataset: latitude[lat], longitude[lon], etc.
  
  Returns:

  sum[lat] -- integrated variable
  """

  # get variables
  v = variables['v']

  # Sum_level dp[level] * x[level, lat, lon] 
  y = np.tensordot(grid['dp'] / constants['g'], x * v, axes=1)

  # y[lat, lon] * dlon[lat, lon]
  result = np.sum(y * grid['dlon'], axis=1)

  return result

def getConstantVariables (infile, convention):
  """Return a hash containing time, latitude, longitude
  Arguments:
  infile -- netCDF4 Dataset
  convention -- naming convention, ['erai']
  """

  conventionERAI = {
    'time': 'time',           # 
    'p': 'level',             # pressure level [Pa]
    'longitude': 'longitude', # [degrees east]
    'latitude': 'latitude',   # [degrees north]
  }
  factorsERAI = {
    'p': 100.0
  }

  if (convention == 'erai'):
    mapping = conventionERAI
    factors = factorsERAI
  else:
    logging.error('Naming convention not supported')
    sys.exit(1)

  variables = {}
  for key in mapping:
    logging.debug('Mapping ' + mapping[key] + 'to ' + key)
    if (key in factors):
      variables[key] = infile.variables[mapping[key]][...] * factors[key]
    else:
      variables[key] = infile.variables[mapping[key]][...]

  return variables

def massCorrection(grid, variables):
  """subract vertically averaged zonal mean velocity

  Arguments:

  grid -- containing 'dp'
  variables -- containing 'v'
  """

  # zonal mean v[level, lat, lon]
  zonal_mean = np.mean(variables['v'], 2)

  # Vertical (mass weighted) v
  # Sum_level dp[level] * v[level, lat] / Sum_level dp[level]
  correction = np.tensordot(grid['dp'], zonal_mean, axes=1) / np.sum(grid['dp'])

  # v[level, lat, lon] , correction[lat]
  for lat in range(len(grid['latitude'])):
    variables['v'][:,lat,:] -= correction[lat]

def getVariables (infile, convention, time):
  """Return a hash containing t, q, u, v
  Arguments:
  infile -- netCDF4 Dataset
  convention -- naming convention, ['erai']
  time -- time step to retreive
  """

  conventionERAI = {
    't': 't',                 # Temperature [K] 
    # 'u': 'u',                 # eastward wind [m/s]
    'v': 'v',                 # northward wind [m/s]
    'q': 'q',                 # Specific humidity [kg/kg]
    'gz': 'z',                # geopotential height [m2/s2]
  }
  factorsERAI = {
  }

  if (convention == 'erai'):
    mapping = conventionERAI
    factors = factorsERAI
  else:
    logging.error('Naming convention not supported')
    sys.exit(1)

  variables = {}
  for key in mapping:
    logging.debug('Mapping ' + mapping[key] + 'to ' + key)
    if (key in factors):
      variables[key] = infile.variables[mapping[key]][time,...] * factors[key]
    else:
      variables[key] = infile.variables[mapping[key]][time,...]

  return variables

if __name__ == '__main__':
  arguments = docopt(__doc__, version='0.1.1rc')

  # open input file
  try:
    infile = cdf.Dataset(arguments['<INPUT>'], 'r')
  except:
    logging.error('Cannot open file ' + arguments['<INPUT>'])
    sys.exit(1) 

  grid = getConstantVariables(infile, arguments['-c'])
  calculateDeltas(grid)

  transport = np.zeros([len(grid['time']), len(grid['latitude'])])

  # cpT
  for time in range(len(grid['time'])):
    print "full, time: ", time
    variables = getVariables(infile, arguments['-c'], time)
    massCorrection(grid, variables)
    h = constants['c_p'] * variables['t']
    h = np.where(variables['gz'] < 0, 0, h)
    transport[time,:] = zonalAndVerticalSum (h, grid, variables)
  
  plt.plot(grid['latitude'], np.mean(transport, 0) * 1e-15)
  plt.suptitle('v(cpT) dlon dp/g')
  axes = plt.gca()
  axes.set_ylim([-10, 10])
  plt.ylabel('Petawatt')
  plt.savefig('all.png')
  plt.close()

  # L v
  for time in range(len(grid['time'])):
    print "full, time: ", time
    variables = getVariables(infile, arguments['-c'], time)
    massCorrection(grid, variables)
    h = constants['L'] * variables['q']
    h = np.where(variables['gz'] < 0, 0, h)
    transport[time,:] = zonalAndVerticalSum (h, grid, variables)
  
  plt.plot(grid['latitude'], np.mean(transport, 0) * 1e-15)
  plt.suptitle('v(Lq) dlon dp/g')
  axes = plt.gca()
  axes.set_ylim([-10, 10])
  plt.ylabel('Petawatt')
  plt.savefig('all.png')
  plt.close()

  # full
  for time in range(len(grid['time'])):
    print "full, time: ", time
    variables = getVariables(infile, arguments['-c'], time)
    massCorrection(grid, variables)
    h = constants['c_p'] * variables['t'] + constants['L'] * variables['q'] + variables['gz']
    h = np.where(variables['gz'] < 0, 0, h)
    transport[time,:] = zonalAndVerticalSum (h, grid, variables)
  
  plt.plot(grid['latitude'], np.mean(transport, 0) * 1e-15)
  plt.suptitle('v(cpT + Lq + gz) dlon dp/g')
  axes = plt.gca()
  axes.set_ylim([-10, 10])
  plt.ylabel('Petawatt')
  plt.savefig('all.png')
  plt.close()

  # clean up
  infile.close()
