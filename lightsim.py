import pynetlogo
import numpy as np
import matplotlib.pyplot as plt
import json

import os,datetime
from tqdm import tqdm
from scipy.special import kl_div
from scipy.stats import relfreq, norm

# Paths

# Netlogo path
nl_path = '/home/lunis/Programs/NetLogo-6.4.0-64'

# Netlogo model
nl_model = "./EC2.9.nlogo"

# Define simulation parameters
# The numbers here MATTER

# Simulation duration
T = 100

# Interval at which to store network configuration
t_rep = 10

# Define simulation parameters
# The numbers here DON'T MATTER
global_vars = {
    
    'N': 1000,
    'beta': 1.,
    'mutrue': 0.,
    'vartrue': 1.,
    'update-type': 2,
    'var-c': 10.,
    'var-d': 2,
    'network-type': "\"scale-free\"",
    'p': 0.01,
    'pref': 1,
    'initial-sampling' : "\"bivariate\"",
    'dist' : 0.

  }

# Could be useful here
mutrue = 0.
vartrue = 1.

# Disable Netlogo gui
nl_gui = False

# Create iterable single simulation function for multiprocessing

def single_sim(N : int, beta : float, dist : float, var_c : float):

    # Set values for simulation

    global_vars['N'] = N
    global_vars['beta'] = beta
    global_vars['var-c'] = var_c
    global_vars['dist'] = dist
    
    global_vars['mutrue'] = mutrue
    global_vars['vartrue'] = vartrue

    # Set path where to save numpy arrays (MAYBE NOT NEEDED)
    outputdir = os.path.join('.','outputs','fisher')
    os.makedirs(outputdir, exist_ok = True)

    netlogo = pynetlogo.NetLogoLink(
          gui = nl_gui,
          netlogo_home = nl_path,
      )

    netlogo.load_model(nl_model)
    

    # Simple function for retrieving nodes variables
    def values(var: str):

        c = netlogo.report(f"map [s -> [{var}] of s] sort nodes")

        return c
    
    

    # Prepare NetLogo enviroment
    netlogo.command('clear-all')

    for name in global_vars:

        netlogo.command(f'set {name} {global_vars[name]}')

    netlogo.command('setup')
 
    iters = T

    # Initialize results arrays
    mus = np.empty((0,N))
    sigma2s = mus

    nets = []

    # Not used for now
    lones = mus
    rewired = mus

    mus = np.concatenate((mus, values('mu0')[np.newaxis, :]), axis = 0)
    sigma2s = np.concatenate((sigma2s, values('var0')[np.newaxis, :]), axis = 0)

    nets.append(netlogo.report("[list ([label] of end1) ([label] of end2)] of edges").astype(int))

    # Run the simulation
    for n in range(1,iters+1):

        netlogo.command('go')
        mus = np.concatenate((mus, values('mu')[np.newaxis, :]), axis = 0)
        sigma2s = np.concatenate((sigma2s, values('var')[np.newaxis, :]), axis = 0)

        if (n%t_rep)==0:
            nets.append(netlogo.report("[list ([label] of end1) ([label] of end2)] of edges").astype(int))


    #Get data at the start
    x0 = mus[0,:]

    # Get data at the end
    x = mus[iters,:]

    #Get KL divergence
    dv0 = gauss_DV(x0, int(N/10))
    dv = gauss_DV(x, int(N/10))

    
    netlogo.kill_workspace()
    return dv0, dv

# Divergence function conveniently wrapped

def gauss_DV(x, nbins: int):

    blind_mean = x.mean()
    blind_dev = vartrue # np.sqrt(x.var()) this method creates problems and changes too much per dataset, results are hard to make sense of

    # Get relative frequencies of the data

    x_freq, llim, wbin, __ = relfreq(x, numbins = nbins)

    # Get gaussian frequencies from the 

    y_freq = np.empty((0))

    for i in range(0, nbins):

        y_freq = np.concatenate((y_freq, norm.pdf( (llim + wbin/2 + i*wbin), loc = blind_mean, scale = blind_dev)[np.newaxis]))

    dv = kl_div(x_freq, y_freq).sum()

    return dv