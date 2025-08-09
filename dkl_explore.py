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
nl_model = "./EC3.0.nlogo"

# Define simulation parameters
# The numbers here MATTER

# Simulation duration
T = 60

# Interval at which to store network configuration
t_rep = 10

# Time decay
tau = 10

# Asintote
asint = 0.5

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

# Propose expected standard deviation decay

# Max outdegree of nodes
friends = 10

def power_law(t, d0: float, dist = float):

    d = (d0 + dist)/np.power(t + 1, (3*friends + 1)/2) + (dist/d0)*np.exp(-t/tau) + asint

    #d = (d0/(np.sqrt(vartrue)))*np.exp(-t/T) + (dist/d0)*np.exp(-t/tau)

    return d

# Divergence function conveniently wrapped

def gauss_DV(x, d: float, nbins: int):

    # Get relative frequencies of the data

    x_freq, llim, wbin, __ = relfreq(x, numbins = nbins)

    # Get gaussian frequencies from the 

    y_freq = []

    for i in range(0, nbins):

        y_freq.append(norm.pdf( (llim + wbin/2 + i*wbin), loc = x.mean(), scale = d))

    y_freq = np.asarray(y_freq)

    dv = kl_div(x_freq, y_freq).sum()

    return dv

# Get agents PDFs and network configuration
def get_internals(N : int, beta : float, dist : float, var_c : float, T = T, t_rep = t_rep):

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
# 

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
    for n in tqdm(range(1,iters+1), desc = 'Running', leave = False):

        netlogo.command('go')
        mus = np.concatenate((mus, values('mu')[np.newaxis, :]), axis = 0)
        sigma2s = np.concatenate((sigma2s, values('var')[np.newaxis, :]), axis = 0)

        if (n%t_rep)==0:
            nets.append(netlogo.report("[list ([label] of end1) ([label] of end2)] of edges").astype(int))

    
    netlogo.kill_workspace()
    return mus, sigma2s, nets



# DEPRECATED use script method pythom -m lightsim

# Create iterable single simulation function for multiprocessing

def rep_sim(N : int, beta : float, dist : float, var_c : float, reps: int, T = T):

    # Set values for simulation

    global_vars['N'] = N
    global_vars['beta'] = beta
    global_vars['var-c'] = var_c
    global_vars['dist'] = dist
    
    global_vars['mutrue'] = mutrue
    global_vars['vartrue'] = vartrue

    netlogo = pynetlogo.NetLogoLink(
          gui = nl_gui,
          netlogo_home = nl_path,
      )

    netlogo.load_model(nl_model)
    

    # Simple function for retrieving nodes variables
    def values(var: str):

        c = netlogo.report(f"map [s -> [{var}] of s] sort nodes")

        return c

    # Results array for multiple simulations
    DV0 = []
    DV = []
 
    # Run multiple simulations
    for r in range(0,reps):

        # Reset NetLogo enviroment
        netlogo.command('clear-all')

        for name in global_vars:

            netlogo.command(f'set {name} {global_vars[name]}')

        netlogo.command('setup')
    
        iters = T

        # Initialize results arrays
        mus = np.empty((0,N))
        sigma2s = mus

        #nets = []

        # Not used for now
        lones = mus
        rewired = mus

        mus = np.concatenate((mus, values('mu0')[np.newaxis, :]), axis = 0)
        #sigma2s = np.concatenate((sigma2s, values('var0')[np.newaxis, :]), axis = 0)

        #nets.append(netlogo.report("[list ([label] of end1) ([label] of end2)] of edges").astype(int))

        # Run the simulation
        for n in range(1,iters+1):

            netlogo.command('go')
            #mus = np.concatenate((mus, values('mu')[np.newaxis, :]), axis = 0)
            #sigma2s = np.concatenate((sigma2s, values('var')[np.newaxis, :]), axis = 0)

            #if (n%t_rep)==0:
                #nets.append(netlogo.report("[list ([label] of end1) ([label] of end2)] of edges").astype(int))

        mus = np.concatenate((mus, values('mu')[np.newaxis, :]), axis = 0)
        #sigma2s = np.concatenate((sigma2s, values('var')[np.newaxis, :]), axis = 0)


        #Get data at the start
        x0 = mus[0,:]

        # Get data at the end
        x = mus[1,:]

        d_0 = power_law(0,d0 = x0.std(), dist = dist)
        d_t = power_law(iters, d0 = x0.std(), dist = dist)

        # Get KL divergence
        dv0 = gauss_DV(x0, d = d_0, nbins = int(N/10))
        dv = gauss_DV(x, d = d_t, nbins =  int(N/10))

        # Append to multiple run results
        DV0.append(dv0)
        DV.append(dv)
    
    netlogo.kill_workspace()

    DV0_var = np.asarray(DV0).var()
    DV_var = np.asarray(DV).var()

    DV0 = np.asarray(DV0).mean()
    DV = np.asarray(DV).mean()
    

    return DV0, DV0_var, DV, DV_var