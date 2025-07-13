import pynetlogo
import numpy as np
import matplotlib.pyplot as plt
import json
import networkx as nx
import os,sys,datetime
from tqdm import tqdm

def main():
  ### SIMULATIONS PARAMETERS
  betas = [0.2,1.,2.]
  varcs = [1.,10.]
  N_iterations = 300
  nl_model = "./EC2.8.nlogo" # modello da caricare in netlogo


  outputdir = os.path.join('.','outputs','update2') # directory dove salvare config e output
  if not os.path.exists(outputdir):
    os.makedirs(outputdir)
  
  nl_gui = False           # false per non far partire la gui di netlogo
  if (sys.platform == 'linux') :
      netlogo = pynetlogo.NetLogoLink(
          gui=nl_gui,
          netlogo_home="/opt/netlogo/",
      )
  else:
      netlogo = pynetlogo.NetLogoLink(
      gui=nl_gui,
  )
  netlogo.load_model(nl_model)

  def values(var: str):
    '''
    Restituisce un array con i valori di una variabile per tutti i nodi.
    '''
    return netlogo.report(f"map [s -> [{var}] of s] sort nodes")
  
  # IMPOSTAZIONE VARIABILI GLOBALI e salvataggio nel file
  global_vars = {
    'N': 1000,
    'beta': 1., # cambierà nei loop
    'mutrue': 0.,
    'vartrue': 1.,
    'update-type': 2,
    'var-c': 10., # cambierà nei loop
    'var-d': 2,
    'network-type': "\"scale-free\"", #random o scale-free
    'p': 0.01, # random network, probabilità di accendere link
    'pref': 1
  }

  with open(os.path.join(outputdir,"global_vars.json"),'w') as f:
    f.write(json.dumps(global_vars,indent=2))
  
  for beta in tqdm(betas,desc='beta'):
    global_vars['beta']=beta
    for varc in tqdm(varcs,desc='varc',leave=False):
      global_vars['var-c']=varc
      # IMPOSTAZIONE VARIABILI GLOBALI E SETUP IN NETLOGO
      netlogo.command('clear-all')
      for name in global_vars:
        netlogo.command(f'set {name} {global_vars[name]}')

      netlogo.command('setup')
      iters = 100
      mus = []
      sigma2s = []
      nets = []
      lones = []
      rewired =[]
      mus.append(values('mu0'))
      sigma2s.append(values('var0'))
      nets.append(netlogo.report("[list ([label] of end1) ([label] of end2)] of edges").astype(int))
      for n in tqdm(range(1,iters+1),desc='iteration',leave=False):
        netlogo.command('go')
        mus.append(values('mu'))
        sigma2s.append(values('var'))
        #lones.append(netlogo.report('lonely'))
        #rewired.append(netlogo.report('rewired'))
        if (n%10)==0:
          nets.append(netlogo.report("[list ([label] of end1) ([label] of end2)] of edges").astype(int))
        #print(f"\r{n}/{iters}",end="",flush=True)
      #print()

      #### Salvataggio
      datestr = datetime.datetime.now().strftime("%m_%d_%H_%M")
      outputsubdir=os.path.join(outputdir,
                          f'beta={global_vars['beta']:.1f}_varc={global_vars['var-c']:.1f}',
                          datestr)

      os.makedirs(outputsubdir)
      np.save(os.path.join(outputsubdir,'mus.npy'),
              mus)
      np.save(os.path.join(outputsubdir,"sigma2s.npy"),
              sigma2s)
      np.save(os.path.join(outputsubdir,"nets.npy"),
              nets)
  netlogo.kill_workspace()
  return
  
if __name__=='__main__':
  main()