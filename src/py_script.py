
import os
import datetime
import subprocess
from shutil import copyfile

def network_motifs(file):
  m = motif_structure(str(file), motifsize=3, nrandomizations=0, stoufferIDs=True)
  print_motif_structure(m, outFile='temp.out')
  f = open('temp.out', 'r')
  return(f.read())
  
def modularity_rgraph(file, seed, iter_factor, cooling_factor, randoms):
      modularity_of_randomizations = 0.0
      sd_mod_of_randomizations = 0.0
      
      tmp_dir = '../temp_'+str(datetime.datetime.now()).replace(' ','')+'/'
      os.mkdir(tmp_dir)
      
      copyfile(str(file), tmp_dir+str(file))
      
      file_out = open(tmp_dir+'temp_out.out', 'w')
      try:
          retcode = subprocess.call(['netcarto_cl'], stdout=file_out, stderr=subprocess.STDOUT)
      except OSError:
          print 'Netcarto is not installed or is not available on the PATH in this computer.'
          print 'Cannot obtain the modularity using the rgraph library.'
          print 'Check your rgraph installation.'
          return 0.0, 0
      
      args = ['netcarto_cl', str(file), str(seed), '-1', str(iter_factor), str(cooling_factor), str(randoms)]
      subprocess.call(args, stdout=file_out, stderr=subprocess.STDOUT, cwd=tmp_dir)
      file_out.close()
      
      file_modules = open(tmp_dir+'modules.dat', 'r')
      modules = []
      for line in file_modules:
          a,b,mods = line.partition('---')
          
          if mods != '':
              modules.append(mods.split())
          else:
              start = line.index('=')
              modularity = line[start+1:].strip()
      
      file_modules.close()
                  
      number_of_modules = len(modules)
      modularity = modularity
      
#       #we assign each node to its corresponding module
#       current_module = 1
#       for module in modules:            
#           for n in module:
#               self.node[dict_numbers_nodes[int(n)]]['module'] = current_module
#       
#           current_module += 1
#       
#       #we also find the role of each node and add this information to the network
#       file_roles = open(tmp_dir+'roles.dat', 'r')
#       roles = dict()
#       for line in file_roles:
#           a,b,nodes = line.partition('---')            
#           role = a[0]
#           roles[role] = nodes.split()
#           
#       file_roles.close()
#       
# #       for r in roles.keys():
# #           for n in roles[r]:
# #               self.node[dict_numbers_nodes[int(n)]]['role'] = r
#       
      if randoms > 0:
          file_randoms = open(tmp_dir+'randomized_mod.dat', 'r')
      
          for line in file_randoms:
              if line.startswith('#'):
                  continue
              
              values = line.split()
              modularity_of_randomizations = values[1]
              sd_mod_of_randomizations = values[2]
                  
          file_randoms.close()
      
      return modularity, number_of_modules, modularity_of_randomizations, sd_mod_of_randomizations
