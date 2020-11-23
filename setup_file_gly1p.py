import sys
import os
import numpy as np
import math
from collections import OrderedDict
import time 
import random


filenum = sys.argv[1]


filename= str(filenum) + ".pdb"
a = np.genfromtxt(filename, dtype=str, skip_header=5, skip_footer=2)
atom_gly = 14
atom_wat = 3
num_gly = 1
num_wat = 3550
num = atom_gly * num_gly
numdivatom_gly = int(num/atom_gly)
norm_vec = []
centre = [48/2, 48/2, 48/2]
pos_cO1 = 1
pos_cO2 = 2
pos_cO3 = 3
pos_cH1 = 4
pos_cH2 = 5
pos_cH3 = 6
pos_cen_atom = 1

# current_milli_time = lambda: int(round(time.time() * 1000))
# random.seed(current_milli_time)
# random_3d = random.random(),random.random(),random.random()
# box_c = 5
# random_center = [(centre[0]-box_c)+random_3d[0]*10,(centre[1]-box_c)+random_3d[1]*10,(centre[2]-box_c)+random_3d[2]*10]
# print(random_center)


print(num)

for mol in range(0,num_gly,1):

  vec = [float(a[mol*atom_gly+1,5]),float(a[mol*atom_gly+1,6]),float(a[mol*atom_gly+1,7])]
  norm_vec.append((int(a[mol*atom_gly+1,1]),int(a[mol*atom_gly+1,4]),np.linalg.norm(np.subtract(vec,centre))))
#  print(mol,a[mol,4],vec)

norm_vec = sorted(norm_vec,key=lambda vector: vector[2])

norm_vec = np.array(norm_vec)
print ("central molecule", norm_vec[0,1],"distance",norm_vec[0,2])
vec_ca1 = [float(a[int(norm_vec[0,0])+pos_cO1,5]),float(a[int(norm_vec[0,0])+pos_cO1,6]),float(a[int(norm_vec[0,0])+pos_cO1,7])]
vec_ca2 = [float(a[int(norm_vec[0,0])+pos_cO2,5]),float(a[int(norm_vec[0,0])+pos_cO2,6]),float(a[int(norm_vec[0,0])+pos_cO2,7])]
vec_ca3 = [float(a[int(norm_vec[0,0])+pos_cO3,5]),float(a[int(norm_vec[0,0])+pos_cO3,6]),float(a[int(norm_vec[0,0])+pos_cO3,7])]
vec_ca4 = [float(a[int(norm_vec[0,0])+pos_cH1,5]),float(a[int(norm_vec[0,0])+pos_cH1,6]),float(a[int(norm_vec[0,0])+pos_cH1,7])]
vec_ca5 = [float(a[int(norm_vec[0,0])+pos_cH2,5]),float(a[int(norm_vec[0,0])+pos_cH2,6]),float(a[int(norm_vec[0,0])+pos_cH2,7])]
vec_ca6 = [float(a[int(norm_vec[0,0])+pos_cH3,5]),float(a[int(norm_vec[0,0])+pos_cH3,6]),float(a[int(norm_vec[0,0])+pos_cH3,7])]

#For nearest glycerol
O_dist = {}
O_dist['dist_vec_a1'] = []
O_dist['dist_vec_a2'] = []
O_dist['dist_vec_a3'] = []
O_dist['dist_vec_a4'] = []
O_dist['dist_vec_a5'] = []
O_dist['dist_vec_a6'] = []

O_vec = {}
O_vec['vec_a1'] = []
O_vec['vec_a2'] = []
O_vec['vec_a3'] = []
O_vec['vec_a4'] = []
O_vec['vec_a5'] = []
O_vec['vec_a6'] = []



print(norm_vec.shape[0])
#norm_vec = np.array(norm_vec)
for i in range(1,7):
  for j in range(1,7):

    key_OO = "dist_vec_a" + str(i)
    key_O = "vec_a" + str(j)
    
    for mol in range(0,numdivatom_gly,1):
      
      if int(norm_vec[0,1]) - 1 != mol :
        O_vec[key_O] = [float(a[mol*atom_gly+(2+j),5]),float(a[mol*atom_gly+(2+j),6]),float(a[mol*atom_gly+(2+j),7])]
          
        O_dist[key_OO].append([int(a[mol*atom_gly+(2+j),1]),int(a[mol*atom_gly+(2+j),4]),np.linalg.norm(np.subtract(O_vec[key_O],vec_ca1))])
#          
        O_dist[key_OO].append([int(a[mol*atom_gly+(2+j),1]),int(a[mol*atom_gly+(2+j),4]),np.linalg.norm(np.subtract(O_vec[key_O],vec_ca2))])
#          
        O_dist[key_OO].append([int(a[mol*atom_gly+(2+j),1]),int(a[mol*atom_gly+(2+j),4]),np.linalg.norm(np.subtract(O_vec[key_O],vec_ca3))])
  
        O_dist[key_OO].append([int(a[mol*atom_gly+(2+j),1]),int(a[mol*atom_gly+(2+j),4]),np.linalg.norm(np.subtract(O_vec[key_O],vec_ca4))])
#          
        O_dist[key_OO].append([int(a[mol*atom_gly+(2+j),1]),int(a[mol*atom_gly+(2+j),4]),np.linalg.norm(np.subtract(O_vec[key_O],vec_ca5))])
#          
        O_dist[key_OO].append([int(a[mol*atom_gly+(2+j),1]),int(a[mol*atom_gly+(2+j),4]),np.linalg.norm(np.subtract(O_vec[key_O],vec_ca6))])
 
  O_dist[key_OO] = np.array(O_dist[key_OO])


  O_dist[key_OO] = sorted(O_dist[key_OO],key=lambda vector: vector[2])
  O_dist[key_OO] = np.array(O_dist[key_OO])

Dist_O = O_dist['dist_vec_a1'] + O_dist['dist_vec_a2'] + O_dist['dist_vec_a3'] + O_dist['dist_vec_a4'] + O_dist['dist_vec_a5'] +O_dist['dist_vec_a6']
Dist_O = sorted(O_dist[key_OO],key=lambda vector: vector[2])

print("Dis_O")
#print(Dist_O)
Dist_O = np.array(Dist_O)

mol_id = []
for i in Dist_O:
  mol_id.append(i[1])

#print("mol_id")
mol_id = np.array(mol_id)
#print(mol_id)

uniq_mol_id = list(OrderedDict.fromkeys(mol_id))
uniq_mol_id = np.array(uniq_mol_id)
#print(uniq_mol_id)

#nn_O1 = int(uniq_mol_id[0])
#nn_O2 = int(uniq_mol_id[1])

#print(nn_O1,nn_O2)
#print(Dist_O[0], Dist_O[1])



f = open( str(filenum) +  "_1gly-3wat_1p_C-edge_pople3z.inp","a")
print(f)

f.write(F'$molecule\n')
f.write(F'0\t1\n')

for mol in range(0,num_gly):
#for mol in range(0,3):

  for glycerol in range(0,atom_gly):
    
    if mol == 0:
        name = a[(int(norm_vec[0,1])-1)*atom_gly+glycerol,2][0]
        f.write(F'{name[0]}\t {a[(int(norm_vec[0,1])-1)*atom_gly+glycerol,5]}\t{a[(int(norm_vec[0,1])-1)*atom_gly+glycerol,6]}\t{a[(int(norm_vec[0,1])-1)*atom_gly+glycerol,7]}\t')
    else:
        name = a[(int(uniq_mol_id[mol-1])-1)*14+glycerol,2][0]
        f.write(F'{name[0]}\t {a[(int(uniq_mol_id[mol-1])-1)*14+glycerol,5]}\t{a[(int(uniq_mol_id[mol-1])-1)*14+glycerol,6]}\t{a[(int(uniq_mol_id[mol-1])-1)*14+glycerol,7]}\t')
   
    if glycerol == 0:
      f.write(F'-1\t{(mol*14)+2}\t{(mol*14)+10}\t{(mol*14)+11}\t{(mol*14)+4}\n')
    elif glycerol == 1:
      f.write(F'-2\t{(mol*14)+1}\t{(mol*14)+3}\t{(mol*14)+6}\t{(mol*14)+12}\n')
    elif glycerol == 2:
      f.write(F'-1\t{(mol*14)+2}\t{(mol*14)+13}\t{(mol*14)+14}\t{(mol*14)+5}\n')
    elif glycerol == 3:
      f.write(F'-6\t{(mol*14)+1}\t{(mol*14)+7}\t{0}\t{0}\n')
    elif glycerol == 4:
      f.write(F'-6\t{(mol*14)+3}\t{(mol*14)+9}\t{0}\t{0}\n')
    elif glycerol == 5:
      f.write(F'-7\t{(mol*14)+2}\t{(mol*14)+8}\t{0}\t{0}\n')
    elif glycerol == 6:
      f.write(F'-5\t{(mol*14)+4}\t{0}\t{0}\t{0}\n')
    elif glycerol == 7:
      f.write(F'-5\t{(mol*14)+6}\t{0}\t{0}\t{0}\n')
    elif glycerol == 8:
      f.write(F'-5\t{(mol*14)+5}\t{0}\t{0}\t{0}\n')
    elif glycerol == 9:
      f.write(F'-3\t{(mol*14)+1}\t{0}\t{0}\t{0}\n')
    elif glycerol == 10:
      f.write(F'-3\t{(mol*14)+1}\t{0}\t{0}\t{0}\n') 
    elif glycerol == 11:
      f.write(F'-4\t{(mol*14)+2}\t{0}\t{0}\t{0}\n')
    elif glycerol == 12:
      f.write(F'-3\t{(mol*14)+3}\t{0}\t{0}\t{0}\n')
    elif glycerol == 13:
      f.write(F'-3\t{(mol*14)+3}\t{0}\t{0}\t{0}\n')


#For nearest water     

wat_dist = {}
wat_dist['dist_vec_a1'] = []
wat_dist['dist_vec_a2'] = []
wat_dist['dist_vec_a3'] = []
wat_dist['dist_vec_a4'] = []
wat_dist['dist_vec_a5'] = []
wat_dist['dist_vec_a6'] = []

wat_vec = {}
wat_vec['vec_a1'] = []
wat_vec['vec_a2'] = []
wat_vec['vec_a3'] = []

print(norm_vec.shape[0])
#norm_vec = np.array(norm_vec)
for i in range(1,7):
  for j in range(1,4):

    key_OO = "dist_vec_a" + str(i)
    key_O = "vec_a" + str(j)
    
    for mol in range(0,num_wat,1):

        index = num + mol * atom_wat 
      
        wat_vec[key_O] = [float(a[index+(j-1),5]),float(a[index+(j-1),6]),float(a[index+(j-1),7])]
          
        wat_dist[key_OO].append([int(a[index+(j-1),1]),int(a[index+(j-1),4]),np.linalg.norm(np.subtract(wat_vec[key_O],vec_ca1))])
#          
        #wat_dist[key_OO].append([int(a[index+(j-1),1]),int(a[index+(j-1),4]),np.linalg.norm(np.subtract(wat_vec[key_O],vec_ca2))])
#          
        #wat_dist[key_OO].append([int(a[index+(j-1),1]),int(a[index+(j-1),4]),np.linalg.norm(np.subtract(wat_vec[key_O],vec_ca3))])
  
        wat_dist[key_OO].append([int(a[index+(j-1),1]),int(a[index+(j-1),4]),np.linalg.norm(np.subtract(wat_vec[key_O],vec_ca4))])
#          
        #wat_dist[key_OO].append([int(a[index+(j-1),1]),int(a[index+(j-1),4]),np.linalg.norm(np.subtract(wat_vec[key_O],vec_ca5))])
#          
        #wat_dist[key_OO].append([int(a[index+(j-1),1]),int(a[index+(j-1),4]),np.linalg.norm(np.subtract(wat_vec[key_O],vec_ca6))])
 
  wat_dist[key_OO] = np.array(wat_dist[key_OO])

  wat_dist[key_OO] = sorted(wat_dist[key_OO],key=lambda vector: vector[2])
  wat_dist[key_OO] = np.array(wat_dist[key_OO])

#Dist_wat = wat_dist['dist_vec_a1'] + wat_dist['dist_vec_a2'] + wat_dist['dist_vec_a3'] + wat_dist['dist_vec_a4'] + wat_dist['dist_vec_a5'] + wat_dist['dist_vec_a6']
Dist_wat = wat_dist['dist_vec_a1'] + wat_dist['dist_vec_a4']
Dist_wat = sorted(wat_dist[key_OO],key=lambda vector: vector[2])

#print("Dis_wat")
#print(Dist_wat)
Dist_wat = np.array(Dist_wat)

mol_id_wat = []
for i in Dist_wat:
  mol_id_wat.append(i[1])

print("mol_id")
mol_id_wat = np.array(mol_id_wat)
print(mol_id_wat)

uniq_mol_id_wat = list(OrderedDict.fromkeys(mol_id_wat))
uniq_mol_id_wat = np.array(uniq_mol_id_wat)
print("uniq_mol_id_wat")
print(uniq_mol_id_wat)

#nn_O1_wat = int(uniq_mol_id_wat[0])
#nn_O2_wat = int(uniq_mol_id_wat[1])

#print("uniq_mol_id_wat length", uniq_mol_id_wat.shape)

for mol in range(0,num_wat):
#for mol in range(0,5):
  
  for water in range(0,atom_wat):

    index = num + mol * atom_wat + water
    #index_wat = num + (int(uniq_mol_id_wat[mol])-1)*atom_wat+water
    index_wat = num + (int(uniq_mol_id_wat[mol])-(num_gly+1))*atom_wat + water
    name = a[index_wat,2][0]
    #print("uniq_mol_id_wat", uniq_mol_id_wat[mol], "index_wat", mol, index_wat)

    f.write(F'{name[0]}\t {a[index_wat,5]}\t{a[index_wat,6]}\t{a[index_wat,7]}\t')


    if water == 0:
      f.write(F'2001\t{index+2}\t{index+3}\t0\t0\n')
    elif water == 1:
      f.write(F'2002\t{index}\t0\t0\t0\n')
    elif water == 2:
      f.write(F'2002\t{index-1}\t0\t0\t0\n')
  #print(mol)
f.write(F'$end')

f.close    
