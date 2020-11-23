import sys
import os
import numpy as np
import math
from collections import OrderedDict
import time 
import random
import glob



def distance(C1,C2):

  q = np.subtract(C2,C1)  #C2-C1
  dist = np.linalg.norm(q)

  return dist

def unit(C1,C2):

  q = np.subtract(C2,C1)  #C2-C1
  u = q/np.sqrt(np.dot(q,q))

  return u


def angle(C1,C2,C3):

  theta = math.acos(np.dot(unit(C2,C1),unit(C2,C3)))

  theta_deg = np.degrees(theta)

  return theta_deg


def nearest_wH(a,num,mol,atom_wat,num_wat,vec_ca):

    wat_dist = {}
    wat_dist['dist_vec_a1'] = []
    

    wat_vec = {}
    wat_vec['vec_a2'] = []
    wat_vec['vec_a3'] = []

    for i in range(1,2):
        for j in range(2,4):

            key_OO = "dist_vec_a" + str(i)
            key_O = "vec_a" + str(j)

            for mol in range(0,num_wat,1):

                index = num + mol * atom_wat
                
                wat_vec[key_O] = [float(a[index+(j-1),5]),float(a[index+(j-1),6]),float(a[index+(j-1),7])]

                wat_dist[key_OO].append([int(a[index+(j-1),1]),int(a[index+(j-1),4]),np.linalg.norm(np.subtract(wat_vec[key_O],vec_ca))])

                

        wat_dist[key_OO] = np.array(wat_dist[key_OO])

        wat_dist[key_OO] = sorted(wat_dist[key_OO],key=lambda vector: vector[2])
        wat_dist[key_OO] = np.array(wat_dist[key_OO])

    Dist_wat = wat_dist["dist_vec_a1"]
    Dist_wat = np.array(Dist_wat)

    return Dist_wat


def nearest_wO(a,num,mol,atom_wat,num_wat,vec_ca):

    wat_dist = {}
    wat_dist['dist_vec_a1'] = []
    

    wat_vec = {}
    wat_vec['vec_a1'] = []

    for i in range(1,2):
        for j in range(1,2):

            key_OO = "dist_vec_a" + str(i)
            key_O = "vec_a" + str(j)

            for mol in range(0,num_wat,1):

                index = num + mol * atom_wat
                
                wat_vec[key_O] = [float(a[index+(j-1),5]),float(a[index+(j-1),6]),float(a[index+(j-1),7])]

                wat_dist[key_OO].append([int(a[index+(j-1),1]),int(a[index+(j-1),4]),np.linalg.norm(np.subtract(wat_vec[key_O],vec_ca))])

                

        wat_dist[key_OO] = np.array(wat_dist[key_OO])

        wat_dist[key_OO] = sorted(wat_dist[key_OO],key=lambda vector: vector[2])
        wat_dist[key_OO] = np.array(wat_dist[key_OO])

    Dist_wat = wat_dist["dist_vec_a1"]
    Dist_wat = np.array(Dist_wat)

    return Dist_wat

def Og_Hw(a,num,mol,atom_wat,num_wat,num_gly,vec_ca,dist_cutoff,angle_cutoff):

    nearest_watH = nearest_wH(a,num,mol,atom_wat,num_wat,vec_ca)

    num_Og_Hw = 0
    for H in range(0,10):

        index_O = num + (int(nearest_watH[H,1])-(num_gly+1))*atom_wat + 0
        #print("nearest_watH[H,0]",nearest_watH[H,0])
        #print("index_O",index_O)
        index_H = nearest_watH[H,0]-1
        coord_O = [float(a[int(index_O),5]),float(a[int(index_O),6]),float(a[int(index_O),7])]
        coord_H = [float(a[int(index_H),5]),float(a[int(index_H),6]),float(a[int(index_H),7])]
        
        #print("coord_H",coord_H)

        dist_Og_Hw = distance(vec_ca,coord_O)
        angle_Og_Hw_Ow = angle(vec_ca,coord_O,coord_H)

        # print('H',H)
        # print('dist_Og_Hw',dist_Og_Hw)
        # print('angle_Og_Hw_Ow',angle_Og_Hw_Ow)
        # print('num_Og_Hw',num_Og_Hw)

        if dist_Og_Hw <= dist_cutoff and angle_Og_Hw_Ow <= angle_cutoff:
            num_Og_Hw += 1

                
    return num_Og_Hw


def Hg_Ow(a,num,mol,atom_wat,num_wat,num_gly,vec_ca,vec_cao,dist_cutoff,angle_cutoff):

    nearest_watO = nearest_wO(a,num,mol,atom_wat,num_wat,vec_ca)

    num_Hg_Ow = 0
    for O in range(0,10):

        #index_Og = norm_vec[0,0]+pos_cO1  #not needed here
       #print("nearest_watH[H,0]",nearest_watH[H,0])
        #print("index_O",index_O)
        index_O = nearest_watO[O,0]-1
        coord_Og = vec_cao
        coord_O = [float(a[int(index_O),5]),float(a[int(index_O),6]),float(a[int(index_O),7])]
        
        #print("coord_H",coord_H)

        dist_Hg_Ow = distance(coord_Og,coord_O)
        angle_Og_Hg_Ow = angle(vec_ca,coord_Og,coord_O)

        # print('O',O)
        # print('dist_Og_Hw',dist_Hg_Ow)
        # print('angle_Og_Hw_Ow',angle_Og_Hg_Ow)
        # print('num_Og_Hw',num_Hg_Ow)

        if dist_Hg_Ow <= dist_cutoff and angle_Og_Hg_Ow <=  angle_cutoff:
            num_Hg_Ow += 1

    return num_Hg_Ow

def Og_Hg(a,vec_ca1,vec_ca2,vec_ca3,vec_ca4,vec_ca5,vec_ca6,dist_cutoff,angle_cutoff,types):

    
    num_Og_Hg = 0

    if types == 1: 
        dist_Og_Hg = distance(vec_ca1,vec_ca3)
        angle_Og_Hg_Og = angle(vec_ca1,vec_ca3,vec_ca5)
    if types == 2:
        dist_Og_Hg = distance(vec_ca3,vec_ca2)
        angle_Og_Hg_Og = angle(vec_ca3,vec_ca2,vec_ca6)
    if types == 3:
        dist_Og_Hg = distance(vec_ca2,vec_ca1)
        angle_Og_Hg_Og = angle(vec_ca2,vec_ca1,vec_ca4)        
    print('dist_Og_Hg',dist_Og_Hg, 'angle_Og_Hg_Og',angle_Og_Hg_Og)

    if dist_Og_Hg <= dist_cutoff and  angle_Og_Hg_Og <=  angle_cutoff:
        num_Og_Hg += 1
       

    if types == 1:
        dist_Og_Hg = distance(vec_ca1,vec_ca2)
        angle_Og_Hg_Og = angle(vec_ca1,vec_ca2,vec_ca6)
    if types == 2:
        dist_Og_Hg = distance(vec_ca3,vec_ca1)
        angle_Og_Hg_Og = angle(vec_ca3,vec_ca1,vec_ca4)
    if types == 3:
        dist_Og_Hg = distance(vec_ca2,vec_ca3)
        angle_Og_Hg_Og = angle(vec_ca2,vec_ca3,vec_ca5)       
    print('dist_Og_Hg2',dist_Og_Hg, 'angle_Og_Hg_Og2',angle_Og_Hg_Og)
        

    if dist_Og_Hg <= dist_cutoff and  angle_Og_Hg_Og <=  angle_cutoff:
        num_Og_Hg += 1
                
    return num_Og_Hg

def Hg_Og(a,vec_ca1,vec_ca2,vec_ca3,vec_ca4,vec_ca5,vec_ca6,dist_cutoff,angle_cutoff,types):

    
    num_Hg_Og = 0
    
    if types == 1: 
        dist_Hg_Og = distance(vec_ca1,vec_ca3)
        angle_Og_Hg_Og = angle(vec_ca4,vec_ca1,vec_ca3)
    if types == 2:
        dist_Hg_Og = distance(vec_ca3,vec_ca2)
        angle_Og_Hg_Og = angle(vec_ca5,vec_ca3,vec_ca2)
    if types == 3:
        dist_Hg_Og = distance(vec_ca2,vec_ca1)
        angle_Og_Hg_Og = angle(vec_ca6,vec_ca2,vec_ca1)        
    print('dist_Hg_Og',dist_Hg_Og, 'angle_Og_Hg_Og',angle_Og_Hg_Og)
    

    if dist_Hg_Og <= dist_cutoff and  angle_Og_Hg_Og <=  angle_cutoff:
        num_Hg_Og += 1

    if types == 1: 
        dist_Hg_Og = distance(vec_ca1,vec_ca2)
        angle_Og_Hg_Og = angle(vec_ca4,vec_ca1,vec_ca2)
    if types == 2:
        dist_Hg_Og = distance(vec_ca3,vec_ca1)
        angle_Og_Hg_Og = angle(vec_ca5,vec_ca3,vec_ca1)
    if types == 3:
        dist_Hg_Og = distance(vec_ca2,vec_ca3)
        angle_Og_Hg_Og = angle(vec_ca6,vec_ca2,vec_ca3)        
    print('dist_Hg_Og',dist_Hg_Og, 'angle_Og_Hg_Og',angle_Og_Hg_Og)
    


    if dist_Hg_Og <= dist_cutoff and  angle_Og_Hg_Og <=  angle_cutoff:
        num_Hg_Og += 1
                
    return num_Hg_Og

def h_bond_calc(filename):

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

    dist_cutoff = 3.7
    angle_cutoff = 50

    types = 1

    num_Og_Hw = 0
    num_Hg_Ow = 0
    num_Og_Hg = 0
    num_Hg_Og = 0


    #print(num)

    for mol in range(0,num_gly,1):

        vec = [float(a[mol*atom_gly+1,5]),float(a[mol*atom_gly+1,6]),float(a[mol*atom_gly+1,7])]
        norm_vec.append((int(a[mol*atom_gly+1,1]),int(a[mol*atom_gly+1,4]),np.linalg.norm(np.subtract(vec,centre))))

    norm_vec = sorted(norm_vec,key=lambda vector: vector[2])

    norm_vec = np.array(norm_vec)
    #print ("central molecule", norm_vec[0,1],"distance",norm_vec[0,2])
    vec_ca1 = [float(a[int(norm_vec[0,0])+pos_cO1,5]),float(a[int(norm_vec[0,0])+pos_cO1,6]),float(a[int(norm_vec[0,0])+pos_cO1,7])]
    vec_ca2 = [float(a[int(norm_vec[0,0])+pos_cO2,5]),float(a[int(norm_vec[0,0])+pos_cO2,6]),float(a[int(norm_vec[0,0])+pos_cO2,7])]
    vec_ca3 = [float(a[int(norm_vec[0,0])+pos_cO3,5]),float(a[int(norm_vec[0,0])+pos_cO3,6]),float(a[int(norm_vec[0,0])+pos_cO3,7])]
    vec_ca4 = [float(a[int(norm_vec[0,0])+pos_cH1,5]),float(a[int(norm_vec[0,0])+pos_cH1,6]),float(a[int(norm_vec[0,0])+pos_cH1,7])]
    vec_ca5 = [float(a[int(norm_vec[0,0])+pos_cH2,5]),float(a[int(norm_vec[0,0])+pos_cH2,6]),float(a[int(norm_vec[0,0])+pos_cH2,7])]
    vec_ca6 = [float(a[int(norm_vec[0,0])+pos_cH3,5]),float(a[int(norm_vec[0,0])+pos_cH3,6]),float(a[int(norm_vec[0,0])+pos_cH3,7])]    

    
    #num_Og_Hw = Og_Hw(a,num,mol,atom_wat,num_wat,num_gly,vec_ca1,dist_cutoff,angle_cutoff)

    num_Hg_Ow = Hg_Ow(a,num,mol,atom_wat,num_wat,num_gly,vec_ca4,vec_ca1,dist_cutoff,angle_cutoff)
   
    #num_Og_Hg = Og_Hg(a,vec_ca1,vec_ca2,vec_ca3,vec_ca4,vec_ca5,vec_ca6,dist_cutoff,angle_cutoff,types)

    #num_Hg_Og = Hg_Og(a,vec_ca1,vec_ca2,vec_ca3,vec_ca4,vec_ca5,vec_ca6,dist_cutoff,angle_cutoff,types)


    total_hb = num_Og_Hw + num_Og_Hg + num_Hg_Ow + num_Hg_Og

    return total_hb


def main():


    f = open( "hbond_Hg1-O.dat","a")

    for name in glob.glob('snapshots*.pdb'):

        print(name)

        hbonds = h_bond_calc(name)
        
        f.write(F'{hbonds}\n')



    f.close




if __name__ == "__main__":
    main()