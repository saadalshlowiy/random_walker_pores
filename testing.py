
from math import *
#print(cos(0.3))
import numpy as np
from plotly.offline import init_notebook_mode, iplot
from PIL import Image , ImageSequence
import pandas as pd
import matplotlib.pyplot as plt


#a=np.array([1,2,3,4,5])


#print(f"i want 3 floors , each floor has 4 arrays , each array has 2 elements {np.ones((3 , 4 , 2))}")
def generate_all_in_one_array():



    numberOfWalkers=19


    # (r , c , z , is_alive , nr , nc , nz , is_alive , col , dead)
    r = np.random.default_rng().integers(low=2 , high=70 , size = numberOfWalkers)
    c = np.random.default_rng().integers(low=2 , high=70 , size = numberOfWalkers)
    z = np.random.default_rng().integers(low=2 , high=70 , size = numberOfWalkers)
    is_alive = np.random.default_rng().integers(low=0, high=2 , size = numberOfWalkers)
    nr = np.random.default_rng().integers(low=2 , high=70 , size = numberOfWalkers)
    nc = np.random.default_rng().integers(low=2 , high=70 , size = numberOfWalkers)
    nz = np.random.default_rng().integers(low=2 , high=70 , size = numberOfWalkers)
    is_collid = np.random.default_rng().integers(low=0, high=2 , size = numberOfWalkers)
    is_dead =  np.random.default_rng().integers(low=0, high=2 , size = numberOfWalkers)
    isAlive = is_alive
    all_in_one = np.column_stack((r,c,z,is_alive,nr,nc,nz,isAlive,is_collid,is_dead))


    ## if no collision , updated next step ...

    r[is_collid == 0] = nr[is_collid == 0 ]
    c[is_collid == 0] = nc[is_collid == 0 ]
    z[is_collid == 0] = nz[is_collid == 0 ]



    ##  if collision and is_dead --> is_alive=0
    before_killing = np.sum(is_alive)
    is_alive[(is_dead == 1) & (is_collid == 1)] = 0
    after_killing = np.sum(is_alive)

    number_of_dead_in_this_bach = before_killing - after_killing
    all_in_one = np.column_stack((r,c,z,is_alive,nr,nc,nz,isAlive,is_collid,is_dead))




def get_array_from_3D_image(path="/home/saad/Desktop/single_pore.tif"):
    img = Image.open(path)
    full_array = np.array([np.array(page) for page in ImageSequence.Iterator(img) ] )
    print(f"3D array has axis {full_array.ndim} with shape {full_array.shape} and has {full_array.size} Voxcels")
    return full_array



def convert_position_to_index(position ):
    (x , y , z) = position
    return (int(x) , int(y) , int(z))

def is_index_in_grain(index, full_array ):
    (i , j , k) = index
    value = full_array[i][j][k]

    if value == 0 :
        # for logging purpose
        #print("walker is in grain")
        return True
    return False

def is_collision(position, full_array):
    index = convert_position_to_index(position)
    return is_index_in_grain(index , full_array)


def simulate_3D_array(array_of_image):
    fig, ax = plt.subplots(subplot_kw = {"projection" : "3d"})
    ax.voxels(array_of_image , edgecolor='k')

    plt.show()


def print_2D_array(array):
    for i in array :
        c=""
        for j in i :
            c+=str(j) + " "
        print(c ,"\n")

    return 0


def print_3D_array(array):
    for two_d in array :
        print_2D_array(two_d)
        print("\n\n\n\n")




#from concurrent.future import ProcessPoolExecutor
