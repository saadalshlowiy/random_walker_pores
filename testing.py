
from math import *
#print(cos(0.3))
import numpy as np
from plotly.offline import init_notebook_mode, iplot
from PIL import Image , ImageSequence
import pandas as pd
import matplotlib.pyplot as plt


#a=np.array([1,2,3,4,5])


#print(f"i want 3 floors , each floor has 4 arrays , each array has 2 elements {np.ones((3 , 4 , 2))}")



def get_array_from_3D_image(path="/home/saad/Desktop/single_pore.tif"):
    img = Image.open(path)
    full_array = np.array([np.array(page) for page in ImageSequence.Iterator(img) ] )
    print(f"3D array has axis {full_array.ndim} with shape {full_array.shape} and has {full_array.size} Voxcels")
    return full_array

help(np.any) 
help(np.all)

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



