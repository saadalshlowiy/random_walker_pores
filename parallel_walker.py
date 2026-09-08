## constants
from math import *
import numpy as np
from PIL import Image , ImageSequence

initial_population_walkers = 10000

current_live_walkers = 10000
# it is much smaller, depending on the resolution of micro-ct
RESOLUTION = 2  # let the lenght of the pixcel be 2 micro-meter
step_distance = RESOLUTION * 0.2
fluid_diffusion_coefficient = 2.5e3 # micro meter^2 / second
surface_relaxivity = 20  # micro-meter/second
p_fraction = []


def temp_sin(delta):
    return sin(delta)
def temp_cos(delta):
    return cos(delta)
def calculate_new_position(current_position ,step_distance,theta, beta) :
    s = step_distance
    (x,y,z) =  current_position
    nx = x + s* temp_sin(beta)* temp_cos(theta)
    ny = y + s* temp_sin(beta)* temp_cos(theta)
    nz = z + s* temp_cos(beta)
    return (nx, ny, nz)

def calculate_increment_time(step_distance ,fluid_diffusion_coefficient ):

    delta_t = step_distance * step_distance / (6 * fluid_diffusion_coefficient)
    return delta_t


def calculate_die_probability(step_distance , fluid_diffusion_coefficient, surface_relaxivity):
    likely = 2 * step_distance * surface_relaxivity / (3 *fluid_diffusion_coefficient )
    return likely



# generatte ramodom threshold , and then compare it to
def is_walker_dead(probability, threshold):
    if probability < threshold :
        return False # he is alive
    return True # he is dead
# SURFACE RELAXATION
def surface_relaxation(p_t):
    # this is mimking T_2 .. starts HIGH , ends LOW
    print("\n\nNOW printing the surface relaxation...\n\n\n")
    for (ratio, _ ) in p_t:
        print(ratio)



# we do not need it , in each pore , we will make one walker in the center ... that is initializing it

def is_initail_walker_in_grain(position):
    # this depends on the geometry of the pore and its boundries
    # for now , just say NO
    return False

    # if we have (3.45 , 2.90 , 12.34) then it is in the voxcel (3,2,12)

#we will store all the position/index of grain voxcels GLOBALLY


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
# from 0 to 2_pi
import math
def random_theta():
    theta = np.random.default_rng().random() * math.pi * 2
    return theta

# from 0 to pie
def random_beta():
    beta = np.random.default_rng().random() * math.pi
    return beta

import time




import threading
def thread_task(boost_number , iterations):
    (p_fraction, coordinates , total_time ) = semi_semi_main(boost_number, iterations)
    final_time = p_fraction[-1][1]
    final_ratio = p_fraction[-1][0]

    surface_relaxation(p_fraction, coordinates)
    os.system(f"cd Figures_Vectorized;mkdir iterations_{iterations}")
    plt.suptitle(f"{iterations} iterations")# edgecolor
    plt.title(f"Final time:{final_time:.5f}sec\nFinal magnetization:{final_ratio:.5f}\nTotal Durations of program:{round(total_time,1)}sec")
    plt.xlabel("Time Duration")
    plt.ylabel("log scale of magnetization")
    plt.tight_layout()
    plt.savefig(f"Figures_Vectorized/iterations_{iterations}/iterations_{iterations}_contraint_{boost_number}_resolution_.2_.png" , dpi = 500 , transparent = False )
    plt.show(block=False)
    plt.pause(6)
    plt.close() # to clear out memory
def master_thread(walker_data,full_array):
    iterations = [76800 , 153600 , 307200 , 614400 , 1228800]
    threads = []
    for i in iterations :
        threads.append(threading.Thread(target = thread_task , args=([0.29] , i  )))
    for thr in threads :
        thr.start()
    for thr in threads :
        thr.join()




def main():
    # initialize walkers then store them
    full_array =  get_array_from_3D_image()
    walker_data = np.ones((initial_population_walkers , 4 )) * ( full_array.shape[0] // 2 )  # so now , they all have the same (x , y , z) which is the center of the image = 121 // 2 = 60
    # now it looks like this ( X Y Z w )    w = 0 , 1
    for walker in walker_data :
        walker[3] = 1


    # now start the simulation

    #(1) first we increment the time
    global current_live_walkers
    global p_fraction

    t = 0

    delta_t = calculate_increment_time(step_distance ,fluid_diffusion_coefficient)

    interations = 20000

    for jj in range(interations):


        # we are iterating through each WALKER
        # every thread has this in their private domain, each process will count how many walkers died !
        number_of_walkers_died_in_one_bash = 0
        for i in  range(len(walker_data)):
            # if this walker DIED , then skip iteration



            if walker_data[i][-1] == 0 :
                continue

            (x , y , z , is_alive) = walker_data[i]
            current_position = (x , y , z)


            theta = random_theta()
            beta = random_beta()
            (nx,ny,nz) = calculate_new_position(current_position ,step_distance,theta, beta)
            if is_collision((nx, ny , nz) , full_array):
            # should be before the two for loops cause it is a constant, but should vary
            # IF IT IN TOUCHING OR BEHOND THE GRAIN, WE CALCULATE TEH LIKELY
                likely = calculate_die_probability(step_distance , fluid_diffusion_coefficient, surface_relaxivity)
                RANDOM_NUMBER = round(np.random.default_rng().random() , 2)  #np.random.random()
                if is_walker_dead(RANDOM_NUMBER , likely):

                    number_of_walkers_died_in_one_bash = number_of_walkers_died_in_one_bash + 1
                    #######        current_live_walkers = current_live_walkers - 1

                    time.sleep(1)
                    print("\n\n!!!!_______WALKER DIED_______!!!!\n\n")

                    is_alive = 0 #    here we mark it DEAD

                    #######        fraction = current_live_walkers / initial_population_walkers # p(t) =  N_1 / N_o
                    #######        p_fraction.append((fraction,t))  # storing( p(t) , t )
                else :
                # we make the walker GO to its previous step !! , if alive
                    (nx,ny,nz) = current_position
            ##
            ## UPDATE THE POSITION OF IT
            walker_data[i] = (nx,ny,nz,is_alive)   # is_alive = 1 :: for being alive      is_alive = 0 for being dead


        # this is done by ONE THREAD ..
        total_number_died= sum(number_of_walkers_died_in_one_bash)  # all THREADS ARE TO BE SUMMED
        current_live_walkers = current_live_walkers - total_number_died
        fraction = current_live_walkers / initial_population_walkers
        p_fraction.append((fraction,t))
        # once we did all walkers in ONE GO , then we update the time
        t = t + delta_t
        print(f"{jj}|T={t}")



    surface_relaxation(p_fraction)

main()
