## constants
from math import *
import numpy as np
from PIL import Image , ImageSequence
import matplotlib.pyplot as plt


initial_population_walkers = 10000
current_live_walkers = initial_population_walkers
placed_walkers = 0
# it is much smaller, depending on the resolution of micro-ct
RESOLUTION = 0.2  # let the lenght of the pixcel be 2 micro-meter

step_distance = RESOLUTION * 0.2  # s = 0.2 X L .... L is lenght of one pixcel in micro meteres
fluid_diffusion_coefficient = 2.5e3 # micro meter^2 / second
surface_relaxivity = 20  # micro-meter/second



RADIUS = 50 # pixcels given from doctor 

RADIUS_IN_MICRO_METERS = RADIUS * RESOLUTION


def analytical_approach(surface_relaxivity , final_time , radius_in_micro_meter, iterations ):
    # we are going to model the equation  M/M = exp(-3pt/r) , were t=0 --> t=final_t || final_t is aquired AFTER the simulation (when did we stop)
    # this method is used as a benshmark towards the simulation ...
    t = np.linspace(0 , final_time ,iterations )

    y =  ((-3 * t * surface_relaxivity) / radius_in_micro_meter )

    #y = [math.exp(e) for e in y]
    # SINCE WE ARE TAKING THE LOG() OF IT , THEN NOO NEED TO TAKE EXP()
    return (t , y )



def temp_sin(delta):
    return sin(delta)
def temp_cos(delta):
    return cos(delta)
def calculate_new_position(current_position ,step_distance,theta, beta) :
    s = step_distance
    (r,c,z) =  current_position
    nr = r + s* temp_sin(beta)* temp_cos(theta)
    nc = c + s* temp_sin(beta)* temp_cos(theta)
    nz = z + s* temp_cos(beta)
    return (nr, nc, nz)





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
def surface_relaxation_for_simulation(p_t):
    # this is mimking T_2 .. starts HIGH , ends LOW
    print("\n\nNOW printing the surface relaxation...\n\n\n")

    number_of_elements = len(p_t)
    # pre-declaration

    x =  np.empty(number_of_elements)
    y =  np.empty(number_of_elements)

    for i in range(number_of_elements):
        (ratio , time) = p_t[i]
        x[i] = time
        y[i] =  math.log(ratio)
        #y[i] =  ratio
        #print(f"x:{time} || y:{ratio}")





    plt.plot(x,y ,'purple' )


def surface_relaxation_for_analytical(coordinates):
    (x , y) = coordinates
    plt.plot(x,y , 'o')


def surface_relaxation(p_t , coordinates):
    plt.figure()
    surface_relaxation_for_simulation(p_t)
    surface_relaxation_for_analytical(coordinates)
    plt.show()





    # if we have (3.45 , 2.90 , 12.34) then it is in the voxcel (3,2,12)

#we will store all the position/index of grain voxcels GLOBALLY


def get_array_from_3D_image(path="/home/saad/Desktop/single_pore.tif"):
    img = Image.open(path)
    full_array = np.array([np.array(page) for page in ImageSequence.Iterator(img) ] )
    print(f"\n3D array has axis {full_array.ndim} with shape {full_array.shape} and has {full_array.size} Voxcels\n")
    return full_array

def get_indexes_of_pore_from_3D_array(array):
    indexes = [] # unknown number   THUS THE INDEX ARE IN FORMAT --> [ (r , c , z) ... ]
    for z in range(len(array)):
        for r in range(len(array[z])):
            for c in range(len(array[z][r])):
                if array[z][r][c] == 1 :    # pores are 1     grains are 0
                    indexes.append((r , c , z))

    return indexes


def convert_position_to_index(position ):
    (r , c , z) = position
    return (int(r) , int(c) , int(z))



        ## walker_data      full_array



def is_index_in_grain(index, full_array ):
    (r , c , z) = index
    if r > 120 or c > 120 or z > 120 :
        return True
    #value = full_array[i][j][k] # (i=z j=x k=y)
    value = full_array[z][r][c]  # cause from walker data (x y z) --> (i j k) ..and array[z][x][y]
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

def convert_index_to_mid_point_position(index):
    (r , c , z) = index
    r+=0.5
    c+=0.5
    z+=0.5
    return (r , c , z)


def get_initiated_walkers(full_array , initial_number_of_walkers):
    pores_indexes = get_indexes_of_pore_from_3D_array(full_array)
    number_of_pores = len(pores_indexes)
    walker_data = []

    for i in range(initial_number_of_walkers):
        # first , for this wakler , we need to find a random PORE
        #     1      we randomly PICK a pore
        i = int(np.random.default_rng().random() * number_of_pores)
        index = pores_indexes[i]

        # 2 we calculate the MIDDLE POSITION OF THIS PORE (4 , 5 , 8) --> (4.5 , 5.5 , 8.5)
        position_for_uninitialized_walker = convert_index_to_mid_point_position(index) + (1,)      # so now it is like this (x , y , z , is_alive)

        walker_data.append(position_for_uninitialized_walker)


    return walker_data


def return_collided_walkers(converted_to_index_new_walker_data , full_array) :
            is_collision = np.zeros(len(converted_to_index_new_walker_data))

            for i in range(len(converted_to_index_new_walker_data)) :

            ## this is a BUG , WE ARE GETTING NEGITIVE VALUES AND WE DO NOT KNOW WHY
                (r , c , z , a) = np.abs(converted_to_index_new_walker_data[i])
                #print("index of pixcels is :" , r , c , z , a)
                if r > 120 or c > 120 or z > 120 :
                    is_collision[i] = 1
                    continue

                if full_array[int(z)][int(r)][int(c)] == 0 :
                    is_collision[i] = 1

            return is_collision

def update_walker_state(walker_data , new_walker_data  , is_dead , collided ) :


        number_of_dead_in_this_bach=0
        for i in range(len(is_dead)):

            # if there is no collision
            if collided[i] == 0 :
                walker_data[i] = new_walker_data[i]
            else :
                #if collided and dead
                #print(f"now it might be dead or prev dead , we will check--> {is_dead[i]} or {walker_data[i][3]}")
                #time.sleep(1.75)
                prev_state = walker_data[i][3] # if 1 , then it was alive .... if zero , it was dead long time ago
                if is_dead[i] == True  :
                    walker_data[i][3] = 0 # mark as dead
                    if prev_state == 1 :
                        number_of_dead_in_this_bach+=1      # if he WAS ALIVE , AND NOW WE MARK HIM DEAD , then we increemnt
                else:
                    #if collided and did not die , then do not do anything about the state
                    pass

        return number_of_dead_in_this_bach

def semi_main():
    # initialize walkers then store them
    full_array =  get_array_from_3D_image()


     # now it looks like this ( X Y Z w )    w = 0 , 1
    walker_data = np.array(get_initiated_walkers(full_array , initial_population_walkers ))




    # now start the simulation

    #(1) first we increment the time
    global current_live_walkers

    global p_fraction

    t = 0

    delta_t = calculate_increment_time(step_distance ,fluid_diffusion_coefficient)

    interations = 1200
    p_fraction = []

    start = time.time()
    for jj in range(interations):


        # we are iterating through each WALKER
        # every thread has this in their private domain, each process will count how many walkers died !


            # if this walker DIED , then skip iteration




            aaa = time.time()

            r = walker_data[: , 0]
            c = walker_data[: , 1]
            z = walker_data[: , 2]
            is_alive =  walker_data[: , 3]


           # theta  = np.ones(len(walker_data))  # [ th , th , th ... th ]
            #beta   = np.ones(len(walker_data))  # [ be , be , be ... be
            theta =  np.random.default_rng().random( len(walker_data) ) * math.pi * 2
            beta  =  np.random.default_rng().random( len(walker_data) ) * math.pi
        # WE EXTRACT X Y Z SEPRATELY

            nr =  r + (step_distance * np.sin(beta) * np.cos(theta))
            nc =  c + (step_distance * np.sin(beta) * np.cos(theta))
            nz =  z + ( step_distance * np.cos(beta) )
            new_walker_data = np.column_stack((nr , nc , nz , is_alive))
            #if nr < 0 or nc < 0 or nz < 0 :
            #    print(f"new position is negitive")
            #r=round(np.random.default_rng().random() , 2)
            #if r  > 0.89 :

             #   print(f"{jj}|~{r}~ new position calculated ({x} {y} {z}) --> ({nx} {ny} {nz})...")
              #  time.sleep(0.25)

            converted_to_index_new_walker_data = new_walker_data // 1       # from POSITION --> Index

            collided = return_collided_walkers(converted_to_index_new_walker_data , full_array)

            # should be before the two for loops cause it is a constant, but should vary
            # IF IT IN TOUCHING OR BEHOND THE GRAIN, WE CALCULATE TEH LIKELY
               # print(f"walker has HIT a grain , now calculating if it is dead.......?")
                #time.sleep(0.5)
            likely = np.ones(len(collided)) * ( 2 * step_distance * surface_relaxivity / (3 *fluid_diffusion_coefficient ) )
            random_number = np.random.default_rng().random(size = len(collided))
            is_dead = random_number > likely ## if random number > likely    then DEAD      # output is [ False , True , False ..... etc ]

                    #print("walker dead!!!!")


            ## O(n) ... we need to loop Through the entire set , to update ..
            # i do not know if we can do it , with just a if else statement on np.array
            number_of_dead_in_this_bach = update_walker_state(walker_data , new_walker_data  , is_dead , collided )

            current_live_walkers = current_live_walkers - number_of_dead_in_this_bach




            fraction = current_live_walkers / initial_population_walkers # p(t) =  N_1 / N_o

            t = t + delta_t



            if current_live_walkers == 0 :
                break


            p_fraction.append((fraction,t))


            bbb = time.time()
            print(f"{jj}|| T:{t} || F:{fraction} || {bbb - aaa}")

            # is_alive = 1 :: for being alive      is_alive = 0 for being dead




    end = time.time()

    print(f"\n\n\nOVERALL TIME {end - start} ")
    print(f"final time that is going to anaylitical approach is {t}")
    coordinates = analytical_approach(surface_relaxivity , t , RADIUS_IN_MICRO_METERS, interations )

    print(f"\n\n\nSTEP DISTANCE USED IS {step_distance/RESOLUTION} X L")
    surface_relaxation(p_fraction, coordinates)




def main():
    # initialize walkers then store them
    full_array =  get_array_from_3D_image()

     # now it looks like this ( X Y Z w )    w = 0 , 1
    walker_data = np.array(get_initiated_walkers(full_array , initial_population_walkers ))




    # now start the simulation

    #(1) first we increment the time
    global current_live_walkers

    global p_fraction

    t = 0

    delta_t = calculate_increment_time(step_distance ,fluid_diffusion_coefficient)

    interations = 1200
    p_fraction = []

    start = time.time()
    for jj in range(interations):


        # we are iterating through each WALKER
        # every thread has this in their private domain, each process will count how many walkers died !
        number_of_walkers_died_in_one_bash = 0

        aaa = time.time()
        for i in  range(len(walker_data)):
            # if this walker DIED , then skip iteration



            if walker_data[i][-1] == 0 :
                continue

            (r , c , z , is_alive) = walker_data[i]
            current_position = (r , c, z)


            theta = random_theta()
            beta = random_beta()
            (nr,nc,nz) = calculate_new_position( current_position , step_distance , theta , beta )
            #r=round(np.random.default_rng().random() , 2)
            #if r  > 0.89 :

             #   print(f"{jj}|~{r}~ new position calculated ({x} {y} {z}) --> ({nx} {ny} {nz})...")
              #  time.sleep(0.25)
            if is_collision((nr, nc , nz) , full_array):
            # should be before the two for loops cause it is a constant, but should vary
            # IF IT IN TOUCHING OR BEHOND THE GRAIN, WE CALCULATE TEH LIKELY
               # print(f"walker has HIT a grain , now calculating if it is dead.......?")
                #time.sleep(0.5)
                likely = calculate_die_probability(step_distance , fluid_diffusion_coefficient, surface_relaxivity)
                RANDOM_NUMBER = round(np.random.default_rng().random() , 2)  #np.random.random()
                if is_walker_dead(RANDOM_NUMBER , likely):
                    #print("walker dead!!!!")
                    number_of_walkers_died_in_one_bash = number_of_walkers_died_in_one_bash + 1
                    current_live_walkers = current_live_walkers - 1


                    #print("\n_______WALKER DIED_______\n")
                    #print(f"{jj}| new position calculated ({r} {c} {z}) --> ({nr} {nc} {nz})...")
                    #time.sleep(0.005)
                    is_alive = 0 #    here we mark it DEAD

            #        fraction = current_live_walkers / initial_population_walkers # p(t) =  N_1 / N_o
             #       p_fraction.append((fraction,t))  # storing( p(t) , t )
                else :
                # we make the walker GO to its previous step !! , if alive
                    #print("walker did not die")

                    (nr,nc,nz) = current_position
            ##
            ## UPDATE THE POSITION OF IT
            walker_data[i] = (nr,nc,nz,is_alive)   # is_alive = 1 :: for being alive      is_alive = 0 for being dead

        if current_live_walkers == 0 :
            break
        fraction = current_live_walkers / initial_population_walkers # p(t) =  N_1 / N_o
        p_fraction.append((fraction,t))  # storing( p(t) , t )
        # once we did all walkers in ONE GO , then we update the time
        t = t + delta_t
        bbb = time.time()
        print(f"{jj}|T={t} | F={fraction} || {bbb - aaa}")

    end = time.time()

    print(f"\nOver all time {end - start} ")
    coordinates = analytical_approach(surface_relaxivity , t , RADIUS_IN_MICRO_METERS, interations )
    print(f"\n\n\nSTEP DISTANCE USED IS {step_distance/RESOLUTION} X L")
    surface_relaxation(p_fraction, coordinates)


semi_main()

#main()
