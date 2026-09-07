import matplotlib.pyplot as plt


import numpy as np

import math
from PIL import Image , ImageSequence

from math import *

print(np.log(np.ones(4)))

import time


def get_indexes_of_pore_from_3D_array(array):
    indexes = [] # unknown number
    for z in range(len(array)):
        for r in range(len(array[z])):
            for c in range(len(array[z][r])):
                if array[z][r][c] == 1 :
                    indexes.append((r , c , z))

    return indexes


def analytical_approach(surface_relaxivity=20 , final_time=0.00001 , radius_in_micro_meter=50*400, iterations=100 ):
    # we are going to model the equation  M/M = exp(-3pt/r) , were t=0 --> t=final_t || final_t is aquired AFTER the simulation (when did we stop)
    # this method is used as a benshmark towards the simulation ...
    t = np.linspace(0 , final_time ,iterations )

    y =  ((-3 * t * surface_relaxivity) / radius_in_micro_meter )

    #y = [math.exp(e) for e in y]
    # SINCE WE ARE TAKING THE LOG() OF IT , THEN NOO NEED TO TAKE EXP()


    plt.figure()
    plt.plot(t,y , 'o')
    plt.show()

    return 0





def convert_index_to_mid_point_position(index):
    (r , c , z) = index
    r+=0.5
    c+=0.5
    z+=0.5
    return (r , c , z)



def get_array_from_3D_image(path="/home/saad/Desktop/single_pore.tif"):
    img = Image.open(path)
    full_array = np.array([np.array(page) for page in ImageSequence.Iterator(img) ] )
    print(f"\n3D array has axis {full_array.ndim} with shape {full_array.shape} and has {full_array.size} Voxcels\n")
    return full_array



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
        # but if we land in edge of map , like (121 , 121 , 121 ) --> ( 121.5 , 121.5 , 121.5 )
        position_for_uninitialized_walker = convert_index_to_mid_point_position(index) + (1,)      # so now it is like this (x , y , z , is_alive)

        walker_data.append(position_for_uninitialized_walker)


    return walker_data



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


def update_walker_state(walker_data , new_walker_data  , is_dead , collided ) :


        number_of_dead_in_this_bach=0

        print(f"w:{walker_data}\nNW:{new_walker_data}\n isdead:{is_dead}\ncollided: {collided}")

        all_in_one = np.column_stack((walker_data , new_walker_data  , is_dead , collided))
        # format is as follows :--> [ (r , c , z , is_alive , nr , nc , nz , is_alive , collide , is_dead) ... ]

        groups = [
            all_in_one[:][8] == 0 ,  # this group did not collide into a grain --> update position
            (all_in_one[:][8] == 1) & (all_in_one[:][-1] == True),  # this group collided and died --> mark it as ZERO -> is_alive
            (all_in_one[:][8] == 1) & (all_in_one[:][-1] == False) # collided and did not die --> do not change it position !
        ]

        action = [
            action_no_collision ,
            action_collision_and_dead ,
            action_collision_and_live
        ]


        new_updated_all_in_one = np.piecewise(all_in_one , groups , action)


        print(f"new updated {new_updated_all_in_one}")
        time.sleep(15)






        for i in range(len(is_dead)):

            # if there is no collision
            if collided[i] == 0 :
                walker_data[i] = new_walker_data[i]
            else :
                #if collided and dead
                if is_dead[i] == 1 and not walker_data[i][3] :
                    walker_data[i][3] = 0 # mark as dead
                    number_of_dead_in_this_bach += 1
                else :
                    #if collided and did not die , then do not do anything about the state
                    pass
        return number_of_dead_in_this_bach


def x():
        #################
        #################
        #################
        #INITIALIZATION #
        #################
        #################
        #################
        initial_population_walkers = 100000
        current_live_walkers = initial_population_walkers
        placed_walkers = 0
        # it is much smaller, depending on the resolution of micro-ct
        RESOLUTION = 0.2  # let the lenght of the pixcel be 2 micro-meter

        step_distance = RESOLUTION * 120.2  # s = 0.2 X L .... L is lenght of one pixcel in micro meteres
        fluid_diffusion_coefficient = 2.5e3 # micro meter^2 / second
        surface_relaxivity = 20  # micro-meter/second



        RADIUS = 50 # pixcels

        RADIUS_IN_MICRO_METERS = RADIUS * RESOLUTION





        #is_alive = np.random.default_rng().integers(low = 0 , high = 2 , size = 8)

        #print("this is new demo walker " , np.column_stack((walker_data , is_alive))  )

        #walker_data = np.column_stack((walker_data , is_alive))

        ## given full walker data ##
        interations = 100
        t = 0
        delta_t = calculate_increment_time(step_distance ,fluid_diffusion_coefficient)
        p_fraction = []
        full_array =   get_array_from_3D_image()
        walker_data = np.array(get_initiated_walkers(full_array , initial_population_walkers ))

        for jj in range(interations):
        # EVERY NEW ITERATION , WE MAKE NEW SET OF VALUES FOR THETA & BETA
            theta  = np.ones(len(walker_data))  # [ th , th , th ... th ]
            beta   = np.ones(len(walker_data))  # [ be , be , be ... be
            theta = theta * np.random.default_rng().random( len(walker_data) ) * math.pi * 2
            beta  = beta  * np.random.default_rng().random( len(walker_data) ) * math.pi
        # WE EXTRACT X Y Z SEPRATELY
            r = walker_data[: , 0]
            c = walker_data[: , 1]
            z = walker_data[: , 2]
            is_alive =  walker_data[: , 3]

        # WE COMPUTE NEW VALUES FOR X Y Z
            nr =  r + (step_distance * np.sin(beta) * np.cos(theta))
            nc =  c + (step_distance * np.sin(beta) * np.cos(theta))
            nz =  z + ( step_distance * np.cos(beta))
        #print("parallel computing:-->\n")
        #print(f"\nr-->{r}\nc-->{c}\nz-->{z}")
        #print(f"\nnr-->{nr}\nnc-->{nc}\nnz-->{nz}")
            new_walker_data = np.column_stack((nr , nc , nz , is_alive))
            converted_to_index_new_walker_data = new_walker_data // 1       # from POSITION --> Index
            collided = return_collided_walkers(converted_to_index_new_walker_data , full_array)
            likely = np.ones(len(collided)) * ( 2 * step_distance * surface_relaxivity / (3 *fluid_diffusion_coefficient ) )
            random_number = np.random.default_rng().random(size = len(collided))
            is_dead = random_number > likely ## if random number > likely    then DEAD      # output is [ False , True , False ..... etc ]
            #  walker_data is gonna get updated  and we get HOW MANY died
            number_of_dead_in_this_bach = update_walker_state(walker_data , new_walker_data  , is_dead , collided )
            current_live_walkers = current_live_walkers - number_of_dead_in_this_bach
            fraction = current_live_walkers / initial_population_walkers # p(t) =  N_1 / N_o
            t = t + delta_t
            p_fraction.append((fraction,t))  # storing( p(t) , t )
            print(f"T:{t} || F:{fraction}")


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

    interations = 120


    p_fraction = np.empty(interations)


    for jj in range(interations):


        # we are iterating through each WALKER
        # every thread has this in their private domain, each process will count how many walkers died !

        ## here is start of optimized ##

        theta  = np.ones(len(walker_data))  # [ th , th , th ... th ]
        beta   = np.ones(len(walker_data))  # [ be , be , be ... be
        theta = theta * np.random.default_rng().random( len(walker_data) ) * math.pi * 2
        beta  = beta  * np.random.default_rng().random( len(walker_data) ) * math.pi



        greeks = np.column_stack((theta , beta))

        r = random_walkers[: , 0]
        c = random_walkers[: , 1]
        z = random_walkers[: , 2]





        number_of_walkers_died_in_one_bash = 0
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
                    print("\n_______WALKER DIED_______\n")
                    #print(f"{jj}| new position calculated ({r} {c} {z}) --> ({nr} {nc} {nz})...")
                    #time.sleep(0.005)
                    is_alive = 0 #    here we mark it DEAD

                    fraction = current_live_walkers / initial_population_walkers # p(t) =  N_1 / N_o
                    p_fraction[jj] = (fraction,t)  # storing( p(t) , t )
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
        print(f"{jj}|T={t} | F={fraction}")


    coordinates = analytical_approach(surface_relaxivity , t , RADIUS_IN_MICRO_METERS, interations )
    print(f"\n\n\nSTEP DISTANCE USED IS {step_distance/RESOLUTION} X L")
    surface_relaxation(p_fraction, coordinates)





def calculate_new_position(current_position ,step_distance,theta, beta) :
    s = step_distance
    (r,c,z) =  current_position
    nr = r + s* temp_sin(beta)* temp_cos(theta)
    nc = c + s* temp_sin(beta)* temp_cos(theta)
    nz = z + s* temp_cos(beta)
    return (nr, nc, nz)




def y():

        walker_data = random_walkers = np.random.default_rng().random( (8,4 ))

        step_distance = 0.2 * 0.2
        theta  = np.ones(len(walker_data))  # [ th , th , th ... th ]
        beta   = np.ones(len(walker_data))  # [ be , be , be ... be
        theta = theta * np.random.default_rng().random( len(walker_data) ) * math.pi * 2
        beta  = beta  * np.random.default_rng().random( len(walker_data) ) * math.pi

        sig_beta =  + (step_distance * np.sin(beta) * np.cos(theta))
        print(sig_beta)



def return_collided_walkers(representation_value_of_grain = 0 ) :

            full_array = get_array_from_3D_image()
            leng = 28
            rng = np.random.default_rng()
            r = rng.integers(low=0 , high=60 , size = leng)
            c  = rng.integers(low=0 , high=60 , size = leng)
            z  = rng.integers(low=0 , high=60 , size = leng)
            is_alive =  rng.integers(low=0 , high=2 , size = leng)

            values = full_array[z , r , c] 
            represent_collision_with_this_number = np.max(values)+1
            #print(f"values before\n{values}\n{values[values == representation_value_of_grain]}\n")
            values[values == representation_value_of_grain]  = represent_collision_with_this_number
            values = values // represent_collision_with_this_number
            converted_to_index_new_walker_data = np.floor(np.column_stack((r , c ,z , is_alive)))
            is_collision = np.zeros(len(converted_to_index_new_walker_data))


            #all_z_index_levels_walkers_at = converted_to_index_new_walker_data.T[2]
            #all_c_index_levels_walkers_at = converted_to_index_new_walker_data.T[1]
            ##all_r_index_levels_walkers_at = converted_to_index_new_walker_data.T[0]#

            #result = full_array[all_z_index_levels_walkers_at ,all_r_index_levels_walkers_at ,all_c_index_levels_walkers_at  ]

            #print(f"z:{all_z_index_levels_walkers_at}\nr:{all_r_index_levels_walkers_at}\nc:{all_c_index_levels_walkers_at}\nresult is :{ result }")
            #print(f"testing\n\n1st:{full_array[all_z_index_levels_walkers_at[0]][all_r_index_levels_walkers_at[0]][all_c_index_levels_walkers_at[0]]}\n2nd:{full_array[all_z_index_levels_walkers_at[1]][all_r_index_levels_walkers_at[1]][all_c_index_levels_walkers_at[1]]}\n3nd:{full_array[all_z_index_levels_walkers_at[2]][all_r_index_levels_walkers_at[2]][all_c_index_levels_walkers_at[2]]}\n4nd:{full_array[all_z_index_levels_walkers_at[3]][all_r_index_levels_walkers_at[3]][all_c_index_levels_walkers_at[3]]}#\n5nd:{full_array[all_z_index_levels_walkers_at[4]][all_r_index_levels_walkers_at[4]][all_c_index_levels_walkers_at[4]]}")

             #database = np.array( [
             # [
             #    [2,3,4] ,
             #    [5,6,7] ,
             #    [8,9,10]
             #] ,
             #[
             #    [12,13,14] ,
             #    [15,16,17] ,
             #    [18,19,20]
             #] ,
             #[
             #    [22,23,24] ,
             #    [25,26,27] ,
             #    [28,29,30]
             #]  ] )
             #
             #z = [0  , 2]
             #r = [2  ,  1]
             #c = [1  ,  2]
            #print(f"database::\n{database[z,r,c]}")


            #print(f"converted array is \n{converted_to_index_new_walker_data}\n\n")
            #print(f"z levels:\n{all_z_index_levels_walkers_at}")

            #all_z_levels_walkers_at = full_array[all_z_index_levels_walkers_at]
            #print(f"\n\nwe have {len(all_z_levels_walkers_at)} levels of z\nfirst is {all_z_levels_walkers_at[0].shape}\nsecond is{all_z_levels_walkers_at[1].shape}")


            for i in range(len(converted_to_index_new_walker_data)) :


                (r , c , z , a) = np.abs(converted_to_index_new_walker_data[i])
                #print("index of pixcels is :" , r , c , z , a)
                if r > 120 or c > 120 or z > 120 :
                    is_collision[i] = 1
                    continue

                if full_array[int(z)][int(r)][int(c)] == representation_value_of_grain :
                    is_collision[i] = 1

            print(f"this is iterations:::\n{np.int16(is_collision)}\n\nAnd this is vector::\n{values}")


return_collided_walkers()
# arr = [ (2,3) , (4,5) , (6 , 9 )]

# p_t = np.array(arr)
# print(p_t)
# x = p_t.T[0]

# print(f"x:{x}")

# y = p_t.T[1]

# print(f"y:{y}")
