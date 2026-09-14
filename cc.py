## constants
from math import *

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageSequence

from time import perf_counter

#0.003---t
# it is much smaller, depending on the resolution of micro-ct
RESOLUTION = 0.1  # let the lenght of the pixcel be 2 micro-meter
# 1000 walker
step_distance = RESOLUTION  * 1 # s = 0.2 X L .... L is lenght of one pixcel in micro meteres
fluid_diffusion_coefficient = 2.5e3 # micro meter^2 / second
surface_relaxivity = 20  # micro-meter/second

RADIUS = 50 # pixcels given from doctor

RADIUS_IN_MICRO_METERS = RADIUS * RESOLUTION




def temp_sin(delta):
    return sin(delta)
def temp_cos(delta):
    return cos(delta)
def calculate_new_position(current_position ,step_distance,theta, beta) :
    s = step_distance
    (r,c,z) =  current_position
    nr = r + s* temp_sin(beta)* temp_cos(theta)
    nc = c + s* temp_sin(beta)* temp_sin(theta)
    nz = z + s* temp_cos(beta)
    #print(f"({r},{c},{z})-->({nr},{nc},{nz})")
    #time.sleep(1.1)

    #(10.5,10.9,2.1)-->(10.490318255253307,10.93067667095346,2.076225945279085)
    #(10.5,11.1,2.1)-->(10.487039624241293,11.10400919399578,2.062370822176171)
    #(10.5,11.3,2.1)-->(10.466011761377334,11.314682668605277,2.08486002384558)
    #(10.5,11.5,2.1)-->(10.467272863944928,11.48053968638885,2.0877431341775914)

    return (nr, nc, nz)





def calculate_increment_time(step_distance ,fluid_diffusion_coefficient ):

    delta_t = step_distance * step_distance / (6 * fluid_diffusion_coefficient)
    return delta_t


def calculate_die_probability(step_distance , fluid_diffusion_coefficient, surface_relaxivity):
    likely = 2 * step_distance * surface_relaxivity / (3 *fluid_diffusion_coefficient )
    return likely



# generatte ramodom threshold , and then compare it to
def is_walker_dead(random_number,  likely) :
    if random_number < likely  :
        return True # he is dead
    return False # he is alive
# SURFACE RELAXATION
def surface_relaxation_for_simulation(p_t):
    # this is mimking T_2 .. starts HIGH , ends LOW
    print("\n\nNOW printing the surface relaxation...\n\n\n")

    p_t = np.array(p_t)
    x = p_t.T[1]
    y = np.log(p_t.T[0])
    print_all_coordinations(x,y)
    plt.plot(x,y ,'purple'  , label = "simulation")

def surface_relaxation_for_simulation_no_log(p_t):
    # this is mimking T_2 .. starts HIGH , ends LOW
    print("\n\nNOW printing the surface relaxation...\n\n\n")

    p_t = np.array(p_t)
    x = p_t.T[1]
    y = p_t.T[0]
    print(f"\n\nSimulations\n")
    print_all_coordinations(x,y)
    plt.plot(x,y ,'purple'  , label = "simulation")


def surface_relaxation_for_analytical(coordinates):
    (x , y) = coordinates
    print(f"\n\nAnalytical\n")
    print_all_coordinations(x,y)
    plt.plot(x,y , 'o' , label = "analytical")




def analytical_approach(surface_relaxivity , final_time , radius_in_micro_meter, iterations ):
    # we are going to model the equation  M/M = exp(-3pt/r) , were t=0 --> t=final_t || final_t is aquired AFTER the simulation (when did we stop)
    # this method is used as a benshmark towards the simulation ...
    t = np.linspace(0 , final_time ,iterations )
    y =  np.exp((-3 * t * surface_relaxivity) / radius_in_micro_meter )
    y = np.log(y)# base _ 10
    #print(f"this is the value of y\nnumbers:{len(y)}\nmin:{ np.min(y) }\nmax:{ np.max(y)} " )
    #OVERALL TIME 135.56942486763
    #final time that is going to anaylitical approach is 0.0002666666666666664
    #this is the value of y
    #numbers:100
    #min:-0.0015999999999999522
    #max:0.0
    return (t , y)
def analytical_approach_no_log(surface_relaxivity , final_time , radius_in_micro_meter, iterations ):
        # we are going to model the equation  M/M = exp(-3pt/r) , were t=0 --> t=final_t || final_t is aquired AFTER the simulation (when did we stop)
        # this method is used as a benshmark towards the simulation ...
        t = np.linspace(0 , final_time ,iterations )
        y =  np.exp((-3 * t * surface_relaxivity) / radius_in_micro_meter )

        #print(f"this is the value of y\nnumbers:{len(y)}\nmin:{ np.min(y) }\nmax:{ np.max(y)} " )
        #OVERALL TIME 135.56942486763
        #final time that is going to anaylitical approach is 0.0002666666666666664
        #this is the value of y
        #numbers:100
        #min:-0.0015999999999999522
        #max:0.0
        return (t , y)

def print_all_coordinations(x , y):
    pre_string = "===============\n"
    for i in range(len(x)):
        pre_string += f"  {x[i]:.10f}  ||  {y[i]:.10f}  \n"
    print(pre_string)
def surface_relaxation(p_t , coordinates):
    plt.figure()
    surface_relaxation_for_simulation_no_log(p_t)
    surface_relaxation_for_analytical(coordinates)
    plt.legend(loc = "upper right")
    #plt.tight_layout()


    plt.show(block=False))


    # if we have (3.45 , 2.90 , 12.34) then it is in the voxcel (3,2,12)

#we will store all the position/index of grain voxcels GLOBALLY


def get_array_from_3D_image(path="/home/saad/Desktop/single_pore.tif"):
    img = Image.open(path)
    full_array = np.array([np.array(page) for page in ImageSequence.Iterator(img) ] )
    print(f"\n3D array has axis {full_array.ndim} with shape {full_array.shape} and has {full_array.size} Voxcels\n")
    return full_array

def get_indexes_of_pore_from_3D_array(array , representation_value_of_pore):
    indexes = [] # unknown number   THUS THE INDEX ARE IN FORMAT --> [ (r , c , z) ... ]
    for z in range(len(array)):
        for r in range(len(array[z])):
            for c in range(len(array[z][r])):
                if array[z][r][c] == representation_value_of_pore  :    # pores are 1     grains are 0
                    indexes.append((r , c , z))

    return indexes





def is_index_in_grain(index, full_array , representation_value_of_grain):
    (r , c , z) = index
    if r > 120 or c > 120 or z > 120 :
        return True
    #value = full_array[i][j][k] # (i=z j=x k=y)
    value = full_array[int(z)][int(r)][int(c)]  # cause from walker data (x y z) --> (i j k) ..and array[z][x][y]
    if value == representation_value_of_grain :
        # for logging purpose
        #print("walker is in grain")
        return True
    return False
def convert_position_to_index_unoptimized_main(position , resolution):
    ( x,y,z) = position
    return (x//resolution , y // resolution , z // resolution )

def is_collision(position, full_array,resolution , representation_value_of_grain):
    index = convert_position_to_index_unoptimized_main(position,resolution )
    return is_index_in_grain(index , full_array , representation_value_of_grain)
# from 0 to 2_pi
import math


def random_theta(rng):
    inclusive_high = np.nextafter(math.pi * 2 , np.inf)
    theta = rng.uniform(0 ,   inclusive_high )
    return theta

# from 0 to pie
def random_beta(rng):
    inclusive_high = np.nextafter(math.pi * 1, np.inf)
    beta = rng.uniform(0 , inclusive_high )
    return beta

import time


def convert_index_to_mid_point_position(index , resolution):
    (r , c , z) = index

    pr = (r) * resolution + (0.5 * resolution)
    pc = (c) * resolution + (0.5 * resolution)
    pz = (z) * resolution + (0.5 * resolution)

    if r == 0 or c == 0 or z == 0 :
        if r == 0 :
            pr = (0.5 * resolution)
        if c == 0 :
            pc = (0.5 * resolution)
        if z == 0 :
            pz = (0.5 * resolution)

    #print(f"\nconverting index to position , for a walker :: ({r},{c},{z})-->  ({pr} , {pc} , {pz}) ")
    #time.sleep(10)
    return (pr , pc , pz)

def get_initiated_walkers(full_array,representation_value_of_pore, resolution, number_of_walkers_initiated = 1000):
    pores_indexes = get_indexes_of_pore_from_3D_array(full_array, representation_value_of_pore)
    number_of_pores = len(pores_indexes)
    #print(f"number of walkers should me {number_of_pores}")
    #time.sleep(10)

    #dt = np.dtype( [('x' , np.float32 ),('y' , np.float32 ),('z' , np.float32 ),('is_alive' , np.int16)] )
    # np.empty(number_of_walkers_initiated , dtype = dt)
    walker_data = [] # LATER PRE-ALLOCATE !!



    # number of walkers should me 523,305
    # we converted (60, 59, 9) to (60.5, 59.5, 9.5, 1) for walker[0]
    # we converted (51, 55, 10) to (51.5, 55.5, 10.5, 1) for walker[1]
    # we converted (51, 56, 10) to (51.5, 56.5, 10.5, 1) for walker[2]
    # we converted (51, 57, 10) to (51.5, 57.5, 10.5, 1) for walker[3]
    # we converted (51, 58, 10) to (51.5, 58.5, 10.5, 1) for walker[4]
    # we converted (51, 59, 10) to (51.5, 59.5, 10.5, 1) for walker[5]
    # we converted (51, 60, 10) to (51.5, 60.5, 10.5, 1) for walker[6]
    # we converted (51, 61, 10) to (51.5, 61.5, 10.5, 1) for walker[7]
    #
    # we loop through every PORE , we put ONE walker it in ...
    for i in range(number_of_walkers_initiated):

        index = pores_indexes[i] # we will get the first pore's pixcel , then second ..
        # 2 we calculate the MIDDLE POSITION OF THIS PORE (4 , 5 , 8) --> (4.5 , 5.5 , 8.5)
        position_for_uninitialized_walker = convert_index_to_mid_point_position(index , resolution) + (1,)  # so now it is like this (x , y , z , is_alive)
     #   print(f"we converted {index} to {position_for_uninitialized_walker} for walker[{i}]")
        walker_data.append(position_for_uninitialized_walker)



    return walker_data


def return_collided_walkers(converted_to_index_new_walker_data , full_array , representation_value_of_grain = 0 ) :
    converted_to_index_new_walker_data=np.int32(converted_to_index_new_walker_data)
    all_z_index_levels_walkers_at = converted_to_index_new_walker_data.T[2]
    all_c_index_levels_walkers_at = converted_to_index_new_walker_data.T[1]
    all_r_index_levels_walkers_at = converted_to_index_new_walker_data.T[0]
    #print(f"type{all_z_index_levels_walkers_at[0]}")
    result = full_array[ all_z_index_levels_walkers_at , all_r_index_levels_walkers_at , all_c_index_levels_walkers_at  ]
    # now we have the VALUES of the pixcels .. there might be 1 , 2, 3 ..etc as values
    # # we do not know which is grain , unless TOLD


    # first detech the GRAIN , then change their value --> unique = MAX() + 1
    represent_collision_with_this_number = np.max(result)+1
    # we replace , each GRAIN value , with this UNIQUE number
    result[result == representation_value_of_grain]  = represent_collision_with_this_number
    result = result // represent_collision_with_this_number  # all other values are 0 , cause they are smaller , except GRAINS --> 1 which is COLLISION
    return result


def convert_position_to_index(walker_data , resolution):
    #the format of walker_data is [[nr , nc , nz , is_alive] () () ..()]

    # if i did    pr =  (r-1) * resolution
    # then to get r --> (pr/res) + 1
    # now let me check

    # get a copy
    indexed_walker_data = np.copy(walker_data)
    # operate on all
    indexed_walker_data = np.floor(indexed_walker_data / resolution).astype(np.int32)
    #print(f"\n\nWalker Data is converted to represent the INDEX of pixcels they represent \n{indexed_walker_data[:5]}")
    #time.sleep(4)

      #Walker Data is converted to represent the INDEX of pixcels they represent
      #[[60 59  9]
      # [51 55 10]
      # [51 56 11]
      # [51 57 10]
      # [51 58 10]]

    return indexed_walker_data

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
    walker_data = np.array(get_initiated_walkers(full_array ))

    seed = 99
    rng = np.random.default_rng(seed)


    # now start the simulation

    #(1) first we increment the time
    global current_live_walkers

    global p_fraction
    current_live_walkers = len(walker_data)
    global initial_population_walkers
    initial_population_walkers = len(walker_data)
    t = 0

    delta_t = calculate_increment_time(step_distance ,fluid_diffusion_coefficient)

    interations = 7500
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
            theta =  np.random.default_rng(seed).random( len(walker_data) ) * math.pi * 2
            beta  =  np.random.default_rng(seed).random( len(walker_data) ) * math.pi
        # WE EXTRACT X Y Z SEPRATELY

            nr =  r + (step_distance * np.sin(beta) * np.cos(theta))
            nc =  c + (step_distance * np.sin(beta) * np.cos(theta))
            nz =  z + ( step_distance * np.cos(beta) )
            number_of_negitive_values = len(nr[nr < 0 ]) + len(nc[nc < 0 ]) + len(nz[nz < 0 ])
            new_walker_data = np.column_stack((nr , nc , nz , is_alive))

            converted_to_index_new_walker_data = convert_position_to_index( new_walker_data , RESOLUTION  )  # from POSITION --> Index
            converted_to_index_new_walker_data = np.column_stack((converted_to_index_new_walker_data , is_alive))
            collided = return_collided_walkers(converted_to_index_new_walker_data, full_array,representation_value_of_grain)

            # should be before the two for loops cause it is a constant, but should vary
            # IF IT IN TOUCHING OR BEHOND THE GRAIN, WE CALCULATE TEH LIKELY
               # print(f"walker has HIT a grain , now calculating if it is dead.......?")
                #time.sleep(0.5)
            likely = np.ones(len(collided)) * ( 2 * step_distance * surface_relaxivity / (3 *fluid_diffusion_coefficient ) )
            random_number = np.random.default_rng(seed).random(size = len(collided))
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
            print(f"{jj}|| T:{t} || F:{fraction} || {bbb - aaa} ||negitive : {number_of_negitive_values}")

            # is_alive = 1 :: for being alive      is_alive = 0 for being dead




    end = time.time()

    print(f"\n\n\nOVERALL TIME {end - start} ")
    print(f"final time that is going to anaylitical approach is {t}")
    coordinates = analytical_approach(surface_relaxivity , t , RADIUS_IN_MICRO_METERS, interations )

    print(f"\n\n\nSTEP DISTANCE USED IS {step_distance/RESOLUTION} X L")
    surface_relaxation(p_fraction, coordinates)


# the walkers that have is_alive = 1 , are kept .. the others are gone !
def garbage_clear_for_walker_data(walker_data):
    return walker_data[walker_data.T[-1] == 1 ]

def semi_semi_main( boost_number = 0.18  , iterations = 600):
    # get the image , and store it in np_array
    full_array =  get_array_from_3D_image()
    boundries = len(full_array)

    representation_value_of_pore = 1
    representation_value_of_grain = 0
     # now it looks like this ( r , c , z , is_alive=1 )
    walker_data = np.array(get_initiated_walkers(full_array,representation_value_of_pore,RESOLUTION ))

    seed = 99

    rng = np.random.default_rng(seed)


    global current_live_walkers

    global p_fraction
    current_live_walkers = len(walker_data)
    global initial_population_walkers
    initial_population_walkers = len(walker_data)
    t = 0

    delta_t = calculate_increment_time(step_distance ,fluid_diffusion_coefficient)
    # should be 900 iterations

    p_fraction = []


    start = time.time()
    for jj in range(iterations):
        # every 500 iterations , we should erase all the zombie walkers
        if jj % 500 == 0 :
           walker_data =  garbage_clear_for_walker_data(walker_data)

        aaa = time.time()
        length_of_walker_data_array = len(walker_data)
        r = walker_data[:, 0]
        c = walker_data[:, 1]
        z = walker_data[:, 2]
        is_alive = walker_data[:, 3]

        # theta  = np.ones(len(walker_data))  # [ th , th , th ... th ]
        # beta   = np.ones(len(walker_data))  # [ be , be , be ... be
        theta = rng.random(length_of_walker_data_array) * math.pi * 2
        beta = rng.random(length_of_walker_data_array) * math.pi
        # WE EXTRACT X Y Z SEPRATELY

        nr = r + (step_distance * np.sin(beta) * np.cos(theta))
        nc = c + (step_distance * np.sin(beta) * np.sin(theta))
        nz = z + (step_distance * np.cos(beta))

        # we want to check here , how many of walkers are outside boundries

        # here we are just CHECKING if we generated NEGITIVE positions !
        #number_of_negitive_values = len(nr[nr < 0]) + len(nc[nc < 0]) + len(nz[nz < 0])

        new_walker_data = np.column_stack((nr, nc, nz, is_alive))


        converted_to_index_new_walker_data = convert_position_to_index( new_walker_data , RESOLUTION  )  # from POSITION --> Index
        converted_to_index_new_walker_data = np.column_stack((converted_to_index_new_walker_data , is_alive))
        collided = return_collided_walkers(converted_to_index_new_walker_data, full_array,representation_value_of_grain)

        lenght_of_collided_array = len(collided)
        # IF IT IN TOUCHING OR BEHOND THE GRAIN, WE CALCULATE TEH LIKELY
        likely = np.ones(lenght_of_collided_array) * (2 * step_distance * surface_relaxivity / (3 * fluid_diffusion_coefficient))

        # from 0.17 ---> 0.25
        # boost_number = 0.18 # we decrease the amplitude of the numbers in random_numbers
        random_number = rng.random(size=lenght_of_collided_array) * boost_number

        # random_number > likely KILLS TOO FAST !
        #is_dead = (random_number >  likely  )

        #a bit too slow it gave F:->0.9800 and it is suppose to be F:0.930
        is_dead =  (random_number  < likely  )
        # we can imaging likely as being a Bar high in the air
        # if you manage to jump over it , then you live
        # if you can't , you die
        # now boost  , makes you a boost as a help to jump,
        # the bigger the boost , the more likely you will live
        # hence LOW BOOST --> more death --> more magnetization --> more steap incline line



        ## O(n) ... we need to loop Through the entire set , to update ..


        ## IF THERE IS NO COLLISION , THEN r <- nr
        r[collided == 0] = nr[collided == 0]
        c[collided == 0] = nc[collided == 0]
        z[collided == 0] = nz[collided == 0]

        ## IF THERE IS COLLISION AND DEAD , THEN is_alive=0
        before_killing = np.sum(is_alive)
        is_alive[(is_dead == 1) & (collided == 1)] = 0
        after_killing = np.sum(is_alive)

        number_of_dead_in_this_bach = before_killing - after_killing
        current_live_walkers = current_live_walkers - number_of_dead_in_this_bach

        fraction = (current_live_walkers / initial_population_walkers )  # p(t) =  N_1 / N_o
        true_fraction = math.exp((-3 * t * surface_relaxivity) / RADIUS_IN_MICRO_METERS )
        t = t + delta_t


        if current_live_walkers == 0:
            break

        p_fraction.append((fraction, t))
        ## [() () () () ()]
        bbb = time.time()
        print(f"{jj}|| T:{t} || F:{fraction} || {bbb - aaa} || diff : {true_fraction - fraction}")

            # is_alive = 1 :: for being alive      is_alive = 0 for being dead




    end = time.time()

    print(f"\nOVERALL TIME {end - start} ")
    time.sleep(2)
    print(f"t -> anaylyical ::is::--> {t}")
    time.sleep(1.3)
    coordinates = analytical_approach(surface_relaxivity , t , RADIUS_IN_MICRO_METERS, iterations )
    # (t , y)
    print(f"\nSTEP DISTANCE USED IS {step_distance/RESOLUTION} X L")
    time.sleep(1.5)
    return (p_fraction, coordinates , end - start)



def un_optimized_main():
    # initialize walkers then store them
    full_array =  get_array_from_3D_image()

     # now it looks like this ( X Y Z w )    w = 0 , 1
    representation_value_of_pore = 1
    representation_value_of_grain = 0
    walker_data = np.array(get_initiated_walkers(full_array ,representation_value_of_pore , RESOLUTION))
    # [(r,c,z , is_alive)]
    current_live_walkers = len(walker_data)

    initial_population_walkers = current_live_walkers

    t = 0

    delta_t = calculate_increment_time(step_distance ,fluid_diffusion_coefficient)

    interations = 50
    p_fraction = [(1.0 , 0.0)]

    start = time.time()
    seed = 99
    rng = np.random.default_rng(seed)
    print(f"{interations} iterations")
    likely = calculate_die_probability(step_distance , fluid_diffusion_coefficient, surface_relaxivity)
    for jj in range(interations):



        aaa = time.time()
        for i in  range(len(walker_data)):
            # if this walker DIED , then skip iteration



            if walker_data[i][-1] :
                (r , c , z , is_alive) = walker_data[i]
                current_position = (r , c, z)
                theta = random_theta(rng)
                beta = random_beta(rng)
                (nr,nc,nz) = calculate_new_position( current_position , step_distance , theta , beta )
                if is_collision((nr, nc , nz) , full_array , RESOLUTION , representation_value_of_grain):
                    random_number = rng.random()  #np.random.random()
                    if is_walker_dead(random_number , likely):
                        current_live_walkers = current_live_walkers - 1
                        is_alive = 0 #    here we mark it DEAD
                    else :
                    # if not dead , go prev position
                        (nr,nc,nz) = current_position
            ##
            ## UPDATE THE POSITION OF IT
                walker_data[i] = (nr,nc,nz,is_alive)   # is_alive = 1 :: for being alive      is_alive = 0 for being dead
                if current_live_walkers == 0 :
                    break

        # FINISHED ALL THE WALKERS ...
        fraction = current_live_walkers / initial_population_walkers # p(t) =  N_1 / N_o
        true_fraction = math.exp((-3 * t * surface_relaxivity) / RADIUS_IN_MICRO_METERS )
        p_fraction.append((fraction,t))  # storing( p(t) , t )
        # once we did all walkers in ONE GO , then we update the time
        t = t + delta_t
        bbb = time.time()
        print(f"{jj}|T={t} | F={fraction} || {bbb - aaa} || diff : {true_fraction - fraction} ")

    end = time.time()

    print(f"\nOver all time {end - start} ")
    coordinates = analytical_approach_no_log(surface_relaxivity , t , RADIUS_IN_MICRO_METERS, interations )
    print(f"\nSTEP DISTANCE USED IS {step_distance/RESOLUTION} X L")
    surface_relaxation(p_fraction, coordinates)



un_optimized_main()
# still have not tryied this
import logging
import os

def main():
    final_time = 0.000005
    factor = 1
    while(final_time < 1) :
        boost_number_list = [ 1 ]
        iterations = 600 * factor
        for boost_number in boost_number_list :
               (p_fraction, coordinates , total_time ) = semi_semi_main( boost_number   , iterations)
               final_time = p_fraction[-1][1]
               final_ratio = p_fraction[-1][0]

               ### here we should retreive the image , name it , store it
               print(f"boost Number :{boost_number}\nIterations:{iterations}\nF_t:{final_time}\nfinal_ration:{final_ratio}")
               time.sleep(1)
               try:
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

               except Exception as exc:
                   logging.exception("Exception occurred")

        factor = factor * 2  ### next time multiply by 2,3,4,5...?


#main()

    # store the OVERAL time
    # the T_final
    # F final


#semi_main() ;print("semi main()")



#print("main()")
