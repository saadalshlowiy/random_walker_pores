## constants
from ast import walk
from math import *

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageSequence

from time import perf_counter


# it is much smaller, depending on the resolution of micro-ct
RESOLUTION = 0.2  # let the lenght of the pixcel be 2 micro-meter
step_distance = RESOLUTION  * 0.2  # s = 0.2 X L .... L is lenght of one pixcel in micro meteres
fluid_diffusion_coefficient = 2.5e3 # micro meter^2 / second
surface_relaxivity = 20  # micro-meter/second
RADIUS = 50 # pixcels given from doctor
RADIUS_IN_MICRO_METERS = RADIUS * RESOLUTION


def calculate_increment_time(step_distance ,fluid_diffusion_coefficient ):

    delta_t = (step_distance * step_distance) / (6 * fluid_diffusion_coefficient)
    #print(f"\ncalcuting Delta T --> {delta_t}")
    return delta_t

def surface_relaxation_for_simulation(p_t):
    # this is mimking T_2 .. starts HIGH , ends LOW
    #print("\n\nNOW #printing the surface relaxation...\n\n\n")
    p_t = np.array(p_t)
    x = p_t.T[1]    # x is TIME ... p_t (ratio , time)
    y = np.log(p_t.T[0])
    plt.plot(x,y ,'purple'  , label = "simulation")
def surface_relaxation_for_analytical(coordinates):
    (x , y) = coordinates
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
def surface_relaxation(p_t , coordinates):
    plt.figure()
    surface_relaxation_for_simulation(p_t)
    surface_relaxation_for_analytical(coordinates)
    plt.legend(loc = "upper right")
    #plt.tight_layout()





    # if we have (3.45 , 2.90 , 12.34) then it is in the voxcel (3,2,12)

#we will store all the position/index of grain voxcels GLOBALLY


def get_array_from_3D_image(path="/home/saad/Desktop/single_pore.tif"):
    img = Image.open(path)
    full_array = np.array([np.array(page) for page in ImageSequence.Iterator(img) ] )
    #print(f"\n3D array has axis {full_array.ndim} with shape {full_array.shape} and has {full_array.size} Voxcels\n")
    return full_array

def get_indexes_of_pore_from_3D_array(array , representation_value_of_pore):
    indexes = [] # unknown number   THUS THE INDEX ARE IN FORMAT --> [ (r , c , z) ... ]
    lz = len(array) ; lr = len(array[0]) ; lc = len(array[0][0])
    for z in range(lz):
        for r in range(lr):
            for c in range(lc):
                if array[z][r][c] == representation_value_of_pore  :    # pores are 1     grains are 0
                    indexes.append((r , c , z))

    return indexes
import time

## the index of a pore MAY be something like (60, 59, 9) ..BUT !!
# 60 does not represent the distance from an AXIS
# it repreesnts the number ... I am the 60th pixcel/pore/apple !
# we need to represent HOW far it is !
# if i am in the 3 poxcel , then i have 2 previous poxcels ! --> 2 X L --> then add 0.5L since we want it to be in HALF
def convert_index_to_mid_point_position(index , resolution):
    (r , c , z) = index

    # instead of 12.3000000000000000004 --> round() ---> 12.300
    pr = round((r) * resolution + (0.5 * resolution) , 3)
    pc = round((c) * resolution + (0.5 * resolution) , 3)
    pz = round((z) * resolution + (0.5 * resolution) , 3)


    #print(f"\nconverting index to position , for a walker :: ({r},{c},{z})-->  ({pr} , {pc} , {pz}) ")
    #time.sleep(10)
    return (pr , pc , pz)


def get_initiated_walkers(full_array,representation_value_of_pore,resolution):
    pores_indexes = get_indexes_of_pore_from_3D_array(full_array, representation_value_of_pore)
    number_of_pores = len(pores_indexes)
    walker_data = []
    #print(f"number of walkers should me {number_of_pores}")
    #time.sleep(10)

    # number of walkers should me 523,305
    #converting index to position , for a walker :: (69,61,108)-->  (13.700000000000001 , 12.1 , 21.500000000000004)
    #converting index to position , for a walker :: (69,62,108)-->  (13.700000000000001 , 12.3 , 21.500000000000004)
    # we loop through every PORE , we put ONE walker it in ...
    for i in range(number_of_pores):
        index = pores_indexes[i] # we will get the first pore's pixcel , then second ..
        # 2 we calculate the MIDDLE POSITION OF THIS PORE (4 , 5 , 8) --> (4.5 , 5.5 , 8.5)
        position_for_uninitialized_walker = convert_index_to_mid_point_position(index , resolution) + (1,)  # so now it is like this (x , y , z , is_alive=1)
        #print(f"we converted {index} to {position_for_uninitialized_walker} for walker[{i}]")
        walker_data.append(position_for_uninitialized_walker)
    return walker_data



# this method did not work , cause it produced a stairs in the simulation !!
def get_initiated_walkers_random(full_array , representation_value_of_pore, resolution , number_of_walkers = 10000):
    pores_indexes = get_indexes_of_pore_from_3D_array(full_array, representation_value_of_pore)
    walker_data = []
    number_of_pores = len(pores_indexes)
    for i in range(number_of_walkers):
        # randomly pick pore index
        index = pores_indexes[np.random.default_rng().integers(low = 0 , high = number_of_pores ) ]
        position_for_uninitialized_walker = convert_index_to_mid_point_position(index , resolution) + (1,)
        walker_data.append(position_for_uninitialized_walker)
        #print(f"walker Born ({index}) --> ({position_for_uninitialized_walker})")
        #time.sleep(0.5)
    return walker_data


def return_collided_walkers(converted_to_index_new_walker_data , full_array , representation_value_of_grain = 1 ) :
    converted_to_index_new_walker_data= np.int32(converted_to_index_new_walker_data)
    #print("\n Representing Walkers index's ... we will take only 5 walkers and show them\n<===============>\n")
    #print(f"\nall walkers here are converted to Index of their voxcels--->{converted_to_index_new_walker_data[0:5]}..\n...")
    all_z_index_levels_walkers_at = converted_to_index_new_walker_data.T[2]
    all_c_index_levels_walkers_at = converted_to_index_new_walker_data.T[1]
    all_r_index_levels_walkers_at = converted_to_index_new_walker_data.T[0]
    #print(f"All Z --> {all_z_index_levels_walkers_at[:5]}..\nAll R --> {all_r_index_levels_walkers_at[:5]}..\nAll C --> {all_c_index_levels_walkers_at[:5]}..\n")
    #print(f"type{all_z_index_levels_walkers_at[0]}")
    #time.sleep(10)
    result = full_array[ all_z_index_levels_walkers_at , all_r_index_levels_walkers_at , all_c_index_levels_walkers_at  ]
    # now we have the VALUES of the pixcels .. there might be 1 , 2, 3 ..etc as values
    # # we do not know which is grain , unless TOLD
    #print(f"\nThese voxcels have grains -->\n")
    pre = len(result)

    # first detech the GRAIN , then change their value --> unique = MAX() + 1
    represent_collision_with_this_number = np.max(result)+1
    # we replace , each GRAIN value , with this UNIQUE number
    result[result == representation_value_of_grain]  = represent_collision_with_this_number
    result = result // represent_collision_with_this_number  # all other values are 0 , cause they are smaller , except GRAINS --> 1 which is COLLISION

    #after = len(result)


    #xxx=len(converted_to_index_new_walker_data)
    #for i in range(xxx) :
    #    #print(f"{converted_to_index_new_walker_data[i]}")
    #    #time.sleep(14)
    #    (r , c , z , a ) = converted_to_index_new_walker_data[i]
    #    #print("index of pixcels is :" , r , c , z , a)
    #    if r > 120 or c > 120 or z > 120 :
        #        continue
        #    if full_array[int(z)][int(r)][int(c)] == representation_value_of_grain :
            #        if result[i] != 1 :
                #            print(f"incorrect calculation of is_collide")
                #            raise(Exception)
    return result


# the walkers that have is_alive = 1 , are kept .. the others are gone !
def garbage_clear_for_walker_data(walker_data):
    return walker_data[walker_data.T[-1] == 1 ]

def convert_position_to_index(walker_data , resolution):
    #the format of walker_data is [[nr , nc , nz , is_alive] () () ..()]

    # if i did    pr =  (r-1) * resolution
    # then to get r --> (pr/res) + 1
    # now let me check

    # get a copy
    indexed_walker_data = np.copy(walker_data)
    # operate on all
    indexed_walker_data = np.floor(
        indexed_walker_data / resolution
    ).astype(np.int32)

    #print(f"\n\nWalker Data is converted to represent the INDEX of pixcels they represent \n{indexed_walker_data[:5]}")
    #time.sleep(4)

      #Walker Data is converted to represent the INDEX of pixcels they represent
      #[[60 59  9]
      # [51 55 10]
      # [51 56 11]
      # [51 57 10]
      # [51 58 10]]

    return indexed_walker_data
import math
#Smooth sphere center: 12.1 11.9 11.9
#Smooth sphere radius: 10.0
#Voxel-center walkers outside smooth sphere: 35
def calculate_discrete_pore_geometry(
    full_array,
    resolution,
    pore_value=1
):
    pore = full_array == pore_value

    # Pore volume
    pore_voxels = np.count_nonzero(pore)
    volume = pore_voxels * resolution**3

    # Count pore faces adjacent to grain
    exposed_faces = 0

    # Internal boundaries in all three dimensions
    # along the z axis, we detect difference in value between z=1 -- z=0   then z=2 -- z=1 ..etc
    exposed_faces += np.sum(pore[1:, :, :] != pore[:-1, :, :])
    # along the row axis
    exposed_faces += np.sum(pore[:, 1:, :] != pore[:, :-1, :])
    # along the column axis
    exposed_faces += np.sum(pore[:, :, 1:] != pore[:, :, :-1])



    surface_area = exposed_faces * resolution**2
    surface_to_volume = surface_area / volume

    equivalent_radius = 3.0 / surface_to_volume
    #print(f"number of pore Voxcels is {pore_voxels}\n\nNumber Of Exposed Faces to Grain Surface is {exposed_faces}\n\ncalculated Vol {volume}\n\nSurface Area of image sphere is {surface_area}\n\nS/V is {surface_to_volume}\n")
    #print(f"number of pore Voxcels is {pore_voxels}") ; time.sleep(1.9)
    #print(f"nNumber Of Exposed Faces to Grain Surface is {exposed_faces}\n")
    print(f"\n=======\ncalculated Vol {volume}\n\n=======\nSurface Area of image sphere is {surface_area}\n\n=======\nS/V is {surface_to_volume}\n")
    return({
        "pore_voxels": pore_voxels,
        "exposed_faces": exposed_faces,
        "volume": volume,
        "surface_area": surface_area,
        "surface_to_volume": surface_to_volume,
        "equivalent_radius": equivalent_radius
    })
def semi_semi_main(
        boost_number=1,
        iterations=600,
        original_method_detecting_grain=0
    ):

        # get the image , and store it in np_array
        full_array = get_array_from_3D_image()
        print(np.unique(full_array))

        boundries = len(full_array)

        representation_value_of_pore = 1
        representation_value_of_grain = 0
        print(f"here is the geometry of the given image ..."); time.sleep(2)
        calculate_discrete_pore_geometry(full_array , RESOLUTION , representation_value_of_pore)
        time.sleep(9)
        # now it looks like this ( r , c , z , is_alive=1 )
        walker_data = np.array(
            get_initiated_walkers(
                full_array,
                representation_value_of_pore,
                RESOLUTION
            ),
            dtype=float
        )

        pore_z, pore_r, pore_c = np.where(
            full_array == representation_value_of_pore
        )

        # find the center of the sphere , so that we can , reference it
        # to see which walker is outside of the sphere


        # since the IMAGE has the sphere in the middle exactly
        # we can use median , if not than mean

        center_r = np.median(
            (pore_r + 0.5) * RESOLUTION
        )

        center_c = np.median(
            (pore_c + 0.5) * RESOLUTION
        )

        center_z = np.median(
            (pore_z + 0.5) * RESOLUTION
        )

        sphere_radius = RADIUS_IN_MICRO_METERS

        if original_method_detecting_grain == 0:

            initial_distance_squared = (
                (walker_data[:, 0] - center_r) ** 2
                + (walker_data[:, 1] - center_c) ** 2
                + (walker_data[:, 2] - center_z) ** 2
            )

            inside_sphere = (
                initial_distance_squared
                <= sphere_radius ** 2
            )
            print(f"the number of walkers outside the sphere is {len(walker_data) - len(walker_data[inside_sphere])}")
            walker_data = walker_data[inside_sphere]
            time.sleep(4)
            print("using radius method")

        else:

            print("using original grain detection method")

        seed = 99
        rng = np.random.default_rng(seed)

        global current_live_walkers
        global p_fraction
        global initial_population_walkers

        current_live_walkers = len(walker_data)
        initial_population_walkers = len(walker_data)

        t = 0

        delta_t = calculate_increment_time(
            step_distance,
            fluid_diffusion_coefficient
        )

        # should be 900 iterations
        p_fraction = [(1.0, 0.0)]

        start = time.time()

        size_z, size_r, size_c = full_array.shape

        for jj in range(iterations):

            # every 500 iterations , we should erase all the zombie walkers
            # if jj % 500 == 0:
            #     walker_data = garbage_clear_for_walker_data(walker_data)

            aaa = time.time()

            is_alive = walker_data[:, 3]

            live_indices = np.flatnonzero(
                is_alive == 1
            )

            if len(live_indices) == 0:
                break

            length_of_walker_data_array = len(
                live_indices
            )

            # this is how to extract all r ! we confirmed it
            r = walker_data[live_indices, 0]
            c = walker_data[live_indices, 1]
            z = walker_data[live_indices, 2]

            theta = rng.uniform(
                0.0,
                2.0 * np.pi,
                length_of_walker_data_array
            )

            cos_beta = rng.uniform(
                -1.0,
                1.0,
                length_of_walker_data_array
            )

            sin_beta = np.sqrt(
                1.0 - cos_beta ** 2
            )

            # WE EXTRACT X Y Z SEPRATELY
            nr = (
                r
                + step_distance
                * sin_beta
                * np.cos(theta)
            )

            nc = (
                c
                + step_distance
                * sin_beta
                * np.sin(theta)
            )

            nz = (
                z
                + step_distance
                * cos_beta
            )

            if original_method_detecting_grain == 1:

                new_walker_data = np.column_stack(
                    (nr, nc, nz)
                )

                converted_to_index_new_walker_data = (
                    convert_position_to_index(
                        new_walker_data,
                        RESOLUTION
                    )
                )

                new_r_index = (
                    converted_to_index_new_walker_data[:, 0]
                )

                new_c_index = (
                    converted_to_index_new_walker_data[:, 1]
                )

                new_z_index = (
                    converted_to_index_new_walker_data[:, 2]
                )

                outside_boundries = (
                    (new_r_index < 0)
                    | (new_r_index >= size_r)
                    | (new_c_index < 0)
                    | (new_c_index >= size_c)
                    | (new_z_index < 0)
                    | (new_z_index >= size_z)
                )

                collided = np.zeros(
                    length_of_walker_data_array,
                    dtype=bool
                )

                collided[outside_boundries] = True

                inside_boundries = ~outside_boundries

                collided[inside_boundries] = (
                    full_array[
                        new_z_index[inside_boundries],
                        new_r_index[inside_boundries],
                        new_c_index[inside_boundries]
                    ]
                    == representation_value_of_grain
                )

            else:

                distance_from_center_squared = (
                    (nr - center_r) ** 2
                    + (nc - center_c) ** 2
                    + (nz - center_z) ** 2
                )

                collided = (
                    distance_from_center_squared
                    > sphere_radius ** 2
                )

            lenght_of_collided_array = len(collided)

            # IF IT IN TOUCHING OR BEHOND THE GRAIN, WE CALCULATE TEH LIKELY
            likely = np.ones(
                lenght_of_collided_array
            ) * (
                (
                    2
                    * step_distance
                    * surface_relaxivity
                )
                / (
                    3
                    * fluid_diffusion_coefficient
                )
            )

            # from 0.17 ---> 0.25
            # boost_number = 0.18
            # we decrease the amplitude of the numbers in random_numbers

            random_number = rng.random(
                size=lenght_of_collided_array
            )

            # random_number > likely KILLS TOO FAST !
            # is_dead = (random_number > likely)

            # a bit too slow it gave F:->0.9800 and it is suppose to be F:0.930
            is_dead = (
                random_number < likely
            )

            # IF THERE IS NO COLLISION , THEN r <- nr
            accepted_move = collided == 0

            accepted_global_indices = live_indices[
                accepted_move
            ]

            walker_data[
                accepted_global_indices,
                0
            ] = nr[accepted_move]

            walker_data[
                accepted_global_indices,
                1
            ] = nc[accepted_move]

            walker_data[
                accepted_global_indices,
                2
            ] = nz[accepted_move]

            # IF THERE IS COLLISION AND DEAD , THEN is_alive=0
            killed_this_step = (
                (is_dead == 1)
                & (collided == 1)
            )

            killed_global_indices = live_indices[
                killed_this_step
            ]

            before_killing = np.sum(
                walker_data[:, 3]
            )

            walker_data[
                killed_global_indices,
                3
            ] = 0

            after_killing = np.sum(
                walker_data[:, 3]
            )

            number_of_dead_in_this_bach = (
                before_killing
                - after_killing
            )

            current_live_walkers = (
                current_live_walkers
                - number_of_dead_in_this_bach
            )

            # p(t) = N_1 / N_o
            fraction = (
                current_live_walkers
                / initial_population_walkers
            )

            t = t + delta_t

            true_fraction = math.exp(
                (
                    -3
                    * t
                    * surface_relaxivity
                )
                / RADIUS_IN_MICRO_METERS
            )

            if current_live_walkers == 0:
                break

            p_fraction.append(
                (fraction, t)
            )

            bbb = time.time()

            print(
                f"{jj}"
                f"|| METHOD:{original_method_detecting_grain}"
                f"|| T:{t}"
                f"|| F:{fraction}"
                f"|| {bbb - aaa}"
                f"|| diff:{true_fraction - fraction}"
                f"|| collisions:{np.sum(collided)}"
                f"|| deaths:{number_of_dead_in_this_bach}"
            )

            # is_alive = 1 :: for being alive
            # is_alive = 0 for being dead

        end = time.time()

        print(
            f"\nOVERALL TIME {end - start}"
        )

        print(
            f"t -> anaylyical ::is::--> {t}"
        )

        coordinates = analytical_approach(
            surface_relaxivity,
            t,
            RADIUS_IN_MICRO_METERS,
            len(p_fraction)
        )

        # (t , y)
        print(
            f"\nSTEP DISTANCE USED IS "
            f"{step_distance / RESOLUTION} X L"
        )

        return (
            p_fraction,
            coordinates,
            end - start
        )
#this is the stairs spher
# the killings were too high , which means only two things ::
# A. Killing probability is too high

#or

#B. Walkers encounter the surface too frequently
# so I looked into this more
# I studied the digital sphera ... cause the analytical was wanting a perfect
# smooth sphere , meaning the surface is sooo smooth .. and thus
# the equation is built on that
# but i wanted to know if digital sphere , had the same smooth surface ,
# that i can find with surface area , and compare
# surface area of digital , with theoritcal surface area that the equation wants
# The digital sphere had a larger surface area , with S = {} and S/V =
# the theritical sphere has S = 1256.64 and S/V = 0.30
#
# To future more make sure that the logic used is accurate , we
# are going to change the logic for COLLIDE , we will detect a
# collision if the walker goes outside the RADIUS
#000001111100000
#000111111110000
#001111111111000
#011111111111100
#011111111111100
#001111111111000
#000111111110000
#000001111100000








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

        factor = factor + 1  ### next time multiply by 2,3,4,5...?


main()

    # store the OVERAL time
    # the T_final
    # F final


#semi_main() ;print("semi main()")



#print("main()")
