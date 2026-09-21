## constants
from ast import walk
from math import *

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageSequence

import test_extract_connectivity_info as yiteng

# it is much smaller, depending on the resolution of micro-ct
RESOLUTION = 0.1  # let the lenght of the pixcel be 2 micro-meter
  # s = 0.2 X L .... L is lenght of one pixcel in micro meteres
fluid_diffusion_coefficient = 2.5e3 # micro meter^2 / second
surface_relaxivity = 20  # micro-meter/second
RADIUS = 50 # pixcels given from doctor
RADIUS_IN_MICRO_METERS = RADIUS * RESOLUTION

pore_value_to_color ={
        '1' : 'black' ,
        '2' : 'red' ,
        '3' : 'green' ,
        '4' : 'blue' ,
        '5' : 'purple' ,
        '6' : 'yellow' ,
        '7' : 'orange'
}



def calculate_increment_time(step_distance ,fluid_diffusion_coefficient ):

    delta_t = (step_distance * step_distance) / (6 * fluid_diffusion_coefficient)
    #print(f"\ncalcuting Delta T --> {delta_t}")
    return delta_t
def plot_all_walkers(walker_data, iteration , time , fraction , acceptance_rate , number_of_immagrants ):
    '''
    Input: walker data [(r,c,z,is_alive)......]
    Take Only r , c and use scatter()
    '''
    array_of_pore_values = (np.unique(walker_data[: , 4])).astype(np.int16)

    # print(f"x\n{x[:4]}\ny\n{y[:4]}")
    plt.figure(figsize=(12,10))
    plt.title(f"{iteration} Iterations\n{time:.5f} current Time\n{fraction:.5f} Magnatization\nAcceptance Rate: {acceptance_rate*100}%\n{number_of_immagrants} immagrents")

    is_alive = (walker_data[: , 3]).astype(np.int16)#;print(f"is_alive is\n{is_alive}")
    mask = is_alive == 1
    x = walker_data[mask][: , 0]
    y = walker_data[mask][: , 1]
    for value in array_of_pore_values: # assuming value is string
        color = pore_value_to_color[str(value)]
        mask_for_this_pore_value = walker_data[mask][: , 4] == value

        plt.scatter(x[mask_for_this_pore_value],y[mask_for_this_pore_value], color = color, s = 2 , alpha = 0.42 )



    plt.xlabel("X-axis");plt.ylabel("Y-axis")
    if iteration > 8000 and np.random.default_rng().random() > 0.6  :
        os.system(f"mkdir walker_scatter_plot_figures ;cd walker_scatter_plot_figures")
        plt.savefig(f"walker_scatter_plot_figures/{iteration}_iterations_{fraction}_fraction.png")
    plt.show(block = False)
    plt.pause(1.8)
    plt.close()

def surface_relaxation_for_simulation(p_t):


    x = p_t[: , 1]    # x is TIME ... p_t (ratio , time)
    y = p_t[:,0]
    #y = np.log10(y)
    plt.plot(x,y ,'purple'  , label = "simulation")


def surface_relaxation_for_simulation_log(p_t):


    x = p_t[: , 1]    # x is TIME ... p_t (ratio , time)
    y = p_t[:,0]
    y = np.log10(y)
    plt.plot(x,y ,'purple'  , label = "simulation")


def surface_relaxation_for_analytical(coordinates):
    (x , y) = coordinates[: , 0] , coordinates[: , 1]
    plt.plot(x,y , 'o' , label = "analytical")
def analytical_approach(surface_relaxivity , final_time , radius_in_micro_meter, iterations ):
    # we are going to model the equation  M/M = exp(-3pt/r) , were t=0 --> t=final_t || final_t is aquired AFTER the simulation (when did we stop)
    # this method is used as a benshmark towards the simulation ...
    t = np.linspace(0 , final_time ,iterations )
    y =  np.exp((-3 * t * surface_relaxivity) / radius_in_micro_meter )
    #y = np.log10(y) # base _ 10
    return np.column_stack((t.T , y.T))

def analytical_approach_log(surface_relaxivity , final_time , radius_in_micro_meter, iterations ):
    # we are going to model the equation  M/M = exp(-3pt/r) , were t=0 --> t=final_t || final_t is aquired AFTER the simulation (when did we stop)
    # this method is used as a benshmark towards the simulation ...
    t = np.linspace(0 , final_time ,iterations )
    y =  np.exp((-3 * t * surface_relaxivity) / radius_in_micro_meter )
    y = np.log10(y) # base _ 10
    return np.column_stack((t.T , y.T))


def surface_relaxation(p_t , coordinates):
    plt.figure()
    surface_relaxation_for_simulation(p_t)
    surface_relaxation_for_analytical(coordinates)
    plt.legend(loc = "upper right")


def get_array_from_3D_image(path="/home/saad/Desktop/single_pore.tif"):
    img = Image.open(path)
    full_array = np.array([np.array(page) for page in ImageSequence.Iterator(img) ] )
    print(f"\n(1)3D array has axis {full_array.ndim} with shape {full_array.shape} and has {full_array.size} Voxcels\n")
    return full_array

def get_indexes_of_pore_from_3D_array(array , representation_value_of_pore):
    indexes = [] # unknown number   THUS THE INDEX ARE IN FORMAT
    '''--> [ (r , c , z , Value) ... ] '''
    lz = len(array) ; lr = len(array[0]) ; lc = len(array[0][0])
    # print(f"lz|lr|lc-->{lz}|{lr}|{lc}"); time.sleep(2)
    # x=0
    for z in range(lz):
        for r in range(lr):
            for c in range(lc):
                if array[z,r,c] == np.int64(representation_value_of_pore)  :    # pores are 1     grains are 0
                    indexes.append( (r , c , z , representation_value_of_pore ) )
                    # print(f"index z r c is {indexes[z,r,c]}")
                    # x+=1
    # print(f"x is {x}")
    # time.sleep(2.9)


    return np.array(indexes)
import time
def get_value_of_index(full_array , index):
    r , c , z = index.T[0] , index.T[1] , index.T[2]
    #print(f"r|c|z have shape of {r.shape}|{c.shape}|{z.shape}")
    #print(f"full_array[z,r,c] -->\n{np.unique(full_array[z , r , c])}")
    #time.sleep(5)
    return full_array[z , r , c]

def get_number_of_walkers_migrated(full_array , new_position , old_position , acceptance_rate):
    #print(f"newp->\n{new_position[:3]}\noldp->\n{old_position[:3]}")
    new_coordinate = get_value_of_index(full_array , new_position) # has less valid walkers , as some might got in grain !
    #print(f"new-->\n{new_coordinate}")
    old_coordinate = get_value_of_index(full_array , old_position)

    # print(f"shape of new_index | old_index --> {new_position.shape} | {old_position.shape} ")
    # shape of new_index | old_index --> (3500, 3) | (3500, 3)
    mask_for_walkers_not_in_grain = new_coordinate != 0


    #print(f"old-->\n{old_coordinate}\n\n new : {len(new_coordinate)}\nold:{len(old_coordinate)}");time.sleep(2)
    number_of_illegal_immigration = np.sum((new_coordinate[mask_for_walkers_not_in_grain] != old_coordinate[mask_for_walkers_not_in_grain] ).astype(np.int16))
    # print(f"Mask for Imigrants --> { (new_coordinate[mask_for_walkers_not_in_grain] != old_coordinate[mask_for_walkers_not_in_grain]).shape }")
    #  Mask for Imigrants --> (3168,)

    # print(f"new{np.unique(new_coordinate[:])} | old{np.unique(old_coordinate[:])}")

    ''' we found the bug ! '''
    #new[0 1] | old[1]

    mask_for_walkers_that_migrated = new_coordinate[mask_for_walkers_not_in_grain] != old_coordinate[mask_for_walkers_not_in_grain]
    # print(f"migration mask\n{mask_for_walkers_that_migrated[:5] }") ; time.sleep(1)
    # migration array
    # [False False False False]
    #   w1     w2    w3   w4
    # Thus if walker has TRUE , then he went to another PORE TYPE
    #if number_of_illegal_immigration != len(mask_for_walkers_that_migrated == True):
    #    print(f"N != len(array)")
    '''
    if NOT allowed to migrate , we need to return to their prev position
    new_position = [ (x,y,z) (x,y,z) (x,y,z) ]
    mask_migrated= [  false    true    true  ]
    disallow_migr= [    X      true     false] --> F: allow migration
    '''
    random_array = np.random.default_rng().random(size = number_of_illegal_immigration )

    disallow_migration_mask = random_array > acceptance_rate # this will result in True False
    prevented = np.sum(disallow_migration_mask.astype(np.int16))
    # print(f"\n\nOld Position is\n{old_position[mask_for_walkers_not_in_grain][mask_for_walkers_that_migrated]}\nnew_position migrated are\n{new_position[mask_for_walkers_not_in_grain][mask_for_walkers_that_migrated]}")
    #   Old Position is
    #   [[59 51 47]]
    #   new_position migrated are
    #   [[59 52 46]]
    prevented_walkers_mask = np.zeros(len(new_position), dtype=bool)
    if prevented >= 1:
        valid_indices = np.flatnonzero(mask_for_walkers_not_in_grain)

        migrated_indices = valid_indices[mask_for_walkers_that_migrated]

        prevented_indices = migrated_indices[disallow_migration_mask]

        new_position[prevented_indices] = old_position[prevented_indices]

        prevented_walkers_mask[prevented_indices] = True

        # print(f"we Prevented {prevented} from migrating !\nnewPosition\n{new_position[mask_for_walkers_not_in_grain][mask_for_walkers_that_migrated]}")
        # print(f"disallowed{disallow_migration_mask.shape}")
        # print(f"the updated position in new walker is-->\n{new_position[mask_for_walkers_not_in_grain][mask_for_walkers_that_migrated][disallow_migration_mask]}")
    # time.sleep(0.8)
    return number_of_illegal_immigration - prevented , new_position , prevented_walkers_mask


def convert_index_to_mid_point_position(index , resolution):
    #print(f"index is {len(index)}")
    (r , c , z) = index

    # instead of 12.3000000000000000004 --> round() ---> 12.300
    pr = round((r) * resolution + (0.5 * resolution) , 3)
    pc = round((c) * resolution + (0.5 * resolution) , 3)
    pz = round((z) * resolution + (0.5 * resolution) , 3)
    return (pr , pc , pz)


def get_initiated_walkers(full_array,representation_value_of_pore,resolution,number_of_walkers = 1000):
    # print(f"we are calling get index of pore ....")
    # print(f"value of pore --> {representation_value_of_pore} checking -problem-")
    pores_indexes = get_indexes_of_pore_from_3D_array(full_array, representation_value_of_pore)
    number_of_pores = len(pores_indexes)
    #print(f"the max index in pores is\n {max(pores_indexes.T[0])}");time.sleep(1)
    #           (520268 , 4)
    # print(f"indexes of all pores has shape of {pores_indexes.shape} ")

    # time.sleep(4)
    # number_of_walkers = 1000
    walker_data = np.zeros( ( number_of_walkers , 5 ) , dtype = np.float64 )

    '''
    NOW WE DETERMINE NUMBER OF WALKER , BUT IF WE WANT OTHERWISE , THEN WE MOD THE CODE !!
    '''

    # THIS IS working as we want it to #
    for i in range(number_of_walkers):
        # this masteriously converted the indexes INSIDE TO STR ?!!?
        # print(f"index {pores_indexes[i][:3]} ")
        index = np.int16(pores_indexes[i][:3])
        # print(f"index to be converted is\n{type(index[0])}");time.sleep(2) # we will get the first pore's pixcel , then second ..
        # 2 we calculate the MIDDLE POSITION OF THIS PORE (4 , 5 , 8) --> (4.5 , 5.5 , 8.5)
        position_for_uninitialized_walker = convert_index_to_mid_point_position(index , resolution) + (1,)  # so now it is like this (x , y , z , is_alive=1)
        # print(f"we converted {index} to {position_for_uninitialized_walker} for walker[{i}]") ; time.sleep(0.4)
        walker_data[i] = position_for_uninitialized_walker + (representation_value_of_pore,)

    return walker_data



# this method did not work , cause it produced a stairs in the simulation !!


def return_collided_walkers(converted_to_index_new_walker_data , full_array , representation_value_of_grain = 0 ) :
    converted_to_index_new_walker_data= np.int32(converted_to_index_new_walker_data)
    #print("\n Representing Walkers index's ... we will take only 5 walkers and show them\n<===============>\n")
    #print(f"\nall walkers here are converted to Index of their voxcels--->{converted_to_index_new_walker_data[0:5]}..\n...")
    all_z_index_levels_walkers_at = converted_to_index_new_walker_data[: , 2]
    all_c_index_levels_walkers_at = converted_to_index_new_walker_data[: , 1]
    all_r_index_levels_walkers_at = converted_to_index_new_walker_data[: , 0]
    #print(f"All Z --> {all_z_index_levels_walkers_at[:5]}..\nAll R --> {all_r_index_levels_walkers_at[:5]}..\nAll C --> {all_c_index_levels_walkers_at[:5]}..\n")
    #print(f"type{all_z_index_levels_walkers_at[0]}")
    #time.sleep(10)
    result = full_array[ all_z_index_levels_walkers_at , all_r_index_levels_walkers_at , all_c_index_levels_walkers_at  ]
    # now we have the VALUES of the pixcels .. there might be 1 , 2, 3 ..etc as values
    # # we do not know which is grain , unless TOLD
    #print(f"\nThese voxcels have grains -->\n")
    # first detech the GRAIN , then change their value --> unique = MAX() + 1
    represent_collision_with_this_number = np.max(result)+1
    # we replace , each GRAIN value , with this UNIQUE number
    result[result == representation_value_of_grain]  = represent_collision_with_this_number
    result = result // represent_collision_with_this_number  # all other values are 0 , cause they are smaller , except GRAINS --> 1 which is COLLISION

    '''
    xxx=len(converted_to_index_new_walker_data)
    for i in range(xxx) :
       #print(f"{converted_to_index_new_walker_data[i]}")
       #time.sleep(14)
       (r , c , z , a ) = converted_to_index_new_walker_data[i]
       #print("index of pixcels is :" , r , c , z , a)
       if r > 120 or c > 120 or z > 120 :
               continue
           if full_array[int(z)][int(r)][int(c)] == representation_value_of_grain :
                   if result[i] != 1 :
                           print(f"incorrect calculation of is_collide")
                           raise(Exception)
                           '''
    return result


# the walkers that have is_alive = 1 , are kept .. the others are gone !

def convert_position_to_index(walker_data , resolution):
    #the format of walker_data is [[nr , nc , nz ] () () ..()]


    indexed_walker_data = np.copy(walker_data)
    # operate on all
    indexed_walker_data = np.floor(indexed_walker_data / resolution).astype(np.int32)
    # print(f"indexed walker data===>\n{indexed_walker_data[:3]}");time.sleep(3)


    return indexed_walker_data
import math
#Smooth sphere center: 12.1 11.9 11.9
#Smooth sphere radius: 10.0
#Voxel-center walkers outside smooth sphere: 35
def get_index_of_each_pore_type(full_array,total_number_of_different_pores,values_of_each_pore):
    index_of_each_pore = {} # pore_1 : np.array([----])
    for i in range(total_number_of_different_pores):
        # print(f"value of pore --> {values_of_each_pore[i]} checking..")
        # print(f"number of different types of pores --> {total_number_of_different_pores}\nValue of pore --> {values_of_each_pore[i]}\n")
        index_of_each_pore[str(i+1)] = get_indexes_of_pore_from_3D_array(full_array, values_of_each_pore[i])
    return index_of_each_pore
def calculate_the_volume_of_the_given_sphere_image(full_array,representation_value_of_pore  , resolution):
    number_of_pore_voxels = len(get_indexes_of_pore_from_3D_array(full_array , representation_value_of_pore))
    volumn_of_one_pore_voxels = resolution * resolution * resolution
    volumn_of_all_pore_voxels = volumn_of_one_pore_voxels * number_of_pore_voxels
    return volumn_of_all_pore_voxels


def calculate_the_surface_area_of_the_sphere(full_array , resolution ):
    all_z_dimentions_last_n_minus_1 = full_array[1: , : , :]
    all_z_dimentions_first_n_minus_1 = full_array[:-1 , : , :]

    all_x_dimentions_last_n_minus_1  = full_array[: , 1: , :]
    all_x_dimentions_first_n_minus_1 = full_array[: , :-1, :]

    all_y_dimentions_last_n_minus_1  = full_array[: , : , 1:]
    all_y_dimentions_first_n_minus_1 = full_array[: , : , :-1]

    number_of_surfaces_facing_grains = 0
    number_of_surfaces_facing_grains += np.sum((all_z_dimentions_last_n_minus_1 != all_z_dimentions_first_n_minus_1 ).astype(np.int16))
    number_of_surfaces_facing_grains += np.sum((all_x_dimentions_last_n_minus_1 != all_x_dimentions_first_n_minus_1 ).astype(np.int16))
    number_of_surfaces_facing_grains += np.sum((all_y_dimentions_last_n_minus_1 != all_y_dimentions_first_n_minus_1 ).astype(np.int16))

    surface_area_of_one_face = resolution * resolution
    surface_area_of_all_pores_facing_grains = surface_area_of_one_face * number_of_surfaces_facing_grains
    print(f"\n===========\nsurface area is {surface_area_of_all_pores_facing_grains}")
    #time.sleep(0.2)
    return surface_area_of_all_pores_facing_grains
import re



# /home/saad/Desktop/programs/slb/tmp/images/ThreePoreSystem_20_50_20_0p2um_shift_5.tif

def semi_semi_main(   iterations = 600):
    # get the image , and store it in np_array
    file_name="ThreePoreSystem_20_50_20_0p2um_shift_5.tif"
    file_name="MultiPoreSystem_20_20_20_50_20_20_20_0p2um_shift_5.tif"
    file_name="TwoPoreSystem_50_20_0p2um_shift_3.tif"
    values_in_file_name = re.findall(r"_(\d+)(?:p(\d+))?" , file_name)
    #    [('50', ''), ('50', ''), ('0', '2'), ('5', '')]
    number_of_voxels_in_the_radius  = int(values_in_file_name[0][0])
    resolution = int(file_name[file_name.find("0p") + 2]) / 10
    radius_in_micro_meters          = number_of_voxels_in_the_radius * resolution
    ''' ^|^   this is used for Analytical Calculation'''
    full_array                      = yiteng.get_full_array(file_name=file_name)
    total_number_of_different_pores = yiteng.get_total_number_of_different_pores(folder_path = "images", file_name=file_name)
    values_of_each_pore             = yiteng.get_the_values_of_each_pore(folder_path = "images" , file_name=file_name)
    step_distance = resolution  *  1
    print(f"we have {total_number_of_different_pores} number of different pores.\nThey are::->\n{values_of_each_pore}")
    print(f"(2)Resolution -> {resolution}");print(f"(3)step_distance -> {step_distance}");print(f"(4)fluid_diffusion_coefficient -> {fluid_diffusion_coefficient} micro_meter^2/sec");print(f"(5)surface_relaxivity -> {surface_relaxivity} micro_meter/sec");print(f"(6)radius -> {number_of_voxels_in_the_radius} voxels")

    '''
    changeable hyper parameter
    number of walkers in each pores
    the pores that we are gonna use
    '''

    # {   '1' : np.array()   '2' : np.array()     '3' : np.array()   }
    indexes_of_pores_dataset = get_index_of_each_pore_type(full_array,total_number_of_different_pores,values_of_each_pore)

    all_walkers_in_one_type_pore = False

    if all_walkers_in_one_type_pore:
        pore_type_to_walker_amount = {
            '1' : 3500 , # pore_1 -> 3500 walkers
            '2' : 0 ,
            '3' : 0
        }
    else :
        di = len(indexes_of_pores_dataset)
        total_number_of_walker = 30000
        pore_type_to_walker_amount = {
        '1' : 1500 , # pore_1 -> 3500 walkers
        '2' : 2000
        }
        for pore_value in indexes_of_pores_dataset :
            pore_type_to_walker_amount[pore_value] = total_number_of_walker // di
    print(f"pores<-->walkers\n{pore_type_to_walker_amount}"); time.sleep(1.4)
    number_of_walkers = sum(pore_type_to_walker_amount.values())
    walker_data = get_initiated_walkers(full_array ,'1' , resolution , pore_type_to_walker_amount['1'] )
    # if we have OTHER PORE TYPES , we ADD more walkers !
    del pore_type_to_walker_amount['1']
    # print(f"walkers data is \n{walker_data[:4]}");time.sleep(4.2)
    if sum(pore_type_to_walker_amount.values()) >= 1 :
        print("Other Pores Will have Walkers in them..");time.sleep(1)
        for pore_type , walker_amount in pore_type_to_walker_amount.items() :
            temp = get_initiated_walkers(full_array , pore_type , resolution , walker_amount)
            # print(f"walkers data is \n{walker_data[-4:-1]}\nporeType {pore_type}\nAmount of walkers {walker_amount}");time.sleep(4.2)
            walker_data = np.vstack((walker_data , temp))
        # print(f"walkers data is \n{walker_data[:4]}");time.sleep(1.2)
    else:
        print("Other Pores Will |- NOT -| have Walkers in them..");time.sleep(1)
    # print(f"shape of walker_data after vstacking is\n{walker_data.shape}")
    # time.sleep(4.2)

    representation_value_of_grain = 0

    seed = 99

    rng = np.random.default_rng(seed)



    current_live_walkers = number_of_walkers

    initial_population_walkers = number_of_walkers
    t = 0

    delta_t = calculate_increment_time(step_distance ,fluid_diffusion_coefficient)
    # should be 900 iterations



    p_fraction = np.zeros( (iterations , 2 ) , dtype = np.float64 )

    length_of_walker_data_array = len(walker_data)
    likely = np.ones(length_of_walker_data_array) * (2 * step_distance * surface_relaxivity / (3 * fluid_diffusion_coefficient))
    start = time.time()
    walkers_migrated = 0
    for jj in range(iterations):

        aaa = time.time()

        r = walker_data[:, 0]
        c = walker_data[:, 1]
        z = walker_data[:, 2]
        is_alive = walker_data[:, 3]

        # walker_data[: , :3 ]
        # theta  = np.ones(len(walker_data))  # [ th , th , th ... th ]
        # beta   = np.ones(len(walker_data))  # [ be , be , be ... be
        theta = rng.random(length_of_walker_data_array) * math.pi * 2
        beta = rng.random(length_of_walker_data_array) * math.pi
        # WE EXTRACT X Y Z SEPRATELY

        nr = r + (step_distance * np.sin(beta) * np.cos(theta))
        nc = c + (step_distance * np.sin(beta) * np.sin(theta))
        nz = z + (step_distance * np.cos(beta))


        new_walker_data = np.column_stack((nr, nc, nz))


        converted_to_index_new_walker_data = convert_position_to_index( new_walker_data , resolution  )  # from POSITION --> Index
        converted_to_index_old_walker_data = convert_position_to_index( walker_data[: , :3 ] , resolution )
        acceptance_rate = 0.5
        # HERE WE prevent migration .. thus Update on NewWalkerData ..
        walkers_migrated, converted_to_index_new_walker_data  , prevented_walkers_mask  = get_number_of_walkers_migrated(full_array , converted_to_index_new_walker_data , converted_to_index_old_walker_data  , acceptance_rate)
        ''' do this func() again ! '''
        # walkers who are not allowed to pass , shall return to their origianl positions
        nr[prevented_walkers_mask] = r[prevented_walkers_mask]
        nc[prevented_walkers_mask] = c[prevented_walkers_mask]
        nz[prevented_walkers_mask] = z[prevented_walkers_mask]
        converted_to_index_new_walker_data = np.column_stack((converted_to_index_new_walker_data , is_alive))
        # print(f"converted to index new walker \n{converted_to_index_new_walker_data[:3]}")
        # time.sleep(5)
        collided = return_collided_walkers(converted_to_index_new_walker_data, full_array,representation_value_of_grain)
        # print(f"collided is {collided[:3]}");time.sleep(9)

        # IF IT IN TOUCHING OR BEHOND THE GRAIN, WE CALCULATE TEH LIKELY
        #likely = np.ones(lenght_of_collided_array) * (2 * step_distance * surface_relaxivity / (3 * fluid_diffusion_coefficient))

        # from 0.17 ---> 0.25
        # boost_number = 0.18 # we decrease the amplitude of the numbers in random_numbers
        random_number = rng.random(size=length_of_walker_data_array)
        is_dead =  (random_number  < likely  )
        ## IF THERE IS NO COLLISION , THEN r <- nr
        r[collided == 0] = nr[collided == 0]
        c[collided == 0] = nc[collided == 0]
        z[collided == 0] = nz[collided == 0]
        ## IF THERE IS COLLISION AND DEAD , THEN is_alive=0

        number_of_walker_before_killing = np.sum(is_alive)
        is_alive[(is_dead == 1) & (collided == 1)] = 0  # here we KILL
        number_of_walkers_after_killing = np.sum(is_alive)

        number_of_dead_in_this_bach = number_of_walker_before_killing - number_of_walkers_after_killing
        current_live_walkers = current_live_walkers - number_of_dead_in_this_bach
        fraction = (current_live_walkers / initial_population_walkers )  # p(t) =  N_1 / N_o
        true_fraction = math.exp((-3 * t * surface_relaxivity) / radius_in_micro_meters )
        t = t + delta_t



        if current_live_walkers == 0:
            break

        p_fraction[jj] = (fraction, t)

        bbb = time.time()

        # this insures that we OUTPUT 4500 ROWS for PRINTING REgardless of # of iterations
        c = 1 + (iterations // 4500)
        if jj % ( c if c > 2 else 2 )  == 0 :
            print(f"{jj}|| T:{t} || F:{fraction} || iteration time: {bbb - aaa} || diff : {true_fraction - fraction}")
            if np.random.default_rng().random() > 0.95 :
                # print(f"walker_data AFTER Update-->\n{walker_data[:4]}");time.sleep(2)
                plot_all_walkers(walker_data , jj , t , fraction , acceptance_rate , number_of_immagrants=walkers_migrated )

            # is_alive = 1 :: for being alive      is_alive = 0 for being dead




    end = time.time()

    print(f"\nOVERALL TIME {end - start} ") ; print(f"\n=================\nWalkers Migrated\n=================\n{walkers_migrated}\n")
    time.sleep(1.2)
    coordinates = analytical_approach(surface_relaxivity , t , radius_in_micro_meters , iterations )
    time.sleep(1.5)

    return (p_fraction, coordinates , end - start , walkers_migrated )



def print_all_coordinations(p_t , coordinates):
    x_s = p_t.T[1]
    y_s = p_t.T[0]
    #print(f"\nCoordinates is ::\n{coordinates}")
    # (x_a , y_a) = coordinates[: , 0] , coordinates[: , 1]
    number_of_points = min(len(p_t), len(coordinates))
    s = "\n||=========simulation==========\n"
    c = 1 + (number_of_points // 4500)
    for i in range(number_of_points):
        if i % ( c if c > 2 else 2 )  == 0 :
            s+=f"||t:{x_s[i]:.6f} -> {y_s[i]:.6f} \n"
    print(s)


# still have not tryied this
import logging
import os

'''
if LOG , then change LOTS of things '''
def main():
    final_time = 0.000005
    factor = 1
    while(final_time < 2):

        iterations = 400000 * factor
        #for boost_number in boost_number_list :
        (p_fraction, coordinates , total_time , walkers_migrated ) = semi_semi_main(   iterations)
        final_time = p_fraction[-1][1]
        final_ratio = p_fraction[-1][0]
        print_all_coordinations(p_fraction,coordinates)
        ### here we should retreive the image , name it , store it
        print(f"Iterations:{iterations}\nF_t:{final_time}\nfinal_ration:{final_ratio}")
        time.sleep(0.1)
        try:
                   surface_relaxation(p_fraction, coordinates)
                   os.system(f"cd Figures_Vectorized;mkdir iterations_{iterations}_migration")
                   plt.suptitle(f"{iterations} iterations\n{walkers_migrated} Migrations ")# edgecolor
                   plt.title(f"Final time:{final_time:.5f}sec\nFinal magnetization:{final_ratio:.5f}\nTotal Durations of program:{round(total_time,1)}sec")
                   plt.xlabel("Time Duration")
                   plt.ylabel("scale of magnetization")
                   plt.tight_layout()
                   # plt.savefig(f"Figures_Vectorized/iterations_{iterations}_migration/iterations_is{iterations}_walkers_is3500_res_0p2.png" , dpi = 500 , transparent = False )
                   plt.show(block=False)
                   plt.pause(6.6)
                   plt.close() # to clear out memory

        except Exception as exc:
                   logging.exception("Exception occurred")

        factor = factor * 2


main()
