# test_extract_connectivity_info.py
import os

import numpy as np
import porespy as ps
import matplotlib.pyplot as plt
import cv2
import scipy.sparse as sp

# from precompiled_module import get_solid_array, local2global

folder_path = "images"
file_name = "ThreePoreSystem_20_50_20_0p2um_shift_5.tif"
# file_name = "MultiPoreSystem_50_50_50_50_50_50_50_0p2um_shift_5.tif"
file_path = os.path.join(folder_path, file_name)
# ***** initialize image
tf, image_list = cv2.imreadmulti(file_path, flags=0)
bw = np.array(image_list, dtype=np.uint8).transpose((0,2,1))
#print(bw.shape)
# print(np.unique(bw))
# fig1 = plt.figure(figsize=(10,6))
# ax1 = fig1.add_subplot(111)
# # ax1.imshow(bw[60])
# last_slice = [dim //2 * 2 for dim in bw.shape]
# bw = bw[:last_slice[0], :last_slice[1], :last_slice[2]]
slices = [slice(dim-1) if dim % 2 == 1 else slice(dim) for dim in bw.shape]
bw = bw[*slices]
# ***** segment pore space by SNOW algorithm
# snow2 has to crop image so that all side lengths become even
snow = ps.networks.snow2(bw, voxel_size=1, boundary_width=0)
label_mask = snow.regions
dims = label_mask.shape
#print(f"labelmaskunique is {4} shape {dims}\nfirst \n{label_mask[0].shape}\nSecond\n{label_mask[1]}\n\nunique\n{np.unique(label_mask)}")


if __name__ == "__main__" :
    fig2 = plt.figure(figsize=(10,6))
    ax2 = fig2.add_subplot(111)
    ax2.imshow(label_mask[60])
# ax2.imshow(label_mask[60]/bw[60], origin="lower", interpolation="none")
#ax2.imshow(label_mask[160], origin="lower", interpolation="none")
    plt.show()
# ***** extract key parameters of segmented pores and throats
conns = snow.network["throat.conns"]
throat_area = snow.network["throat.cross_sectional_area"]
pore_volume = snow.network["pore.region_volume"]
# calculate coordination number further
m_pore = len(pore_volume)
zcoord = np.bincount(conns.ravel(), minlength=m_pore)
# ***** create sparse matrix
coords = np.vstack((conns, conns[:,::-1]))
ones = np.ones(coords.shape[0], dtype=np.float64)
adj_matrix = sp.coo_array((ones, (coords[:, 0], coords[:, 1])), shape=(m_pore, m_pore), dtype=np.float64)
csr_matrix = adj_matrix.tocsr()
# data = csr_matrix.data
indices = csr_matrix.indices
indptr = csr_matrix.indptr
# an alternative way to calculate the coordination number
# zcoord = np.array([indptr[kk+1]-indptr[kk] for kk in range(m_pore)], dtype=np.int64)
# create a numpy array to store bi-directional connectivity info
data = np.zeros((coords.shape[0], 3), dtype=np.float64)
# data = np.zeros((coords.shape[0], 4), dtype=np.float64)
data[:, 0] = zcoord[coords[:, 1]]                  # coordination number
data[:, 1] = np.hstack((throat_area, throat_area)) # cross-sectional area
data[:, 2] = pore_volume[coords[:, 1]]             # adjacent pore volume
# data[:, 3] = pore_volume[coords].max(axis=1)     # maximal volume between conencted pores
# # ***** locate the solid voxels along solid-pore interface
# !!!!! before we assumed that the pore voxels are zeros and solid voxels are one in the binarized data
# TODO: if we input the labelled data, the pore voxels are nonzero while the solid voxels are zero (this has to be taken into account)
# arr = np.pad(label_mask, 1, "wrap")
# _local = get_solid_array(arr, 1)
# block_size = np.array(_local.shape, dtype=np.int64)[:, np.newaxis]
# solid_index = local2global(_local, 0, 0, 0, block_size, 1, dims)
#print("done")
def base(folder_path,file_name):
    file_path = os.path.join(folder_path, file_name)
    tf, image_list = cv2.imreadmulti(file_path, flags=0)
    bw = np.array(image_list, dtype=np.uint8).transpose((0,2,1))
    slices = [slice(dim-1) if dim % 2 == 1 else slice(dim) for dim in bw.shape]
    bw = bw[*slices]
    snow = ps.networks.snow2(bw, voxel_size=1, boundary_width=0)
    label_mask = snow.regions
    return label_mask
# if there is 5 differen pores , we return 5 !
def get_total_number_of_different_pores(folder_path,file_name):
    label_mask = base(folder_path,file_name)
    total_number_of_all_values_with_grain = len(np.unique(label_mask))
    number_of_values_of_only_pores = total_number_of_all_values_with_grain - 1
    return number_of_values_of_only_pores



def get_the_values_of_each_pore():
    label_mask = base(folder_path,file_name)
    return np.unique(label_mask)[1:]


def get_full_array(folder_path = "images" , file_name="TwoPoreSystem_50_50_0p2um_shift_5.tif"):
    label_mask = base(folder_path,file_name)
    return label_mask
