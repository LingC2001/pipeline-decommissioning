import random as rnd
import numpy as np
import gzip
import os
import time
# Generate a matrix
# Compress the matrix
# Store in a database
# Read and send the matrix

DIRECTROY = 'pipe-recognisation/Architecture/Sensors/'

def generate_matrix (a, b, dim):
    return np.random.random_integers(a, b, dim)

def compress_matrix (arr) :
    np.save(DIRECTROY + 'data.npy', arr)
    print( f" Original size : {os.path.getsize(DIRECTROY + 'data.npy')}" )
    st = time.time()
    np.savez_compressed(DIRECTROY + 'data.npy', a=arr)
    print("Compression time :", time.time() - st, 'seconds')
    print( f"Compressed size : {os.path.getsize(DIRECTROY + 'data.npy.npz')}")
    print( f"Ratio : {os.path.getsize(DIRECTROY + 'data.npy.npz') / os.path.getsize(DIRECTROY + 'data.npy')}")
    


def main () :
    compress_matrix(generate_matrix(0, 100, (1000, 1000)))

main()
