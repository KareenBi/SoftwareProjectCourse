
import argparse
import sys
from numpy.core.fromnumeric import shape
from numpy.core.records import array
import pandas as pd
import numpy as np
from pandas.io.parsers import read_csv
import spkmeans

def main(K, goal, file_name):
    T = spkmeans.fit_spkmeans(file_name, goal, K)
    T = np.array(T)
    if(goal != "spk"):
        return None
    K = T.shape[1] #d=k
    N = T.shape[0]
    centroids = initialize_centroids(T, K, K, N)
    T = T.tolist()
    centorids_list = centroids.tolist()
    centroids =  spkmeans.fit_kmeans(T, centorids_list, N, K, K)
    return None


#initializes the centroids using kmeans++ algorithm
def initialize_centroids(data, K, d, N):
    indices = np.ndarray(K, int)
    centroids = np.ndarray((K,d), float)
    index = np.random.choice(N)
    indices[0] = index
    centroids[0] = data[index]
    Z = 1
    D = np.zeros((N))
    P = np.zeros((N))
    while Z<K:
        for i in range(0, N):
            min = float("inf")
            for j in range(0, Z):
                distance = np.sum(np.power(data[i] - centroids[j], 2))
                if(distance<min):
                    min = distance
            D[i] = min
        sum_D = np.sum(D)
        P = np.divide(D, sum_D)
        index = np.random.choice(N, p=P)
        indices[Z] = index
        centroids[Z] = data[index]
        Z+=1
    print(','.join(str(i) for i in indices), flush=True) 
    return centroids


def start(): #gets arguments and starts the algorethim.
    parser = argparse.ArgumentParser()
    parser.add_argument("K", type=int, help="K is the number of clusters")
    parser.add_argument("goal", type=str, help="Can get the following values: spk, wam, ddg, lnorm, jacobi")
    parser.add_argument("file_name", type=str, help="The path to file which contains N observations, the file extension is .txt/.csv")
    args = parser.parse_args()
    K = args.K 
    goal = args.goal
    file_name = args.file_name
    #assertions
    if K == None:
        print("Invalid Input!")
        return -1
    if(K<0):
        print("Invalid Input!")
        return -1
    if file_name == None:
        print("Invalid Input!")
        return -1
    if (goal == None or goal not in ["spk", "wam", "ddg", "lnorm", "jacobi"]):
        print("Invalid Input!")
        return -1
    main(K, goal, file_name)

np.random.seed(0)
start()