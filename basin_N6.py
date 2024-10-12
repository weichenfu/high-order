import numpy as np
import time
import math
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
from PIL import Image
from scipy.integrate import solve_ivp
from matplotlib import rcParams
from scipy import integrate
import random
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from concurrent.futures import ProcessPoolExecutor
import os
import copy

### generate_n_dimensional_grid data
def generate_n_dimensional_grid(N, num):
    linspace = np.linspace(0.0001, np.pi-0.0001, num)
    grid = np.meshgrid(*[linspace] * N, indexing='ij')
    grid_reshaped = np.stack(grid, axis=-1).reshape(-1, N)
    return grid_reshaped

def main(N):
    n_dimensional_grid =generate_n_dimensional_grid(N, num)
    n_dimensional_grid[n_dimensional_grid<=0.01]=0.1
    n_dimensional_grid[n_dimensional_grid >= math.pi-0.01] = math.pi-0.1
    np.savetxt('6_dimensional_grid.txt', n_dimensional_grid, delimiter=',', header=str(n_dimensional_grid.shape))
if __name__ == "__main__":
    N = 6
    num=10
    main(N)

N=6
k = 16
k_s = 7.0
delta = 0.001

def diff_equation(y, k, k_s):
    sin_y = np.sin(y[:, np.newaxis] - y)
    dy = k / N * np.sum(sin_y, axis=0) - k_s * np.sin(2 * y)
    return dy

def euler_method(y0, k, k_s, xmax=200, node=40000):
    t = np.linspace(0, xmax, num=node)
    h = t[1] - t[0]  # 计算步长
    y = np.array(y0)
    results = np.zeros((node, len(y0)))
    results[0] = y0
    # result1=y0
    for i in range(1, node):
        y += h * diff_equation(y, k, k_s)
        results[i] = y
        if np.sum(np.sin(np.round(results[i], 5)))<0.0001:
            results[-1] = results[i]
            # print(i)
            break
    return results[-1]

def main(N):
    n_dimensional_grid = np.loadtxt('6_dimensional_grid.txt', delimiter=',')
    count1 = count2 = count3 = count4 = count5 = count6 = count7 = 0
    a = len(n_dimensional_grid)
    # print(a)
    for i in range(a):
        results = euler_method(n_dimensional_grid[i, :], k, k_s)
        print(i)
        # print(np.round(results, 5))
        if np.abs(np.mean(results)) <= delta and np.abs(np.std(results)) <= delta:
            count1 = count1 + 1
        if np.abs(np.mean(results)-np.pi) <= delta and np.abs(np.std(results)) <= delta:
            count2= count2+1
        contains_zero = np.count_nonzero(results <= delta)
        contains_pi = np.count_nonzero(results >= np.pi-delta)
        if np.all(contains_zero | contains_pi) and contains_zero!=N and contains_pi!=N and contains_zero > contains_pi:  #### 0 more than pi
            count3 = count3 + 1
            # print(np.round(results, 5))
            # print(n_dimensional_grid[i, :])
        if np.all(contains_zero | contains_pi) and contains_zero!=N and contains_pi!=N and contains_zero < contains_pi:   #### pi more than 0
            count4 = count4 + 1
            # print(np.round(results, 5))
        if np.all(contains_zero | contains_pi) and contains_zero != N and contains_pi != N and contains_zero == contains_pi:  #### pi equal 0
            count5 = count5 + 1
            # print(np.round(results, 5))
        if np.all(np.abs(np.mean(k / N * np.sum(np.sin(results[:, np.newaxis] - results), axis=0) - k_s * np.sin(2 * results)))<0.0000001 and contains_zero==0 and contains_pi==0): #### other equilibrium
            ####(theta,pi-theta)
            t = copy.deepcopy(list(results))
            # three max value
            max_number = []
            for _ in range(3):
                number = max(t)
                max_number.append(number)
            # print(max_number)
            # three min value
            min_number = []
            for _ in range(3):
                number = min(t)
                min_number.append(number)
            # print(min_number)
            if max_number[0] == max_number[1]==max_number[2] and min_number[0] == min_number[1]==min_number[2] and abs(
                    np.sum(np.cos(results), axis=0)) <= delta:
                count6 = count6 + 1
                print(results)
            if abs(np.sum(np.cos(results), axis=0)) > delta:
                count7 = count7 + 1
                print("other equilibrium", np.round(results, 5))

    print("N=:", N, "k=:", k, "k_s=:", k_s)
    print("Number of initial conditions converging to (0,0):", count1)
    print("Number of initial conditions converging to (pi,pi):", count2)
    print("Number of initial conditions converging to (0,pi) (0 more than pi):", count3)
    print("Number of initial conditions converging to (0,pi) (pi more than 0):", count4)
    print("Number of initial conditions converging to (0,pi) (pi equal 0):", count5)
    print("Number of initial conditions converging to (theta,pi-theta):", count6)
    print("Number of initial conditions converging to others:", count7)

if __name__ == "__main__":
    main(N)


