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


### generate_n_dimensional_grid
def generate_n_dimensional_grid(N, num):
    linspace = np.linspace(0, math.pi, num)
    grid = np.meshgrid(*[linspace] * N, indexing='ij')
    grid_reshaped = np.stack(grid, axis=-1).reshape(-1, N)
    return grid_reshaped

def main(N):
    n_dimensional_grid =generate_n_dimensional_grid(N, num)
    n_dimensional_grid[n_dimensional_grid<=0.01]=0.1
    n_dimensional_grid[n_dimensional_grid >= math.pi-0.01] = math.pi-0.1
    np.savetxt('3_dimensional_grid.txt', n_dimensional_grid, delimiter=',', header=str(n_dimensional_grid.shape))
if __name__ == "__main__":
    N = 3
    num=20
    main(N)

N = 3
k = 8
k_s = 4.1
delta = 0.001

def diff_equation(y, k, k_s):
    sin_y = np.sin(y[:, np.newaxis] - y)
    dy = k / N * np.sum(sin_y, axis=0) - k_s * np.sin(2 * y)
    return dy

def euler_method(y0, k, k_s, xmax=100, node=10000):
    t = np.linspace(0, xmax, num=node)
    h = t[1] - t[0]
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
    n_dimensional_grid = np.loadtxt('3_dimensional_grid.txt', delimiter=',')
    count1 = count2 = count3 = count4 = count5 = 0
    result_matrix1=np.empty((0,N))
    result_matrix2=np.empty((0,N))
    result_matrix3=np.empty((0,N))
    result_matrix4= np.empty((0,N))
    result_matrix5=np.empty((0,N))
    a = len(n_dimensional_grid)
    # print(a)
    for i in range(a):
        results = euler_method(n_dimensional_grid[i, :], k, k_s)
        print(i)
        # print(np.round(results, 5))
        if np.abs(np.mean(results)) <= delta and np.abs(np.std(results)) <= delta:
            count1 = count1 + 1
            result_matrix1=np.vstack((result_matrix1, n_dimensional_grid[i, :]))
        if np.abs(np.mean(results)-np.pi) <= delta and np.abs(np.std(results)) <= delta:
            count2= count2+1
            result_matrix2=np.vstack((result_matrix2, n_dimensional_grid[i, :]))
        contains_zero = np.count_nonzero(results <= delta)
        contains_pi = np.count_nonzero(results >= np.pi-delta)
        if np.all(contains_zero | contains_pi) and contains_zero!=N and contains_pi!=N and contains_zero > contains_pi:  #### 0 more than pi
            count3 = count3 + 1
            result_matrix3 = np.vstack((result_matrix3, n_dimensional_grid[i, :]))
        if np.all(contains_zero | contains_pi) and contains_zero!=N and contains_pi!=N and contains_zero < contains_pi:  #### pi more than 0
            count4 = count4 + 1
            result_matrix4 = np.vstack((result_matrix4, n_dimensional_grid[i, :]))
        if np.all(np.abs(np.mean(k / N * np.sum(np.sin(results[:, np.newaxis] - results), axis=0) - k_s * np.sin(2 * results)))<0.0000001 and contains_zero==0 and contains_pi==0): #### other equilibrium
            count5 = count5 + 1
            result_matrix5 = np.vstack((result_matrix5, n_dimensional_grid[i, :]))
            print(results)

    print("N=:", N, "k=:", k, "k_s=:", k_s)
    print("Number of initial conditions converging to (0,0):", count1)
    print("Number of initial conditions converging to (pi,pi):", count2)
    print("Number of initial conditions converging to (0,pi) (0 more than pi):", count3)
    print("Number of initial conditions converging to (0,pi) (pi more than 0):", count4)
    print("Number of initial conditions converging to others:", count5)

    print(np.shape(result_matrix1))
    print(np.shape(result_matrix2))
    print(np.shape(result_matrix3))
    print(np.shape(result_matrix4))

    np.savetxt('zero_k8_ks4point1.txt', result_matrix1, delimiter=',', header=str(result_matrix1.shape))
    np.savetxt('pi_k8_ks4point1.txt', result_matrix2, delimiter=',', header=str(result_matrix2.shape))
    np.savetxt('first0pi_k8_ks4point1.txt', result_matrix3, delimiter=',', header=str(result_matrix3.shape))
    np.savetxt('second0pi_k8_ks4point1.txt', result_matrix4, delimiter=',', header=str(result_matrix4.shape))


if __name__ == "__main__":
    main(N)


#### figure
zero_rows = np.loadtxt('zero_k8_ks4point1.txt', delimiter=',')
pi_rows = np.loadtxt('pi_k8_ks4point1.txt', delimiter=',')
first_rows = np.loadtxt('first0pi_k8_ks4point1.txt', delimiter=',')
second_rows = np.loadtxt('second0pi_k8_ks4point1.txt',delimiter=',')

fig = plt.figure(figsize=(7, 7))
ax = fig.add_subplot(111, projection='3d')

# Plot the points with different colors
ax.scatter(zero_rows[:, 0], zero_rows[:, 1], zero_rows[:, 2], c='r')   #, label=r'$\varphi^{ps}$')
ax.scatter(pi_rows[:, 0], pi_rows[:, 1], pi_rows[:, 2], c='g')  #, label=r'$\varphi^{ps}_{*}$')
ax.scatter(first_rows[:, 0], first_rows[:, 1], first_rows[:, 2], c='b')  #, label=r'$\varphi^{bp}_{\text{I}}$')
ax.scatter(second_rows[:, 0], second_rows[:, 1], second_rows[:, 2], c='y')  #, label=r'$\varphi^{bp}_{\text{II}}$')

ax.set_xlabel(r'$\varphi_1$',fontsize=18,fontname="Times New Roman")
ax.set_ylabel(r'$\varphi_2$',fontsize=18,fontname="Times New Roman")
ax.set_zlabel(r'$\varphi_3$',fontsize=18,fontname="Times New Roman")

# plt.savefig('k8_ks4.jpg', dpi=300)
plt.show()

