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

N=2
k=np.sqrt(2)
k_s=1
xmax=20
node=2000

z0=[5*math.pi/16,5*math.pi/8]
# z0=[math.pi/16,5*math.pi/8]
# z0=[math.pi/16,19*math.pi/32]
# z0=[3*math.pi/8,9*math.pi/16]
# print(z0)

def diff_equation(t,y):
    dy=np.zeros(shape=(N))
    for i in list(range(N)):
        dy[i]=k/N*(np.sum(np.sin(y-y[i]),axis=0))-k_s*np.sin(2*y[i])
    return np.array(dy)
t=np.linspace(0,xmax,num=node)
result=solve_ivp(diff_equation,t_span=(0,xmax),y0=z0,t_eval=t,method='RK45')
result = result.y.T
# result=result.round(6)
print("result",result[-1,:])

config = {
    "font.family": 'Times New Roman',
    "font.size": 20,
    "mathtext.fontset": 'stix',
    "font.serif": ['simhei'],
}
rcParams.update(config)

fig2, ax1 = plt.subplots(figsize=(9, 7))
#plt.axes(yscale = "log")
#plt.plot(t, m1, label=r'$m(\xi^{1})$', color='#00BFFF', linewidth=1.5, linestyle='-')
#plt.plot(t, m2, label=r'$m(\xi^{2})$', color='r', linewidth=1.5, linestyle='-')

for i in list(range(N)):
    plt.plot(t, result[:,i],linewidth=1.5, linestyle='-',label=r'$\varphi$'+str(i+1))
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
#ax1.set_xlabel(r'$m$', fontsize=23)
#ax1.set_ylabel(r'$t$', fontsize=23)

ax1.legend(loc='best', fontsize=18)
plt.xlim([0,xmax])
plt.ylim([-0.1,3.24])

x_major_locator=MultipleLocator(5)
y_major_locator=MultipleLocator(0.5)

ax=plt.gca()
ax.xaxis.set_major_locator(x_major_locator)
ax.yaxis.set_major_locator(y_major_locator)
# plt.savefig('fig1.png')
plt.show()




