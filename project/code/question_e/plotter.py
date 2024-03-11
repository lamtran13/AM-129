import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('sol_n8.dat')
data16 = np.loadtxt('sol_n16.dat')
data32 = np.loadtxt('sol_n32.dat')

t = data[:,0]
Tf = len(data[:,0])
timesteps = [np.ceil(Tf/4), np.ceil(Tf/2), np.ceil((3/4) * Tf), Tf-1]

x = []
x16 = [] 
x32 = []
for i in timesteps:
    x.append(data[np.int_(i),4])
    x16.append(data16[np.int_(i),4])
    x32.append(data32[np.int_(i),4])

fig, ax = plt.subplots()

plt.plot(timesteps,x,'-k.',label='N=8', linewidth=1)
plt.plot(timesteps,x16,'-m.',label='N=16',linewidth=1)
plt.plot(timesteps,x32,'-r.',label='N=32',linewidth=1)
plt.title('Question e')

ax.legend(loc='upper right')

plt.figure()
plt.plot(data16[:,0], data16[:,8])
plt.title('figure 2')

plt.show()

