import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('dist.dat')
x = data[:,0]

def S(x):
    ps = 0.75*np.exp(-np.abs(x-0.25))
    return ps/np.sum(ps)

points = data[:,0]
s = S(points)
y = data[:,1]

plt.plot(points, y, 'k')
plt.plot(points, s,'r', linestyle= "--")


plt.legend(('Ppart(x)','S(x)'),loc='best')
plt.title('Stationary PMF of pparticle position (N=200, BIAS=0.001)')
plt.xlabel('$x$')
plt.ylabel('$Probability$')
plt.grid('both')
plt.show()