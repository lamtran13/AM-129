import numpy as np
from numpy.random import default_rng
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from scipy.integrate import solve_ivp
from scipy.stats import linregress


rng = default_rng()

# Define right hand side of Lorenz equations
def lorenz(t,y):
    f = np.zeros(y.shape)
    f[0] = 10*(y[1] - y[0])
    f[1] = y[0]*(28 - y[2]) - y[1]
    f[2] = y[0]*y[1] - 8*y[2]/3
    return f


# Solve the Lorenz equations for each initial condition
def solveEnsemble(y0_ens,t0,tf):
    t_eval = np.linspace(t0,tf,2000)
    sol_list = []

    for n in range(y0_ens.shape[0]):
        sol = solve_ivp(lorenz,(t0,tf),y0_ens[n, :], t_eval=t_eval, rtol=1.e-8, atol=1.e-10 )
        sol_list.append(sol.y)
    return t_eval,sol_list

# Find mean and standard deviations of the ensemble as functions of time
def ensembleStats(sol_list):
    nEnsemble = len(sol_list)

    # Find mean trajectory
    y_mean = np.zeros_like(sol_list[0])
    for y in sol_list:
        y_mean = y_mean + y
    y_mean = y_mean/nEnsemble

    # Find standard deviation of ensemble
    y_std = np.zeros_like(sol_list[0])
    for y in sol_list:
        y_std = y_std + (y - y_mean)**2
    y_std = np.sqrt(y_std)/(nEnsemble - 1)
    return y_mean,y_std

   
# Generate distances from each pair of solutions
def ensembleDistance(sol_list):
    nEnsemble = len(sol_list)
    dist = []
    
    for i in range(nEnsemble):
        for j in range(i):
            diff = sol_list[i] - sol_list[j]
            dist.append(np.sqrt(np.sum(diff**2,axis=0)))
    return dist

    
# Fit an exponential function and report the maximal LE
# Don't modify this function, it is already complete
def expfit(t,y):
    lr = linregress(t,np.log(y))
    print('Approximate maximal Lyapunov exponent: %f' % lr.slope)
    fit = np.exp(lr.intercept + lr.slope*t)
    return fit

# Plot x/y/z solution components of all members in ensemble
# Don't modify this function, it is already complete
def plotEnsemble(t_eval,sol_list):
    # Get mean and std deviation of ensemble
    y_mean,y_std = ensembleStats(sol_list)

    # 3x2 plots, full ensemble on left, mean/std on right
    fig, axs = plt.subplots(3,2,figsize=(12,10))
    
    # Plot time series of each component
    for y in sol_list:
        for nd in range(3):
            axs[nd,0].plot(t_eval,y[nd,:])
    
    # Plot time series of mean and std dev
    for nd in range(3):
        axs[nd,1].plot(t_eval,y_mean[nd,:],'-b')
        ax_twin = axs[nd,1].twinx()
        ax_twin.semilogy(t_eval,y_std[nd,:],'--r')

    # Annotate plot
    axs[0,0].set_title('Full Ensemble',fontsize=16)
    axs[0,0].set_ylabel('x',fontsize=14); axs[0,0].grid('both')
    axs[1,0].set_ylabel('y',fontsize=14); axs[1,0].grid('both')
    axs[2,0].set_ylabel('z',fontsize=14); axs[2,0].grid('both')
    axs[2,0].set_xlabel('t',fontsize=14)
    
    axs[0,1].set_title('Ensemble Mean and Std. Deviation',fontsize=16)
    axs[0,1].set_ylabel('x',fontsize=14); axs[0,1].grid('both')
    axs[1,1].set_ylabel('y',fontsize=14); axs[1,1].grid('both')
    axs[2,1].set_ylabel('z',fontsize=14); axs[2,1].grid('both')
    axs[2,1].set_xlabel('t',fontsize=14)

    # Small hack to get a legend with twinx axes
    bline = mlines.Line2D([],[],color='blue',label='Mean Solution')
    rline = mlines.Line2D([],[],color='red',linestyle='dashed',label='Std. Deviation')
    axs[2,1].legend(handles=[bline,rline],loc='lower right')

    fig.tight_layout()
    plt.savefig('solutions.png')

# Plot distance between solutions
# Don't modify this function, it is already complete
def plotDistances(t_eval,dist):
    # Get mean and std deviation of ensemble
    y_mean,y_std = ensembleStats(dist)

    # Fit line to first portion of average separation
    t_sub,y_sub = t_eval[:1500],y_mean[:1500]
    fit = expfit(t_sub,y_sub)

    # 1x2 plots, full ensemble on top, mean/std on bottom
    fig, axs = plt.subplots(1,2,figsize=(12,8))
    
    # Plot time series of full ensemble
    for y in dist:
        axs[0].semilogy(t_eval,y)
    
    # Plot time series of mean and std dev
    axs[1].semilogy(t_eval,y_mean[:],'-b',t_sub,fit,'--k')

    # Annotate plot
    axs[0].set_title('Ensemble Distances',fontsize=16)
    axs[0].set_xlabel('t',fontsize=14)
    axs[0].set_ylabel('Distance',fontsize=14)
    axs[0].grid('both')
    
    axs[1].set_title('Mean Distance and Fit',fontsize=16)
    axs[1].set_xlabel('t',fontsize=14)
    axs[1].set_ylabel('Mean',fontsize=14)
    axs[1].grid('both')
    axs[1].legend(['Mean Distance','Exponential Fit'],fontsize=14,loc='lower right')

    fig.tight_layout()
    plt.savefig('distances.png')


def plotDisplacement(sol_list):
    # Form all displacements
    disp = []
    for i in range(nEnsemble):
        for j in range(i):
            disp.append(sol_list[i] - sol_list[j])
    
    ax = plt.figure().add_subplot(projection='3d')

    for n in range(5):
        ax.plot(disp[n][0,:],disp[n][1,:],disp[n][2,:])

    ax.set_xlabel('x',fontsize=15)
    ax.set_ylabel('y',fontsize=15)
    ax.set_zlabel('z',fontsize=15)
    ax.set_title('Ensemble Displacements',fontsize=18)

    plt.show()
    
if __name__ == '__main__':
    # Number of solutions we'll generate
    nEnsemble = 10
    
    # Set initial conditions
    y0_mean = np.array([2,3,-14])
    y0_var = (1.e-9)**2 # std dev is 1e-9, variance is square of that
    y0_cov = np.diag(y0_var*np.ones(3))
    y0_ens = rng.multivariate_normal(y0_mean,y0_cov,nEnsemble)

    # Call solver
    t_eval,sol_list = solveEnsemble(y0_ens,0.0,35.0)

    # Plot time series as ensemble
    plotEnsemble(t_eval,sol_list)

    # Get distances between all solutions
    dist = ensembleDistance(sol_list)
    plotDistances(t_eval,dist)

   # Additional plot for funsies
    plotDisplacement(sol_list)
