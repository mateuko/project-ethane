import sys
import numpy as np
import csv
import math
from cantera import *
from matplotlib.pylab import *
#################################################################
#Mechanism used for the process
gas = Solution('gri30.xml')
#Enter general parameters
Tmin = 0.8 # Kelvin
Tmax = 1.0 # Kelvin
Pmin = 0.5e5
Pmax = 1e6
npoints = 6
#Specify the number of time steps and the time step size
nt = 100000
dt = 1.e-5 #s
#Storage
#Temperature storage variables
Ti = np.zeros(npoints,'d')
Ti2 = np.zeros(npoints,'d')
Pi = np.zeros(npoints,'d')
#The initial storage variable become case dependent
tim = np.zeros(nt,'d')
temp_cas = np.zeros((npoints,npoints,nt),'d')
dtemp_cas = np.zeros(nt-1,'d')
#Additional storage variables are needed to differentiate each case
Autoignition_cas = np.zeros((npoints,npoints),'d')
FinalTemp_cas = np.zeros((npoints,npoints),'d')
mfrac_cas = np.zeros([npoints,gas.n_species],'d')
tt=[]
#################
#Loop over initial conditions
for k in range(npoints):
    Pi[k]=Pmin+(Pmax-Pmin)*k/(npoints-1)
    for j in range(npoints):
        Ti2[j] = Tmin + (Tmax - Tmin)*j/(npoints - 1)
        Ti[j] = 1000/Ti2[j]
    #Set gas state, always at stoichiometry
        gas.TPX = Ti[j], Pi[k], 'C2H6:0.286,O2:1,N2:3.76'
    #Create the batch reactor
        r = IdealGasReactor(gas)
    # Now create a reactor network consisting of the single batch reactor
        sim = ReactorNet([r])
    #Run the simulation
    # Initial simulation time
        time = 0.0
    #Loop for nt time steps of dt seconds.
        for n in range(nt):
            time += dt
            sim.advance(time)
            tim[n] = time
            temp_cas[k][j][n] = r.T
        mfrac_cas[j][:] = r.thermo.Y
        # Get autoignition timing
        Dtmax = [0,0.0]
        for n in range(nt-1):
            dtemp_cas[n] = (temp_cas[k][j][n+1]-temp_cas[k][j][n])/dt
            if (dtemp_cas[n]>Dtmax[1]):
                Dtmax[0] = n
                Dtmax[1] = dtemp_cas[n]
        # Local print
        Autoignition = (tim[Dtmax[0]]+tim[Dtmax[0] + 1])/2.
        print 'For '+str(Pi[k])+'Pa and' +str(Ti[j]) +'K, Autoignition time = (s) ' + str(Autoignition)
        # Posterity
        Autoignition_cas[k][j] = Autoignition*1000
        FinalTemp_cas[k][j] = temp_cas[k][j][nt-1]

# write output CSV file
csv_file = 'T_P_t.csv'
with open(csv_file, 'w') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(['(Ti, Pi, Autoignition_time'])
    for k in range(npoints):
        for i in range(npoints):
            writer.writerow(['(' +str(Ti[i]), Pi[k], str(Autoignition_cas[k][i])+ ')'])

time=0.0
for n in range(nt):
    time += dt
    tt.append(time)

for i in range(npoints):
    for j in range(npoints):
        plot(tt, temp_cas[i][j], '-', color='orange')

    xlabel(r'Time [s]', fontsize=20)
    ylabel("Temperature [K]")
    title(r'Autoignition of $C_{2}H_{6}$ + Air mixture at $\Phi$ = 1, and Pi='+str(Pi[i])+'Pa',
          fontsize=22, horizontalalignment='center')
    axis([0, 0.5, 1000.0, 3000.0])
    grid()
    # show()
    name='P'+str(i)+'_T_Trange_UV.png'
    savefig(name, bbox_inches='tight')
    clf()

for i in range(npoints):
    for j in range(npoints):
        plot(tt, temp_cas[j][i], '-', color='orange')

    xlabel(r'Time [s]', fontsize=20)
    ylabel("Temperature [K]")
    title(r'Autoignition of $C_{2}H_{6}$ + Air mixture at $\Phi$ = 1, and Ti='+str(Ti[i])+'K',
          fontsize=22, horizontalalignment='center')
    axis([0, 0.5, 1000.0, 3000.0])
    grid()
    # show()
    name='T'+str(i)+'_T_Trange_UV.png'
    savefig(name, bbox_inches='tight')
    clf()

for i in range(npoints):
    for j in range(npoints):
        plot(Pi[j], Autoignition_cas[j][i], 'o', color='red')
    xlabel(r'Initial pressure [Pa]', fontsize=20)
    ylabel("Autoignition time [ms]")
    title(r'Autoignition of $C_{2}H_{6}$ + Air mixture at $\Phi$ = 1, and Ti='+str(Ti[i])+'K',
          fontsize=22, horizontalalignment='center')
    axis([0, 1.1e6, 0.0, 500])
    grid()
    # show()
    name='T'+str(i)+'_Autoignition.png'
    savefig(name, bbox_inches='tight')
    clf()

for i in range(npoints):
    plot(Ti, Autoignition_cas[i][range(npoints)], 'o', color='red')
    xlabel(r'Initial temperature [K]', fontsize=20)
    ylabel("Autoignition time [ms]")
    title(r'Autoignition of $C_{2}H_{6}$ + Air mixture at $\Phi$ = 1, and Pi='+str(Pi[i])+'K',
          fontsize=22, horizontalalignment='center')
    axis([900, 1300, 0.0, 500])
    grid()
    # show()
    name='P'+str(i)+'_Autoignition.png'
    savefig(name, bbox_inches='tight')
    clf()