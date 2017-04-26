import sys
import numpy as np
from cantera import *
from matplotlib.pylab import *

#################################################################
# Mechanism used for the process
gas = Solution('gri30.xml')
# Enter general parameters
Pmin = 1e4
Pmax = 1e6
npoints = 11
# Specify the number of time steps and the time step size
nt = 100000
dt = 1.e-5  # s
# Storage
# Temperature storage variables
Pi = np.zeros(npoints, 'd')
# The initial storage variable become case dependent
tim = np.zeros(nt, 'd')
temp_cas = np.zeros(nt, 'd')
dtemp_cas = np.zeros(nt - 1, 'd')
# Additional storage variables are needed to differentiate each case
Autoignition_cas = np.zeros(npoints, 'd')
FinalTemp_cas = np.zeros(npoints, 'd')
mfrac_cas = np.zeros([npoints, gas.n_species], 'd')
#################
# Loop over initial conditions
for j in range(npoints):
    Pi[j] = Pmin + (Pmax - Pmin) * j / (npoints - 1)
    # Set gas state, always at stoichiometry
    gas.TPX = 1200, Pi[j], 'C2H6:0.286,O2:1,N2:3.76'
    # Create the batch reactor
    r = IdealGasReactor(gas)
    # Now create a reactor network consisting of the single batch reactor
    sim = ReactorNet([r])
    # Run the simulation
    # Initial simulation time
    time = 0.0
    # Loop for nt time steps of dt seconds.
    # print 'time [s] , T [K] , p [Pa] , u [J/kg]'
    for n in range(nt):
        time += dt
        sim.advance(time)
        tim[n] = time
        temp_cas[n] = r.T
    # for sp in range(m):
    # mfrac_cas[j][sp] = r.massFraction(j)
    # print '%10.3e %10.3f %10.3f %14.6e' % (tim[n], r.temperature(),
    # r.pressure(), r.intEnergy_mass())
    mfrac_cas[j][:] = r.thermo.Y
    # Get autoignition timing
    Dtmax = [0, 0.0]
    for n in range(nt - 1):
        dtemp_cas[n] = (temp_cas[n + 1] - temp_cas[n]) / dt
        if (dtemp_cas[n] > Dtmax[1]):
            Dtmax[0] = n
            Dtmax[1] = dtemp_cas[n]
    # Local print
    Autoignition = (tim[Dtmax[0]] + tim[Dtmax[0] + 1]) / 2.
    print 'For ' + str(Pi[j]) + ', Autoignition time = (s) ' + str(Autoignition)
    # Posterity
    Autoignition_cas[j] = Autoignition * 1000  # ms
    FinalTemp_cas[j] = temp_cas[nt - 1]

# write output CSV file for importing into Excel
# csv_file = 'Phi-1_P-1_Trange_UV.csv'
# with open(csv_file, 'w') as outfile:
#    writer = csv.writer(outfile)
#    writer.writerow(['Auto ignition time','Final Temperature'] + gas.species_names)
#    for i in range(npoints):
#        writer.writerow(['cas Ti = ', Ti[i]])
#        writer.writerow([Autoignition_cas[i], FinalTemp_cas[i]] + list(mfrac_cas[i,:]))
# print 'output written to '+csv_file
# create plot
plot(Pi, Autoignition_cas, '^', color='orange')
xlabel(r'Pressure [Pa]', fontsize=20)
ylabel("Autoignition [ms]")
title(r'Autoignition of $C_{6})H_{6}$ + Air mixture at $\Phi$ = 1, and T=1200K',
      fontsize=22, horizontalalignment='center')
axis([1e4, 1e6, 0.0, 20.0])
grid()
# show()
savefig('Phi-1_T-1500_Trange_UV.png', bbox_inches='tight')