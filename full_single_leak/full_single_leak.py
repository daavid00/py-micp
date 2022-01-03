#Setting up and simulating a 3D flow system with a straight leakage path, where
#CO2 is assessed before and after MICP treatment. Information about the model
#and parameters can be found in:

#Landa-Marbán, D., Kumar, K., Tveit, S., Gasda, S.E. Numerical studies of CO2
#leakage remediation by micp-based plugging technology. In: Røkke, N.A. and
#Knuutila, H.K. (Eds) Short Papers from the 11th International Trondheim CCS
#conference, ISBN: 978-82-536-1714-5, 284-290.

#Set the full path to the flow executable
flow = '~/Opm/master/build/opm-simulators/bin/flow'

#Import python dependencies
import os
import meshio
import numpy as np
from mako.template import Template
from matplotlib import pyplot as plt

#Set the model parameters
rhob    = 35.       #Density (biofilm) [kg/m^3]
rhoc    = 2710.     #Density (calcite) [kg/m^3]
rhow    = 1045      #Density (water) [kg/m^3]
kstr    = 2.6e-10   #Detachment rate [m/(Pa s)]
crit    = 0.1       #Critical porosity [-]
eta     = 3.0       #Fitting factor [-]
ko      = 2e-5      #Half-velocity coefficient (oxygen) [kg/m^3]
ku      = 21.3      #Half-velocity constant (urea) [kg/m^3]
mu      = 4.17e-5   #Maximum specific growth rate [1/s]
muu     = 0.0161    #Maximum rate of urease utilization [1/s]
ka      = 8.51e-7   #Microbial attachment rate [1/s]
kd      = 3.18e-7   #Microbial death rate [1/s]
kmin    = 1e-20     #Minimum permeability [m^2]
F       = 0.5       #Oxygen consumption factor [-]
Y       = 0.5       #Yield growth coefficient [-]
muw     = 2.54e-4   #Water viscoscity [Pa s]
Ka      = 1e-14     #Aquifer permeability [m^2]
Kl      = 2e-14     #Leakage permeability [m^2]
phi     = 0.15      #Porosity [-]
md      = 9.87e-16  #Milli darcy [mD]
day     = 86400     #Seconds in a day [s]
tolclg  = 1e-4      #Tolerance before clogging to stop the simulation
QCO2    = 4e-4      #Injection rate (CO2) [m^3/s]
dtCO2   = 10        #Time step to print CO2 results [d]
NCO2    = 50        #Number of time steps to run the simulation [-]

#Delete previous inputs and results
os.system('rm -rf inputs & rm -rf results & wait')
os.system('mkdir inputs & mkdir results & wait')

#Define the porosity-permeability functions (aquifer and leakage)
def perma(bio, cal):
     return (Ka * ((phi - bio - cal - crit) / (phi - crit)) ** eta + kmin) \
                                                       * Ka / (md * (Ka + kmin))
def perml(bio, cal):
     return (Kl * ((phi - bio - cal - crit) / (phi - crit)) ** eta + kmin) \
                                                       * Kl / (md * (Kl + kmin))

#Create the grid
H  = 30                             #Height of the domain [m]
ht = 5                              #Height of the top aquifer [m]
hl = 5                              #Height of the lower aquifer [m]
L  = 100                            #Length of the domain [m]
ux = 14                             #X-gap between aquifer and leakage [m]
uy = 14                             #Y-gap between aquifer and leakage [m]
Wi = 100                            #Width of the domain [m]
X0 = np.linspace(0, -1.5, 16)       #""                         -L*exp(X0)
X1 = np.linspace(-21, 21, 44)       #Discretization of the x-dir X1
X2 = np.linspace(-1.5, 0, 16)       #""                          L*exp(X2)
Y0 = np.linspace(0, -1.5, 16)       #""                         -Wi*exp(Y0)
Y1 = np.linspace(-21, 21, 44)       #Discretization of the y-dir Y1
Y2 = np.linspace(-1.5, 0, 16)       #""                          Wi*exp(Y2)
Z0 = np.linspace(0, ht, ht + 1)     #Discretization of the z-dir Z0
Z1 = np.linspace(ht + 1.0, H - hl + 1., H - hl - ht + 1) #""     Z1
Z2 = np.linspace(H - hl + 2., H, 4) #""                          Z2
Wc = 1                              #Cells for each MICP well injection
XX = []; YY = []; ZZ = []
for i in range(len(X0)):
    XX.append(-L * np.exp(X0[i]))
for i in range(len(X1)):
    XX.append(X1[i])
for i in range(len(X2)):
    XX.append(L * np.exp(X2[i]))
for j in range(len(Y0)):
    YY.append(-Wi * np.exp(Y0[j]))
for j in range(len(Y1)):
    YY.append(Y1[j])
for j in range(len(Y2)):
    YY.append(Wi * np.exp(Y2[j]))
for k in range(len(Z0)):
    ZZ.append(Z0[k])
for k in range(len(Z1)):
    ZZ.append(Z1[k])
for k in range(len(Z2)):
    ZZ.append(Z2[k])
xx = len(X0) + len(X1) + len(X2) - 1
yy = len(Y0) + len(Y1) + len(Y2) - 1
zz = len(Z0) + len(Z1) + len(Z2) - 1
xyz= xx * yy * zz
ix = abs(np.transpose(XX)).argmin() + 1
iy = abs(np.transpose(YY)).argmin() + 1
it = len(np.transpose(ZZ)[tuple([np.transpose(ZZ) < ht])])
il = len(np.transpose(ZZ)[tuple([np.transpose(ZZ) > H - hl])])
os.system(f"octave --eval 'full_single_leak({H}, {ht}, {hl}, {L}, {ux}, " + \
                                             f"{uy}, {Wi}, {XX}, {YY}, {ZZ})' ")

#Create variables used in this script
biof = [0] * xyz ; calc = [0] * xyz ; injestra = []; perm = [Ka / md] * xyz
Lg = []; t0 = []; t1 = []; ACTNUM = []; VOLUM = []

#Set and run the CO2 simulation before MICP treatment
perm[xx * yy * it  : xx * yy * (zz - il)] = [Kl / md] * \
                                                      (xx * yy * (zz - il - it))
mytemplate = Template(filename = 'co2.mako')
var = {'biof' : biof, 'calc' : calc, 'dtCO2' : dtCO2, 'il' : il, 'ix' : ix, \
              'it' : it, 'iy' : iy, 'NCO2' : NCO2, 'perm' : perm, 'phi' : phi, \
                                'QCO2' : QCO2, 'xx' : xx, 'yy' : yy , 'zz' : zz}
FilledTemplate = mytemplate.render(**var)
with open('inputs/co2.data', 'w') as f:
    f.write(FilledTemplate)
os.system(f"{flow} inputs/co2.data --output-dir=results " + \
        " --enable-vtk-output=true --initial-time-step-in-days=0.0001 & wait\n")

#Define injection strategy
#[injection time [h], injection rate [m^3/s], microbial concentration [kg/m^3],
#oxygen concentration [kg/m^3], and urea concentration [kg/m^3]]. Before each
#injection of the solutions in the phases we inject only water for a very short
#time [0.01 hour] to ease the convergence of the solution.
injestra.append([0.01, 8e-2, 0, 0, 0])
injestra.append([15., 8e-2, 0.01, 0, 0])
injestra.append([2., 8e-2, 0, 0, 0])
injestra.append([100., 0, 0, 0, 0])
injestra.append([0.01, 8e-2, 0, 0, 0])
injestra.append([30., 8e-2, 0, 0.04, 0])
injestra.append([2., 8e-2, 0, 0, 0])
injestra.append([30., 0, 0, 0, 0])
injestra.append([0.01, 8e-2, 0, 0, 0])
injestra.append([30., 8e-2, 0, 0, 60.])
injestra.append([2., 8e-2, 0, 0, 0])
injestra.append([30., 0, 0, 0, 0])
injestra.append([0.01, 8e-2, 0, 0, 0])
injestra.append([30., 8e-2, 0, 0.04, 60.])
injestra.append([2., 8e-2, 0, 0, 0])
injestra.append([30., 0, 0, 0, 0])
injestra.append([0.01, 8e-2, 0, 0, 0])
injestra.append([30., 8e-2, 0, 0.04, 60.])
injestra.append([2., 8e-2, 0, 0, 0])
injestra.append([30., 0, 0, 0, 0])
injestra.append([0.01, 8e-2, 0, 0, 0])
injestra.append([30., 8e-2, 0, 0.04, 60.])
injestra.append([2., 8e-2, 0, 0, 0])
injestra.append([30., 0, 0, 0, 0])

# We store the maximum injected oxygen and urea concentrations. They are use
# so that the computed oxygen and urea values are within these during the
# solution step
comax = 0; cumax = 0
for i in range(len(injestra)):
    comax = max(injestra[i][3], comax)
    cumax = max(injestra[i][4], cumax)

#Read the input file
mytemplate = Template(filename = 'micp.mako')

#Set and run the MICP simulation
var = {'biof' : biof, 'calc' : calc, 'crit' : crit, 'comax' : comax, 'cumax' : \
    cumax, 'eta' : eta, 'F' : F, 'il' : il, 'ix' : ix, 'iy' : iy, 'injestra' : \
      injestra, 'ka' : ka, 'kd' : kd, 'kmin' : kmin, 'ko' : ko, 'kstr' : kstr, \
                'ku' : ku, 'mu' : mu, 'muu' : muu, 'muw' : muw, 'perm' : perm, \
                     'phi' : phi, 'rhob' : rhob, 'rhoc' : rhoc, 'rhow' : rhow, \
                                      'tolclg' : tolclg, 'Wc' : Wc, 'xx' : xx, \
                                                  'Y' : Y, 'yy' : yy, 'zz' : zz}
FilledTemplate = mytemplate.render(**var)
with open('inputs/micp.data', 'w') as f:
    f.write(FilledTemplate)
os.system(f"{flow} inputs/micp.data --output-dir=results" + \
   " --initial-time-step-in-days=0.0001 --solver-max-time-step-in-days=0.01" + \
                                           " --enable-vtk-output=true & wait\n")

#Update the porosity and permeability after MICP treatment
with open("inputs/full_single_leak.grdecl", "r") as f:
    list_of_lines = f.readlines()
act = list_of_lines[-int(xyz / 20.) -1 : ]
for k in range(len(act)):
    for i in range(len(act[k])):
      if act[k][i] == '0' or act[k][i] == '1':
          ACTNUM.append(int(act[k][i]))
with open("results/MICP.pvd", "r") as f:
    list_of_lines = f.readlines()
ss = list_of_lines[-3]
mesh = meshio.read(f"results/MICP-{int(ss[-13:-8]):05d}.vtu")
for row in mesh.cell_data['calcite fraction']:
    c = row
for row in mesh.cell_data['biofilm fraction']:
    b = row
biof[ : xx * yy * it] = b[ : xx * yy * it]
calc[ : xx * yy * it] = c[ : xx * yy * it]
perm[ : xx * yy * it] = perma(b[ : xx * yy * it], c[ : xx * yy * it])
j = 0
for k in range(it, zz - il):
    for l in range(0, xx * yy):
        if ACTNUM[xx * yy * k + l] == 1:
            biof[xx * yy * k + l] = b[xx * yy * it + j]
            calc[xx * yy * k + l] = c[xx * yy * it + j]
            perm[xx * yy * k + l] = perml(b[xx * yy * it + j], \
                                                            c[xx * yy * it + j])
            j = j + 1
i = xx * yy * zz - xx * yy * il
I = len(b) - xx * yy * il
biof[i : ] = b[I : ]
calc[i : ] = c[I : ]
perm[i : ] = perma(b[I : ], c[I : ])

#Set and run the CO2 simulation after MICP treatment
mytemplate = Template(filename = 'co2.mako')
var = {'biof' : biof, 'calc' : calc, 'dtCO2' : dtCO2, 'il' : il, 'ix' : ix, \
              'it' : it, 'iy' : iy, 'NCO2' : NCO2, 'perm' : perm, 'phi' : phi, \
                                'QCO2' : QCO2, 'xx' : xx, 'yy' : yy , 'zz' : zz}
FilledTemplate = mytemplate.render(**var)
with open('inputs/co2micp.data', 'w') as f:
    f.write(FilledTemplate)
os.system(f"{flow} inputs/co2micp.data --output-dir=results" + \
        " --enable-vtk-output=true --initial-time-step-in-days=0.0001 & wait\n")

#Obtain the times for the .vtu results and save them in 't' in days
with open("results/CO2.pvd", "r") as f:
    list_of_lines = f.readlines()
nn = 0; k = 0; T = []
for i in range(6, len(list_of_lines)):
    a = []
    for j in range(22, len(list_of_lines[i])):
        a.append(list_of_lines[i][j])
        if list_of_lines[i][j + 1] == '"':
            T.append(a)
            break
for k in range(len(T)):
    aaa = str(T[k]).lstrip('[').rstrip(']')
    bbb = aaa.replace(', ',"")
    xxx = 0; n = 0
    for i in range(len(bbb)):
        if bbb[i] != "'":
            if bbb[len(bbb) - i - 1] == ".":
                xxx = xxx * pow(10, -n)
                n = 0
            else:
                xxx = xxx + pow(10, n) * float(bbb[len(bbb) - i - 1])
                n = n + 1
    t0.append(xxx / day)
with open("results/CO2MICP.pvd", "r") as f:
    list_of_lines = f.readlines()
nn = 0; k = 0; T = []
for i in range(6, len(list_of_lines)):
    a = []
    for j in range(22, len(list_of_lines[i])):
        a.append(list_of_lines[i][j])
        if list_of_lines[i][j + 1] == '"':
            T.append(a)
            break
for k in range(len(T)):
    aaa = str(T[k]).lstrip('[').rstrip(']')
    bbb = aaa.replace(', ',"")
    xxx = 0; n = 0
    for i in range(len(bbb)):
        if bbb[i] != "'":
            if bbb[len(bbb) - i - 1] == ".":
                xxx = xxx * pow(10, -n)
                n = 0
            else:
                xxx = xxx + pow(10, n) * float(bbb[len(bbb) - i - 1])
                n = n + 1
    t1.append(xxx / day)

#Compute the mass percentage of leaked CO2 to the upper aquifer before and after
#MICP treatment
for k in range(it):
    z = ZZ[k + 1] - ZZ[k]
    for j in range(yy):
        y = YY[j + 1] - YY[j]
        for i in range(xx):
            VOLUM.append([(XX[i + 1] - XX[i]) * y * z])
a = []
II = xx * yy * it - 1
for j in range(len(t0)):
    mesh = meshio.read(f"results/CO2-{j:05d}.vtu")
    for row in mesh.cell_data['saturation_gas']:
        CO2 = row
    for row in mesh.cell_data['porosity']:
        p = row
    CO2M = QCO2 * t0[j] * day
    if t0[j] > 0:
        a.append(100. * sum(CO2[: II]  * p[: II] * VOLUM[: II]) / CO2M)
    else:
        a.append(sum(CO2[: II]  * p[: II] * VOLUM[: II]))
Lg.append(a)
a = []
for j in range(len(t1)):
    mesh = meshio.read(f"results/CO2MICP-{j:05d}.vtu")
    for row in mesh.cell_data['saturation_gas']:
        CO2 = row
    for row in mesh.cell_data['porosity']:
        p = row
    CO2M = QCO2 * t1[j] * day
    if t1[j] > 0:
        a.append(100. * sum(CO2[: II]  * p[: II] * VOLUM[: II]) / CO2M)
    else:
        a.append(sum(CO2[: II]  * p[: II] * VOLUM[: II]))
Lg.append(a)
print(Lg[0][-1])
print(Lg[1][-1])

#Plot the results over time and save them as co2mass_comparison.eps
lw = 5
plt.figure(figsize = (5, 4), dpi = 256)
plt.rc('font', size = 9)
plt.rc('legend', fontsize = 9)
axes = plt.subplot(1, 1, 1)
plt.plot(t0, Lg[0][:], \
 color = [1, 0.2, 0.2], linewidth = lw, linestyle = "-", label = "Without MICP")
plt.plot(t1, Lg[1][:], \
     color = [1, 0.5, 0], linewidth = lw, linestyle = "-", label = "After MICP")
plt.xlabel('t [days]')
plt.ylabel('Leaked $CO_{2}$ mass/injected $CO_{2}$ mass [%]')
plt.grid()
plt.legend(loc = 'upper left')
plt.savefig('results/co2mass_comparison.eps', format='eps')
plt.show()

#{
#Copyright 2021, NORCE Norwegian Research Centre AS, Computational
#Geosciences and Modeling.

#This file is part of the py-micp module.

#py-micp is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#ad-wa is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this file.  If not, see <http://www.gnu.org/licenses/>.
#}
