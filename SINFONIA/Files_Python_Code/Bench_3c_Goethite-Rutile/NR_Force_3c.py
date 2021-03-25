"""
SINFONIA, a Python based code for Speciation and Interparticle Forces for Nanoscale Interactions

Copyright (c) 2021
 
Authors: Frank Heberling (frank.heberling@kit.edu) and Teba Gil-Diaz
 
Permission to use and redistribute the source code or binary forms of this software and its documentation, with or without modification is 
hereby granted provided that the above notice of copyright, these terms of use, and the disclaimer of warranty below appear in the source code
and documentation. The names of the authors, or their institutions, may not be used to endorse or promote products derived from this software 
without specific prior written permission from all parties.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
THE USE OR OTHER DEALINGS IN THIS SOFTWARE.
"""
################################################################

import h5py
import numpy as num
import time
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

t0 = time.clock()
################################################################
########### USER INPUT #########################################
filename = 'Benchmark_3c_' ##filename from the library, must be equal to filename defined in "NR_chem_3c.py"

# Choose the geometry for the force-distance calculations:
sphere_sphere = False                # approaching spheres (True)
sphere_plane = True                  # sphere-plane, equivalent to approaching parallel cylinders (True)
plane_plane = False                  # approaching planes (True)
Normalized_force = True              # to obtain Force-distance curves with normalized forces to tip radius

Overview_File = True                 #produce general output file with overview of results

################################################################
########### EXECUTION OF MODEL CALCULATION #####################

# some constants
Temp = 298.15   
teta = (Temp - 273.16)/100
eps_r = 88.15-(41.4+(13.1-4.6*teta)*teta)*teta
eps_0 = 8.85419e-12
ew = eps_r*eps_0
Far = 9.64867e4
kb = 1.38062e-23
e = 1.60219e-19

H = 0.0         # Hamaker constant (J), if 0.0 then no vdW included
r1 = 100        # 1st particle radius in nm 
r2 = 100        # 2nd particle radius in nm

# we need the same information about pH_range and IS as used in the chem_solver, retrieved from HDF5
base = h5py.File(filename+'dictionary.hdf5', 'r')
print base.keys()
pH_range = num.array(base.get('pH_range'))
Initial_IS = num.array(base.get('Initial_IS'))
Distlist = num.array(base.get('Distlist'))
Labels = num.array(base[u'Labels'])
print Labels

# we can have a look at the force-distance curves and save the figure
fig1 = plt.figure(constrained_layout=True)
gs1 = GridSpec(ncols = 1,nrows = len(Initial_IS), figure = fig1)
colors = ['gold', 'orange', 'red', 'sienna', 'olivedrab','lawngreen', 'deeppink', 'turquoise', 'deepskyblue', 'blue', 'blueviolet', 'grey', 'black']
# colors in https://matplotlib.org/3.1.0/gallery/color/named_colors.html


#collecting the data from the library for all IS, pH and distances 
for s in range(len(Initial_IS)):
    CR_list = []
    energylist = []
    all_distlist = []
    all_kappa_dist = []

    for p in range(len(pH_range)):
        osmolist_TmpH = []
        distlist = []
        kappa_dist = []        
        Forcelist_CR = []
        
        for d in range(len(Distlist)):
            #subscript x is used to refer to TiOH and z to FeOH.
            dist = num.array(base.get('IS_'+str(Initial_IS[s])+'/pH_'+str(pH_range[p])+'/distance'+str(Distlist[d])+'/dist'))
            osmo = num.array(base.get('IS_'+str(Initial_IS[s])+'/pH_'+str(pH_range[p])+'/distance'+str(Distlist[d])+'/osmo'))
            kappa_Zhmud = num.array(base.get('IS_'+str(Initial_IS[s])+'/pH_'+str(pH_range[p])+'/distance'+str(Distlist[d])+'/kappa'))
            kappa_distances = num.array(base.get('IS_'+str(Initial_IS[s])+'/pH_'+str(pH_range[p])+'/distance'+str(Distlist[d])+'/KH_2'))
            psi_dx = num.array(base.get('IS_'+str(Initial_IS[s])+'/pH_'+str(pH_range[p])+'/distance'+str(Distlist[d])+'/psidx'))
            psi_dz = num.array(base.get('IS_'+str(Initial_IS[s])+'/pH_'+str(pH_range[p])+'/distance'+str(Distlist[d])+'/psidz'))
            osmolist_TmpH = num.append(osmolist_TmpH, osmo)
            distlist = num.append(distlist, -dist)
            kappa_dist = num.append(kappa_dist, kappa_distances)

            # Calculating interacting forces                                        
            if sphere_sphere:
                Fvdw = (-H/(6*(dist*1e-9)**2))*(r1*r2*1e-9/(r1+r2))                 # vdW force as defined in Table 13.1 from Israelachvili 2011
                if d == 0:
                    Wosmo = osmo
                    force_CR = (2*num.pi*r1*r2*1e-9/(r1+r2))*Wosmo + Fvdw*1e9      # Force in nN
                else:
                    Wosmo = num.trapz(osmolist_TmpH, distlist)
                    force_CR = (2*num.pi*r1*r2*1e-9/(r1+r2))*Wosmo + Fvdw*1e9 

            elif sphere_plane:
                Fvdw = -H*r1/(6*(dist*1e-9)**2)                                     # vdW force as defined in Table 13.1 from Israelachvili 2011
                if d == 0:
                    Wosmo = osmo
                    force_CR = 2*num.pi*r1*1e-9*Wosmo + Fvdw                          #  from Trefalt et al. 2016, force in nN
                else:
                    Wosmo = num.trapz(osmolist_TmpH, distlist)
                    force_CR = 2*num.pi*r1*1e-9*Wosmo + Fvdw   

            elif plane_plane:
                Fvdw = -H/(6*num.pi*(dist*1e-9)**3)                 # in force formulation, Zhao et al. 2015, also Israelachvili Table 13.1
                if d == 0:
                    Wosmo = osmo                                    # it is Fosmo, writen as Wosmo for practical reasons along the script
                    force_CR = num.pi*(r1*1e-9)**2*(Wosmo/dist + Fvdw)*1e9      # Fosmo and Fvdw in Pa, force in nN
                else:
                    Wosmo = num.trapz(osmolist_TmpH, distlist)
                    force_CR = num.pi*(r1*1e-9)**2*(Wosmo/dist + Fvdw)*1e9 
                    
            if Normalized_force:
                Forcelist_CR = num.append(Forcelist_CR, force_CR/r1)        # in nN/nm
            else:
                Forcelist_CR = num.append(Forcelist_CR, force_CR)

            energylist.append(Wosmo*1000)                                   # in mJ/m**2

        #now we produce some output
        CR_list = num.append(CR_list,Forcelist_CR)
        all_distlist = num.append(all_distlist, distlist)
        all_kappa_dist = num.append(all_kappa_dist, kappa_dist)

        #and generate the figure
        a = fig1.add_subplot(gs1[s, 0])
        a.plot(-distlist, Forcelist_CR, color = colors[p], label = str(pH_range[p]))
        a.legend(loc='lower right', title='pH values', ncol=str(len(pH_range)))

        if sphere_sphere:
            a.set_title('IS_'+str(Initial_IS[s]) + '_ss_')
        elif sphere_plane:
            a.set_title('IS_'+str(Initial_IS[s]) + '_sp_')
        elif plane_plane:
            a.set_title('IS_'+str(Initial_IS[s]) + '_pp_')
        else:
            a.set_title('IS_'+str(Initial_IS[s]) + 'No Force')       
        a.set_xlabel('separation (nm)')
        a.set_xlim((0,20))
        if Normalized_force:
            a.set_ylabel('Force (N/m)')
        else:
            a.set_ylabel('Force (nN)')
        a.set_ylim((-0.03,0.03))

    #write overarching results file for all IS, distances and pHs
    if Overview_File:
        if sphere_sphere:
            app1 = '_ss'
        elif sphere_plane:
            app1 = '_sp'
        elif plane_plane:
            app1 = '_pp'
        else:
            app1 = ''
        if Normalized_force:
            app2 = '_norm'
            units = 'nN/nm'
        else:
            app2 = ''
            units = 'nN'
        if H == 0.0:
            app3 = '_NOvdw_'
        else:
            app3 = '_'
        filenamex = filename+str(Initial_IS[s])+'_pH_Forces_'+app1+app2+app3+'general_output.dat'
        g = file((filenamex), 'w')
        g.write(filename+'\n\n')
                
        line = "%12s  %12s  %12s  %12s  %12s  %12s"
        labellist = ['distance','KH/2','pH','ionic_str.', 'interfacial_energy','Force_CR']
        line = line % tuple(labellist)
        line = line +"\n"
        g.write(line)
        line = "%12s  %12s  %12s  %12s  %12s  %12s"
        labellist = ['(nm)', ' ',' ','(mol/L)','mJ/m**2',units]
        line = line % tuple(labellist)
        line = line +"\n\n"
        g.write(line)

        for j in range(len(pH_range)):
            for i in range(len(distlist)):
                k = j*len(distlist)+i
                line = "%.6E  %.6E  %.6E  %.6E  %.6E  %.6E  "
                labellist = [-all_distlist[k],all_kappa_dist[k],pH_range[j],\
                            Initial_IS[s], energylist[k], CR_list[k]]
                line = line % tuple(labellist)
                line = line +"\n" 
                g.write(line)
        g.close()
fig1.savefig('fig_Force.png', dpi = 300, bbox_inches='tight', pad_inches = 0)
plt.show()  
base.close()  
print 'time needed (s): ' + str(time.clock() - t0)
############################### THE END #######################################
