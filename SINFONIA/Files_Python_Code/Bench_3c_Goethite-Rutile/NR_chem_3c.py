"""
SINFONIA, a Python based code for Speciation and Interparticle Forces for Nanoscale Interactions

Copyright (c) 2021
 
Authors: Frank Heberling (frank.heberling@kit.edu) and Teba Gil-Diaz
 
Permission to use and redistribute the source code or binary forms of this software and its documentation, with or without modification is 
hereby granted provided that the above notice of copyright, these terms of use, and the disclaimer of warranty below appear in the source code
and documentation. The names of the authors, or their institutions, may not be used to endorse or promote products derived from this software 
without specific prior written permission from all parties.
If the code, or parts of it are published or used as a part of a scientific publication, proper credit to the original publication: Gil-Diaz et al. (2021) must be given.
 
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
THE USE OR OTHER DEALINGS IN THIS SOFTWARE.
"""
################################################################

import h5py
import numpy as num
import time
from matplotlib import pyplot
from NR_solver_3c import calc_IS, simplex_point, NR_fit


t0 = time.clock()
################################################################
########### USER INPUT #########################################
filename = 'Benchmark_3c_' #filename for the ouput hdf library

pH_range = num.array([3.8, 5.0, 6.3], float)
#list of pH values at which calculations are performed

Initial_IS = [0.001]    # define working ionic strength (M), will be used as total Na+ and Cl- concentration in this example

#pH corrections for e.g. T1m = Na+, T2m = Cl-. Rows are in range(len(Initial_IS)) and columns in pH_range
#For Benchmark 3c no pH corrections were required.
#otherwise the data represent acid (HCl) or base (NaOH) additions needed to adjust the pH to values specified above.
pH_corr_T1m = num.array([[0.0, 0.0, 0.0]], float)
pH_corr_T2m = num.array([[0.0, 0.0, 0.0]], float)

#defining distances between surfaces (nm)
D = num.concatenate((num.linspace(50,30,5),num.linspace(29,10,10),num.linspace(9,0.1,50),[0.08,0.05,0.01]))

# Activate calculations with electrostatic components and conversion of concentrations (used in Kc definitions) into activities (required for Kint) 
Activities = True                   # if True = calculations are performed in activities

# Create library
base = h5py.File(filename+'dictionary.hdf5', 'w')
base.create_dataset('pH_range', data=pH_range)
base.create_dataset('Initial_IS', data=Initial_IS)
# pH_range and IS_initial are required as inputs for the NR_Force_3c.py thus added to the library as independent datasets

################################################################
########### MODEL DEFINITIONS ##################################
tr = 1.0e-13                        # target residual, implemented as a percentage of each component
Bolz = 1                            # Bolzano's approach used in the Newton-Raphson calculation (>=1)

Labels = ['H+','Na+','Cl-','OH-','>TiOH-0.5','>TiOH2+0.5','>FeOH-0.5','>FeOH2+0.5']  #Defining the species involved in the chemical system

Zspec = num.array([1, 1, -1, -1, 0, 0, 0, 0], float)  # Zspec: Vector with solution species charges 
#note: surface species have no activity and therefore no charge

K = num.array([0.0, 0.0, 0.0, -14, 0.0, 5.8, 0.0, 9.5], float)   #K: Vector with log10_k values for reactions considered in LMA calculations
#K: K=(H+, Na+, Cl-, OH-, >TiOH, >TiO-,>FeOH, >FeO-)

#cap: array with capacitances for each surface (F/m**2)
cap_Ti = num.array([1.33, 1e4, 1e4],float)
cap_Fe = num.array([1.10, 1e4, 1e4],float)

#surf: surface definition in sites/nm**2, m**2/g, g/L.
surf_Ti = num.array([12.2, 1.0, 1.0],float)
surf_Fe = num.array([6.15, 1.0, 1.0],float)

#A and B Matrices, defining stoichiometries of reactions (A) and charge balances (B)
#H+ = 0 + H+
#Na+ = 0 + Na+
#Cl- = 0 + Cl-
#OH- = -14 -H+
#TiOH = TiOH
#TiOH2+ = H+ + TiOH
#FeOH = FeOH
#FeOH2+ = H+ + FeOH

A = num.array([\
 [ 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],\
 [ 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],\
 [ 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],\
 [-1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],\
 # surface species definition
 [ 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],\
 [ 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],\
 [ 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],\
 [ 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],\
 ], float)

B = num.array([\
 [ 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],\
 [ 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],\
 [ 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],\
 [-1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],\
 # surface species definition
 [ 0.0, 0.0, 0.0, 1.0, 0.0, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],\
 [ 1.0, 0.0, 0.0, 1.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],\
 [ 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, -0.5, 0.0, 0.0, 0.0],\
 [ 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0],\
 ], float)

################################################################
########### EXECUTION OF MODEL CALCULATIONS ####################

# some constants
Temp = 298.15
teta = (Temp - 273.16)/100
eps_r = 88.15-(41.4+(13.1-4.6*teta)*teta)*teta
eps_0 = 8.85419e-12
ew = eps_r*eps_0
e = 1.60219e-19
Far = 9.64853e4
kb = 1.38062e-23
R = 8.314463

#loop through ionic strengths (in this example only one value is given)
for s in range(len(Initial_IS)):
    
    grp1=base.create_group('IS_'+str(Initial_IS[s]))    #starting the dictionary to collect data per ionic strength
    
    T1 = Initial_IS[s]       #general concentration of component 2 (e.g. Na+)
    T2 = Initial_IS[s]       #general concentration of component 3 (e.g. Cl-)

    pH_corr_T1 = pH_corr_T1m[s,:]     #collecting pH corrections for component 2  
    pH_corr_T2 = pH_corr_T2m[s,:]     #collecting pH corrections for component 3 

    #loop through pH values
    for pH in range(len(pH_range)):
        sbgrp1 = grp1.create_group('pH_'+str(pH_range[pH]))
        
        # T: Vector with the total concentrations of components
        # T: T = (H, Na, Cl, TiOH, FeOH); total H is defined as activity because pH will be fixed
        T = num.array([10**(-pH_range[pH]), T1+pH_corr_T1[pH], T2+pH_corr_T2[pH], 0.0, 0.0], float)
        Zcomp = num.array([1, 1, -1, 0, 0], float) #charges of components used to calculate activity corrections

        T[-2] = surf_Ti[0]/6.02214086e5 * surf_Ti[1] * surf_Ti[2] #set the total surface site concentration based on surface definition given above 
        T[-1] = surf_Fe[0]/6.02214086e5 * surf_Fe[1] * surf_Fe[2] #set the total surface site concentration based on surface definition given above 
            
        # X: Starting values for log10 master-species concentrations set to the corresponding total component concentrations, H+ is always first with fixed activity
        X = num.log10(T)
        #append 8 zeros to T and X, starting values for surface charges and potentials in FLM for respective surfaces (TiOH and FeOH)
        T = num.append(T, [0.,0.,0.,0.,0.,0.,0.,0.])
        X = num.append(X, [0.,0.,0.,0.,0.,0.,0.,0.])
        Zcomp = num.append(Zcomp, [0.,0.,0.,0.,0.,0.,0.,0.])

        # IS ionic strength
        IS = calc_IS(Zcomp, X) # initial guess for the ionic strength

        #loop through distances defined in nm
        for d in range(len(D)):
            dist = D[d]
            result = simplex_point(K,A,B,T,Zcomp,Zspec,IS,X,cap_Ti,cap_Fe,surf_Ti, surf_Fe, dist, eps_r, eps_0, Temp, kb, e, Activities = Activities)
            result.calc_MB(False) #starting point for equilibrium calculation

            result = NR_fit(tr, Bolz, result, True, eps_r, eps_0) #start Newton-Raphson fitting routine
            #param 1: target residual

            #calculating the corresponding KH/2
            kappa_Zhmud = 1e-9*((2*e*Far*1000*result.IS)/(ew*kb*Temp))**0.5
            kappa_dist = (kappa_Zhmud*dist)/2
            
            #grab all output from the final result
            psi, feld, ibal, sig_diff_bvp, osmo = result.calc_EB(True)                                           
                
            #now we produce some output
            activities = 10**(result.C)
            concentrations = 10**(result.C-result.G2)
            charge = num.sum(result.Zspec*concentrations)
            
            #first print with summary of each condition to the shell
            print "\nOUTPUT: \n"
            print "original totals were: " + str( result.T[1:-8] )
            print "Resulting totals are: "+ str( result.Ttmp[1:-8])+"\n"
            print "mass balance residuals are: "+ str(num.fabs(result.Y[1:-8]))
            print "electrostatic balance residuals are: " + str(num.fabs(result.Y[-8:]))
            
            print "electrical balance in bulk solution (mol(eq)): " + str(charge)
            print "electric balance in EDL (C/m**2): " +str((num.sum(result.sig[0:3]*result.fact2_xy))*(num.sum(result.sig[-4:-1]*result.fact2_z))+sig_diff_bvp) +'\n'
            print 'osmotic pressure (Pa): ' + str(round(osmo,2))
            print 'psi_0xy (V): ' + str(round(result.psi[0],5))
            print 'psi_1xy (V): ' + str(round(result.psi[1],5))
            print 'psi_2xy (V): ' + str(round(result.psi[2],5))
            print 'psi_dxy (V): ' + str(round(result.psi[3],5))
            print 'psi_0z (V): ' + str(round(result.psi[4],5))
            print 'psi_1z (V): ' + str(round(result.psi[5],5))
            print 'psi_2z (V): ' + str(round(result.psi[6],5))
            print 'psi_dz (V): ' + str(round(result.psi[7],5))
            print 'sig_0xy (mol, eq): ' + str(result.sig[0])
            print 'sig_1xy (mol, eq): ' + str(result.sig[1])
            print 'sig_2xy (mol, eq): ' + str(result.sig[2])
            print 'sig_dxy (mol, eq): ' + str(result.sig[3])
            print 'sig_0z (mol, eq): ' + str(result.sig[4])
            print 'sig_1z (mol, eq): ' + str(result.sig[5])
            print 'sig_2z (mol, eq): ' + str(result.sig[6])
            print 'sig_dz (mol, eq): ' + str(result.sig[7])
            print 'sig_0xy (C/m**2): ' + str(result.sig[0]*result.fact2_xy)
            print 'sig_1xy (C/m**2): ' + str(result.sig[1]*result.fact2_xy)
            print 'sig_2xy (C/m**2): ' + str(result.sig[2]*result.fact2_xy)
            print 'sig_dxy (C/m**2): ' + str(result.sig[3]*result.fact2_xy)
            print 'sig_0z (C/m**2): ' + str(result.sig[4]*result.fact2_z)
            print 'sig_1z (C/m**2): ' + str(result.sig[5]*result.fact2_z)
            print 'sig_2z (C/m**2): ' + str(result.sig[6]*result.fact2_z)
            print 'sig_dz (C/m**2): ' + str(result.sig[7]*result.fact2_z)+"\n"

            print "for pH : " + str(pH_range[pH]) + " and separation : " + str(round(dist,2)) + " nm"
            print "ionic strength is: " + str( result.IS )
            print "species concentrations, activities, and log10(activity coefficients): \n"
            for i in range(len(Labels)):
                print "%10s   %.3E   %.3E  %.3E" % (Labels[i], concentrations[i], activities[i], result.G2[i]) 
            print '\n'
            
            x = num.linspace(0,dist,1001)

            #set starting conditions for next calculation
            X = result.X
            IS = result.IS
            
            # fill in the library with the collected output: 
            sbbgrp1 = sbgrp1.create_group('distance'+str(D[d]))
            sbbgrp1.create_dataset('dist',data=[dist])
            sbbgrp1.create_dataset('KH_2',data=[kappa_dist])
            sbbgrp1.create_dataset('kappa',data=[kappa_Zhmud])
            sbbgrp1.create_dataset('ionic_str',data=[result.IS])
            sbbgrp1.create_dataset('sum_resid',data=[result.MB])            
            sbbgrp1.create_dataset('ch_bal_aq',data=[charge])

            for c in range(len(Labels)):
                sbbgrp1.create_dataset(str(Labels[c]),data=[concentrations[c]])

            sbbgrp1.create_dataset('psi0x',data=[result.psi[0]])
            sbbgrp1.create_dataset('psi1x',data=[result.psi[1]])
            sbbgrp1.create_dataset('psi2x',data=[result.psi[2]])
            sbbgrp1.create_dataset('psidx',data=[result.psi[3]])
            sbbgrp1.create_dataset('psi0z',data=[result.psi[4]])
            sbbgrp1.create_dataset('psi1z',data=[result.psi[5]])
            sbbgrp1.create_dataset('psi2z',data=[result.psi[6]])
            sbbgrp1.create_dataset('psidz',data=[result.psi[7]])
            sbbgrp1.create_dataset('sig0x',data=[result.sig[0]*result.fact2_xy])
            sbbgrp1.create_dataset('sig1x',data=[result.sig[1]*result.fact2_xy])
            sbbgrp1.create_dataset('sig2x',data=[result.sig[2]*result.fact2_xy])
            sbbgrp1.create_dataset('sigdx',data=[result.sig[3]*result.fact2_xy])
            sbbgrp1.create_dataset('sig0z',data=[result.sig[4]*result.fact2_z])
            sbbgrp1.create_dataset('sig1z',data=[result.sig[5]*result.fact2_z])
            sbbgrp1.create_dataset('sig2z',data=[result.sig[6]*result.fact2_z])
            sbbgrp1.create_dataset('sigdz',data=[result.sig[7]*result.fact2_z])
            sbbgrp1.create_dataset('sig0_eqx',data=[result.sig[0]])
            sbbgrp1.create_dataset('sig1_eqx',data=[result.sig[1]])
            sbbgrp1.create_dataset('sig2_eqx',data=[result.sig[2]])
            sbbgrp1.create_dataset('sigd_eqx',data=[result.sig[3]])
            sbbgrp1.create_dataset('sig0_eqz',data=[result.sig[4]])
            sbbgrp1.create_dataset('sig1_eqz',data=[result.sig[5]])
            sbbgrp1.create_dataset('sig2_eqz',data=[result.sig[6]])
            sbbgrp1.create_dataset('sigd_eqz',data=[result.sig[7]])
            sbbgrp1.create_dataset('osmo',data=[osmo])            
            sbbgrp1.create_dataset('ch_bal_EDL',data=[num.sum(result.sig[:4]*result.fact2_xy)+num.sum(result.sig[4:]*result.fact2_z)+sig_diff_bvp])
                
            sbbgrp1.create_dataset('xlist', data=x,compression="gzip")
            sbbgrp1.create_dataset('feld', data=feld,compression="gzip")    
            sbbgrp1.create_dataset('iball', data =ibal,compression="gzip")
            sbbgrp1.create_dataset('psi', data =psi,compression="gzip")

labels_ALL = ['dist','KH_2','kappa','ionic_str','psi0x','psi1x','psi2x','psidx','psi0z','psi1z','psi2z','psidz','sig0x','sig1x','sig2x','sigdx','sig0z','sig1z','sig2z','sigdz','sig0_eqx','sig1_eqx','sig2_eqx','sigd_eqx','sig0_eqz','sig1_eqz','sig2_eqz','sigd_eqz','osmo','ch_bal_EDL','ch_bal_aq', 'sum_resid']
units_ALL = ['(nm)', ' ','(nm**-1)','(mol/L)','(V)','(V)','(V)','(V)','(V)','(V)','(V)','(V)','(C/m**2)','(C/m**2)','(C/m**2)','(C/m**2)','(C/m**2)','(C/m**2)','(C/m**2)','(C/m**2)','(mol(eq))','(mol(eq))','(mol(eq))','(mol(eq))','(mol(eq))','(mol(eq))','(mol(eq))','(mol(eq))','Pa','(C/m**2)','(mol(eq)/L)',' ']
for i in range(len(Labels)):
    labels_ALL.append(Labels[i])
    units_ALL.append('(mol/L)')
base.create_dataset('Labels', data=labels_ALL)
base.create_dataset('Units', data=units_ALL)
base.create_dataset('Distlist', data=D)

base.close()        
print 'time needed (s): ' + str(time.clock() - t0)
############################### THE END #######################################
