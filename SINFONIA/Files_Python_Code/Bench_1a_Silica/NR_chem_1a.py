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

from NR_solver_1a import calc_IS, simplex_point, NR_fit
import numpy as num
from matplotlib import pyplot
import time

t0 = time.clock()

################################################################
########### USER INPUT #########################################
filename = 'Benchmark_1a_'  #filename output in .dat extension
Individual_Files = False    #produce individual output file for each pH and distance
Overview_File = True        #produce general output file with overview of results

pH_range = num.arange(3.0,9.0,1.0) #define the pH range and increment step

T1 = 0.01       #general concentration (M) of electrolyte 1 (e.g. Na+)
T2 = 0.01       #general concentration (M) of electrolyte 2 (e.g. Cl-)

#pH corrections, data extracted from equivalent simulations of the chemical system in PHREEQC.
pH_corr_T1 = [0., 0., 0., 0., 3.400e-07, 1.834e-06]                 #representing addition of NaOH (M)
pH_corr_T2 = [1.114e-03, 1.109e-04, 1.104e-05, 9.500e-07, 0., 0.]   #representing addition of HCl (M)

#List of interparticle distances (in nm):
distances = num.array([28., 25., 23., 20., 18., 15., 13., 10., 7., 5., 4.5, 4., 3.5, 3., 2.5, 2., 1.75, 1.5, 1.25, 1.1, 1., 0.9, 0.8, 0.7, 0.6, 0.55, 0.5, 0.45, 0.4, 0.3, 0.2, 0.1], float) 

# Activate calculations with conversion of concentrations (used in Kc definitions) into activities (required for Kint) 
Activities = True                   # if True = calculations are performed in activities

# Choose the geometry for the force-distance calculations:
cylinder_cylinder = False       # approaching crossed cylinders (True)
sphere_sphere = False           # approaching spheres (True)
sphere_plane = True             # sphere-plane, equivalent to approaching parallel cylinders (True)

H = 1.2e-20                     # Hamaker constant (J), if equal 0 then no VdW taken into account, otherwise in the order of 1e-20
r1 = 1000000000                 # 1st particle radius in nm
r2 = 10000000000                # 2nd particle radius in nm (for cylinder-cylinder or sphere-sphere geometries)

################################################################
########### MODEL DEFINITIONS ##################################
tr = 1.0e-13                    # target residual, implemented such that error in the mass balance of each component must become < tr
Bolz = 1                        # Bolzano's approach used in the Newton-Raphson calculation (>=1)

Labels = ['H+','Na+','Cl-','OH-','>SOH','>SO-']  #Defining labels for the species involved in the chemical system 

Zspec = num.array([1, 1, -1, -1, 0, 0], float)   # Zspec: Vector with product species charges 
#note: surface species have no activity and therefore no charge is specified here

K = num.array([0.0, 0.0, 0.0, -14, 0.0, -6.81], float)  #K: Vector with log10_k values for reactions considered in LMA calculations
#K: K=(H+, Na+, Cl-, OH-, >SOH, >SO-) 

cap = num.array([1e4],float)  #cap: array with capacitances. High values such as 1e4 convert the present script from a BSM to a DLM.

surf = num.array([4.6, 1.0, 1.0],float)  #surf: surface definition in sites/nm**2, m**2/g, g/L.

#A: Matrix, defining stoichiometries of reactions, presented in the following order:
#H+ = 0 + H+
#Na+ = 0 + Na+
#Cl- = 0 + Cl-
#OH- = -14 -H+
#>SOH = >SOH
#>SO- = - H+ + >SOH

A = num.array([\
 [ 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],\
 [ 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],\
 [ 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],\
 [-1.0, 0.0, 0.0, 0.0, 0.0, 0.0],\
 # surface species >SO
 [ 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],\
 [-1.0, 0.0, 0.0, 1.0, -1.0, 0.0],\
 ], float)

################################################################
########### EXECUTION OF MODEL CALCULATIONS ####################

#prepare arrays to collect full output data
psi_stern_list = num.ndarray((0,2),float)
sig_stern_list = num.ndarray((0,2),float)
sig_sterneq_list = num.ndarray((0,2),float)
concentrationlist = num.ndarray((0,len(Labels)),float)
osmolist = []
osmolist_TmpH = []
sig_diff_list = []
cballist = []
MBlist = []
ISlist = []
kappa_distances = []
distlist = []
Forcelist = []
int_osmo = []

# some constants
teta = (298.15 - 273.16)/100
eps_r = 88.15-(41.4+(13.1-4.6*teta)*teta)*teta
eps_0 = 8.85419e-12
e = 1.60219e-19
Far = 9.64853e4
kb = 1.38062e-23
Temp = 298.15

#loop through pH values
for pH in range(len(pH_range)):
    
    # T: Vector with the total concentrations of components
    # T: T = (H, Na, Cl, SOH); total H is defined as activity because pH will be fixed
    T = num.array([10**(-pH_range[pH]), T1+pH_corr_T1[pH], T2+pH_corr_T2[pH], 0.0], float)
    Zcomp = num.array([1, 1, -1, 0], float) #charges of components for activity calculations
    #calculate total concentration of surface sites from surf
    T[-1] = surf[0]/6.02214086e5 * surf[1] * surf[2]

    # X: Starting values for log10 master-species concentrations, H+ is always first with known activity
    X = num.log10(T)
    #append 2 zeros to T and X, starting values for surface charges and potentials for a BSM
    T = num.append(T, [0.,0.])
    X = num.append(X, [0.,0.])
    Zcomp = num.append(Zcomp, [0.,0.])

    # IS ionic strength
    IS = calc_IS(Zcomp, X) # initial guess for the ionic strength

    #loop through distances 
    for dist in distances:
        xlist = num.ndarray((0,1001),float)
        psilist = num.ndarray((0,1001),float)
        feldlist = num.ndarray((0,1001),float)
        iballist = num.ndarray((0,1001),float)
        
        #create an instance of simplex_point, the class holding all necessary values and functions to calculate the chemical equilibrium under charge regulation consitions
        result = simplex_point(K,A,T,Zcomp,Zspec,IS,X,cap,surf, dist, Activities = Activities)
        result.calc_MB(False) #starting point for simplex

        result = NR_fit(tr, Bolz, result, True) #start Newton-Raphson optimization routine, solving the chemical equilibrium system
        #param 1: target residual

        #calculating equivalent kappa distances and KH/2 conversions
        kappa_Zhmud = 1e-9*((2*e*Far*1000*result.IS)/(eps_0*eps_r*kb*Temp))**0.5
        kappa_dist = (kappa_Zhmud*dist)/2
        kappa_distances = num.append(kappa_distances,kappa_dist)
        
        #grab all output from the final result
        psi, feld, ibal, sig_diff_bvp, osmo = result.calc_EB(True)

        # Calculating interacting forces
        distlist.append(-dist)
        osmolist_TmpH.append(osmo)
        
        if cylinder_cylinder:
            Fvdw = -H*num.sqrt(r1*r2*1e-18)/(6*(dist*1e-9)**2)                 # vdW force as defined in Table 13.1 from Israelachvili 2011 
            if dist == distances[0]:
                Wosmo = osmolist_TmpH[0]
                force = 2*num.pi*num.sqrt(r1*r2*1e-18)*Wosmo*dist + Fvdw*1e9   # Force in nN
            else:
                Wosmo = num.trapz(osmolist_TmpH, distlist)
                force = 2*num.pi*num.sqrt(r1*r2*1e-18)*Wosmo + Fvdw*1e9

        if sphere_sphere:
            Fvdw = (-H/(6*(dist*1e-9)**2))*(r1*r2*1e-9/(r1+r2))                 # vdW force as defined in Table 13.1 from Israelachvili 2011 
            if dist == distances[0]:
                Wosmo = osmolist_TmpH[0]
                force = (2*num.pi*r1*r2*1e-9/(r1+r2))*Wosmo + Fvdw*1e9          # Force in nN
            else:
                Wosmo = num.trapz(osmolist_TmpH, distlist)
                force = (2*num.pi*r1*r2*1e-9/(r1+r2))*Wosmo + Fvdw*1e9
             
        if sphere_plane:
            Fvdw = -H*r1/(6*(dist*1e-9)**2)                                     # vdW force as defined in Table 13.1 from Israelachvili 2011 
            if dist == distances[0]:
                Wosmo = osmolist_TmpH[0]
                force = 2*num.pi*r1*1e-9*Wosmo + Fvdw                           # from Trefalt et al. 2016, force in nN
            else:
                Wosmo = num.trapz(osmolist_TmpH, distlist)
                force = 2*num.pi*r1*1e-9*Wosmo + Fvdw

        Forcelist = num.append(Forcelist, force)
        int_osmo.append(Wosmo)
            
        #now we produce some output
        activities = 10**(result.C)
        concentrations = 10**(result.C-result.G2)
        charge = num.sum(result.Zspec*concentrations)
        
        #first print with summary of each condition to the shell
        print "\nOUTPUT: \n"
        print "original totals were: " + str( result.T[1:-2] )
        print "Resulting totals are: "+ str( result.Ttmp[1:-2])+"\n"
        print "mass balance residuals are: "+ str(num.fabs(result.Y[1:-2]))
        print "electrostatic balance residuals are: " + str(num.fabs(result.Y[-2:]))
        
        print "electrical balance in bulk solution (mol(eq)): " + str(charge)
        print "electric balance in EDL (C/m**2): " +str(2*(num.sum(result.sig*result.fact2))+sig_diff_bvp) +'\n'
        
        print 'osmotic pressure (Pa): ' + str(round(osmo,2))
        print 'psi_0 (V): ' + str(round(result.psi[0],5))
        print 'psi_b (V): ' + str(round(result.psi[1],5))
        print 'sig_0 (mol, eq): ' + str(result.sig[0])
        print 'sig_b (mol, eq): ' + str(result.sig[1])
        print 'sig_0 (C/m**2): ' + str(result.sig[0]*result.fact2)
        print 'sig_b (C/m**2): ' + str(result.sig[1]*result.fact2)
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

        #collect results
        xlist    = num.append(xlist,       [x], axis = 0)
        psilist  = num.append(psilist,   [psi], axis = 0)
        feldlist = num.append(feldlist, [feld], axis = 0)
        iballist = num.append(iballist, [ibal], axis = 0)
        psi_stern_list = num.append(psi_stern_list, [result.psi], axis = 0)
        sig_stern_list = num.append(sig_stern_list, [result.sig*result.fact2], axis = 0)
        sig_sterneq_list = num.append(sig_sterneq_list, [result.sig], axis = 0)
        osmolist.append(osmo)
        sig_diff_list.append(2*(num.sum(result.sig*result.fact2))+sig_diff_bvp)
        cballist.append(charge)
        concentrationlist = num.append(concentrationlist, [concentrations], axis = 0)
        MBlist.append(result.MB)
        ISlist.append(result.IS)

        if Individual_Files:
            #second produce an outputfile with all information for one specific distance/pH
            filenamex = filename+'_'+'_pH_'+str(round(pH_range[pH],1))+'_' + str(round(dist,1))+'nm.dat'
            g = file((filenamex), 'w')
            g.write('Kappa-DISTANCE = '+str(round(kappa_dist,1))+'\n')
            g.write('pH = '+str(round(pH_range[pH],1))+'\n\n')
            g.write("fitting residual = "+str(result.MB)+" \n\n")
            g.write("original totals were: " + str( result.T[1:-2] ) +'\n')
            g.write("Resulting totals are: "+ str( result.Ttmp[1:-2])+"\n" +'\n')
            g.write("mass balance residuals are: "+ str(num.fabs(num.log10(result.T[1:-2])-num.log10(result.Ttmp[1:-2]))) + '\n')
            g.write("electrostatic balance residuals are: " + str(result.Etmp) +'\n')
            g.write("ionic strength is (mol/L): " + str( result.IS ) +'\n')
            g.write("electrical balance in bulk solution (mol(eq)): " + str(charge) +'\n')
            g.write("electric balance in EDL (C/m**2): " +str(2*(num.sum(result.sig*result.fact2))+sig_diff_bvp) +'\n' +'\n')
            g.write('osmotic pressure (Pa): ' + str(osmo) +'\n')
            g.write('psi_0 (V): ' + str(result.psi[0]) +'\n')
            g.write('psi_b (V): ' + str(result.psi[1]) +'\n')
            g.write('sig_0 (mol, eq): ' + str(result.sig[0]) +'\n')
            g.write('sig_b (mol, eq): ' + str(result.sig[1]) +'\n')
            g.write('sig_0 (C/m**2): ' + str(result.sig[0]*result.fact2) +'\n')
            g.write('sig_b (C/m**2): ' + str(result.sig[1]*result.fact2) +'\n')
            g.write("species, concentrations, activities,  log10_conc.,  log10_act,  and log10_act_coeff: \n")
            for i in range(len(Labels)):
                line = "%10s   %.3E   %.3E   %.3E   %.3E  %.3E\n" % (Labels[i], concentrations[i], activities[i],result.C[i],(result.C[i]-result.G2[i]), result.G2[i])
                g.write(line)

            g.write("\n\n distance (nm), potential (V), electric field (V/m), ion_charge_balance (mol(eq)/L) \n")
            for i in range(len(xlist[0])):
                line = "%.6E  %.6E  %.6E  %.6E\n" % (xlist[0,i], psilist[0,i], feldlist[0,i], iballist[0,i])
                g.write(line)
            g.close()
        
    #we have to clean the temporal lists for integrating the osmotic pressure over distances for every pH value
    distlist = []
    osmolist_TmpH = []
    
#write overarching results file for all distances and pHs
if Overview_File:
    if cylinder_cylinder:
        filenamex = filename+str('Force_cc')+'_general_output.dat'
    if sphere_sphere:
        filenamex = filename+str('Force_ss')+'_general_output.dat'
    if sphere_plane:
        filenamex = filename+str('Force_sp')+'_general_output.dat'
    g = file((filenamex), 'w')
    g.write(filename+'\n\n')      

    line = "%12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s"
    labellist = ['distance','kappa-distance','pH','psi0','psib','sig0','sigb','sig0_eq','sigb_eq','osmotic_pr','ch._bal_EDL','ch._bal_aq','ionic_str.', 'Integral_Osmo', 'Force', 'sum_resid.']
    for i in range(len(Labels)):
        line = line +"%12s  "
        labellist.append(Labels[i])
    line = line % tuple(labellist)
    line = line +"\n"
    g.write(line)
    line = "%12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s"
    labellist = ['(nm)', ' ',' ','(V)','(V)','(C/m**2)','(C/m**2)','(mol(eq))','(mol(eq))','Pa','(C/m**2)','(mol(eq)/L)','(mol/L)','Pa.nm', 'nN',' ']
    for i in range(len(Labels)):
        line = line +"%12s  "
        labellist.append('(mol/L)')
    line = line % tuple(labellist)
    line = line +"\n\n"
    g.write(line)
         
    for j in range(len(pH_range)):
        for i in range(len(distances)):
            k = j*len(distances)+i
            line = "%.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  "
            labellist = [distances[i],kappa_distances[i],pH_range[j], psi_stern_list[k][0],psi_stern_list[k][1],\
                           sig_stern_list[k][0],sig_stern_list[k][1],sig_sterneq_list[k][0],sig_sterneq_list[k][1],\
                           osmolist[k], sig_diff_list[k], cballist[k], ISlist[k], int_osmo[k], Forcelist[k], MBlist[k]]
            for a in range(len(Labels)):
                line = line +"%.6E  "
                labellist.append(concentrationlist[k][a])
            line = line % tuple(labellist)
            line = line +"\n" 
            g.write(line)
    g.close()
    
print 'time needed (s): ' + str(time.clock() - t0)
############################### THE END #######################################
