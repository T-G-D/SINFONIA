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

t0 = time.clock()
################################################################
########### USER INPUT #########################################
filename = 'Benchmark_3c_' #filename from the library, must be equal to filename defined in "NR_chem_3c.py"

base = h5py.File(filename+'dictionary.hdf5', 'r')
print base.keys()
pH_range = num.array(base.get('pH_range'))
Initial_IS = num.array(base.get('Initial_IS'))
Distlist = num.array(base.get('Distlist'))
Labels = num.array(base[u'Labels'])
print Labels

find = num.char.find(Labels,'sum_resid')
for i in range(len(Labels)):
    if find[i] == 0:
        Labels_sp = Labels[i+1:]
        break

for s in range(len(Initial_IS)):
    distlist = []
    kappa_dist_list = []
    psi0x_list = []
    psi1x_list = []
    psi2x_list = []
    psidx_list = []
    psi0z_list = []
    psi1z_list = []
    psi2z_list = []
    psidz_list = []
    sig0x_list = []
    sig1x_list = []
    sig2x_list = []
    sigdx_list = []
    sig0z_list = []
    sig1z_list = []
    sig2z_list = []
    sigdz_list = []
    sig0_eqx_list = []
    sig1_eqx_list = []
    sig2_eqx_list = []
    sigd_eqx_list = []
    sig0_eqz_list = []
    sig1_eqz_list = []
    sig2_eqz_list = []
    sigd_eqz_list = []
    osmolist = []
    sig_diff_list = []
    cballist = []
    ISlist = []
    MBlist = []
    concentrations = num.ndarray((0,len(Labels_sp)),float)

    for p in range(len(pH_range)):
        for d in range(len(Distlist)):
            dist = num.array(base.get('IS_'+str(Initial_IS[s])+'/pH_'+str(pH_range[p])+'/distance'+str(Distlist[d])+'/dist'))
            kappa_distances = num.array(base.get('IS_'+str(Initial_IS[s])+'/pH_'+str(pH_range[p])+'/distance'+str(Distlist[d])+'/KH_2'))
            psi0x = num.array(base.get('IS_'+str(Initial_IS[s])+'/pH_'+str(pH_range[p])+'/distance'+str(Distlist[d])+'/psi0x'))
            psi1x = num.array(base.get('IS_'+str(Initial_IS[s])+'/pH_'+str(pH_range[p])+'/distance'+str(Distlist[d])+'/psi1x'))
            psi2x = num.array(base.get('IS_'+str(Initial_IS[s])+'/pH_'+str(pH_range[p])+'/distance'+str(Distlist[d])+'/psi2x'))
            psidx = num.array(base.get('IS_'+str(Initial_IS[s])+'/pH_'+str(pH_range[p])+'/distance'+str(Distlist[d])+'/psidx'))
            psi0z = num.array(base.get('IS_'+str(Initial_IS[s])+'/pH_'+str(pH_range[p])+'/distance'+str(Distlist[d])+'/psi0z'))
            psi1z = num.array(base.get('IS_'+str(Initial_IS[s])+'/pH_'+str(pH_range[p])+'/distance'+str(Distlist[d])+'/psi1z'))
            psi2z = num.array(base.get('IS_'+str(Initial_IS[s])+'/pH_'+str(pH_range[p])+'/distance'+str(Distlist[d])+'/psi2z'))
            psidz = num.array(base.get('IS_'+str(Initial_IS[s])+'/pH_'+str(pH_range[p])+'/distance'+str(Distlist[d])+'/psidz'))
            sig0x = num.array(base.get('IS_'+str(Initial_IS[s])+'/pH_'+str(pH_range[p])+'/distance'+str(Distlist[d])+'/sig0x'))
            sig1x = num.array(base.get('IS_'+str(Initial_IS[s])+'/pH_'+str(pH_range[p])+'/distance'+str(Distlist[d])+'/sig1x'))
            sig2x = num.array(base.get('IS_'+str(Initial_IS[s])+'/pH_'+str(pH_range[p])+'/distance'+str(Distlist[d])+'/sig2x'))
            sigdx = num.array(base.get('IS_'+str(Initial_IS[s])+'/pH_'+str(pH_range[p])+'/distance'+str(Distlist[d])+'/sigdx'))
            sig0z = num.array(base.get('IS_'+str(Initial_IS[s])+'/pH_'+str(pH_range[p])+'/distance'+str(Distlist[d])+'/sig0z'))
            sig1z = num.array(base.get('IS_'+str(Initial_IS[s])+'/pH_'+str(pH_range[p])+'/distance'+str(Distlist[d])+'/sig1z'))
            sig2z = num.array(base.get('IS_'+str(Initial_IS[s])+'/pH_'+str(pH_range[p])+'/distance'+str(Distlist[d])+'/sig2z'))
            sigdz = num.array(base.get('IS_'+str(Initial_IS[s])+'/pH_'+str(pH_range[p])+'/distance'+str(Distlist[d])+'/sigdz'))
            sig0_eqx = num.array(base.get('IS_'+str(Initial_IS[s])+'/pH_'+str(pH_range[p])+'/distance'+str(Distlist[d])+'/sig0_eqx'))
            sig1_eqx = num.array(base.get('IS_'+str(Initial_IS[s])+'/pH_'+str(pH_range[p])+'/distance'+str(Distlist[d])+'/sig1_eqx'))
            sig2_eqx = num.array(base.get('IS_'+str(Initial_IS[s])+'/pH_'+str(pH_range[p])+'/distance'+str(Distlist[d])+'/sig2_eqx'))
            sigd_eqx = num.array(base.get('IS_'+str(Initial_IS[s])+'/pH_'+str(pH_range[p])+'/distance'+str(Distlist[d])+'/sigd_eqx'))
            sig0_eqz = num.array(base.get('IS_'+str(Initial_IS[s])+'/pH_'+str(pH_range[p])+'/distance'+str(Distlist[d])+'/sig0_eqz'))
            sig1_eqz = num.array(base.get('IS_'+str(Initial_IS[s])+'/pH_'+str(pH_range[p])+'/distance'+str(Distlist[d])+'/sig1_eqz'))
            sig2_eqz = num.array(base.get('IS_'+str(Initial_IS[s])+'/pH_'+str(pH_range[p])+'/distance'+str(Distlist[d])+'/sig2_eqz'))
            sigd_eqz = num.array(base.get('IS_'+str(Initial_IS[s])+'/pH_'+str(pH_range[p])+'/distance'+str(Distlist[d])+'/sigd_eqz'))
            osmo = num.array(base.get('IS_'+str(Initial_IS[s])+'/pH_'+str(pH_range[p])+'/distance'+str(Distlist[d])+'/osmo'))
            ch_bal_EDL = num.array(base.get('IS_'+str(Initial_IS[s])+'/pH_'+str(pH_range[p])+'/distance'+str(Distlist[d])+'/ch_bal_EDL'))
            ch_bal_aq = num.array(base.get('IS_'+str(Initial_IS[s])+'/pH_'+str(pH_range[p])+'/distance'+str(Distlist[d])+'/ch_bal_aq'))
            IS_calc = num.array(base.get('IS_'+str(Initial_IS[s])+'/pH_'+str(pH_range[p])+'/distance'+str(Distlist[d])+'/ionic_str'))
            sum_resid = num.array(base.get('IS_'+str(Initial_IS[s])+'/pH_'+str(pH_range[p])+'/distance'+str(Distlist[d])+'/sum_resid'))
            
            for c in range(len(Labels_sp)):
                conc = str(Labels_sp[c])
                sp = num.array(base.get('IS_'+str(Initial_IS[s])+'/pH_'+str(pH_range[p])+'/distance'+str(Distlist[d])+'/'+str(conc)))
                concentrations = num.append(concentrations,[sp])
                        
            distlist = num.append(distlist, dist)
            kappa_dist_list = num.append(kappa_dist_list, kappa_distances)
            psi0x_list = num.append(psi0x_list, psi0x)
            psi1x_list = num.append(psi1x_list, psi1x)
            psi2x_list = num.append(psi2x_list, psi2x)
            psidx_list = num.append(psidx_list, psidx)
            psi0z_list = num.append(psi0z_list, psi0z)
            psi1z_list = num.append(psi1z_list, psi1z)
            psi2z_list = num.append(psi2z_list, psi2z)
            psidz_list = num.append(psidz_list, psidz)
            sig0x_list = num.append(sig0x_list, sig0x)
            sig1x_list = num.append(sig1x_list, sig1x)
            sig2x_list = num.append(sig2x_list, sig2x)
            sigdx_list = num.append(sigdx_list, sigdx)
            sig0z_list = num.append(sig0z_list, sig0z)
            sig1z_list = num.append(sig1z_list, sig1z)
            sig2z_list = num.append(sig2z_list, sig2z)
            sigdz_list = num.append(sigdz_list, sigdz)
            sig0_eqx_list = num.append(sig0_eqx_list, sig0_eqx)
            sig1_eqx_list = num.append(sig1_eqx_list, sig1_eqx)
            sig2_eqx_list = num.append(sig2_eqx_list, sig2_eqx)
            sigd_eqx_list = num.append(sigd_eqx_list, sigd_eqx)
            sig0_eqz_list = num.append(sig0_eqz_list, sig0_eqz)
            sig1_eqz_list = num.append(sig1_eqz_list, sig1_eqz)
            sig2_eqz_list = num.append(sig2_eqz_list, sig2_eqz)
            sigd_eqz_list = num.append(sigd_eqz_list, sigd_eqz)
            osmolist = num.append(osmolist, osmo)
            sig_diff_list = num.append(sig_diff_list, ch_bal_EDL)
            cballist = num.append(cballist, ch_bal_aq)
            ISlist = num.append(ISlist, IS_calc)
            MBlist = num.append(MBlist, sum_resid)
    
    #write overarching results file for all distances and pHs
    filenamex = filename+str(Initial_IS[s])+'general_output.dat'
    g = file((filenamex), 'w')
    g.write(filename+'\n\n')
            
    line = "%12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s"
    labellist = ['distance','KH/2','pH','psi0x','psi1x','psi2x','psidx','psi0z','psi1z','psi2z','psidz','sig0x','sig1x','sig2x','sigdx','sig0z','sig1z','sig2z','sigdz','sig0_eqx','sig1_eqx','sig2_eqx','sigd_eqx','sig0_eqz','sig1_eqz','sig2_eqz','sigd_eqz','osmotic_pr','ch._bal_EDL','ch._bal_aq','ionic_str.', 'sum_resid.']
    for i in range(len(Labels_sp)):
        line = line +"%12s  "
        labellist.append(Labels_sp[i])
    line = line % tuple(labellist)
    line = line +"\n"
    g.write(line)
    line = "%12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s"
    labellist = ['(nm)', ' ',' ','(V)','(V)','(V)','(V)','(V)','(V)','(V)','(V)','(C/m**2)','(C/m**2)','(C/m**2)','(C/m**2)','(C/m**2)','(C/m**2)','(C/m**2)','(C/m**2)','(mol(eq))','(mol(eq))','(mol(eq))','(mol(eq))','(mol(eq))','(mol(eq))','(mol(eq))','(mol(eq))','Pa','(C/m**2)','(mol(eq)/L)','(mol/L)',' ']
    for i in range(len(Labels_sp)):
        line = line +"%12s  "
        labellist.append('(mol/L)')
    line = line % tuple(labellist)
    line = line +"\n\n"
    g.write(line)

    for j in range(len(pH_range)):
        for i in range(len(Distlist)):
            k = j*len(Distlist)+i
            line = "%.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  %.6E  "                
            labellist = [distlist[k],kappa_dist_list[k],pH_range[j], psi0x_list[k],psi1x_list[k],psi2x_list[k],psidx_list[k],\
                             psi0z_list[k],psi1z_list[k],psi2z_list[k],psidz_list[k],sig0x_list[k],sig1x_list[k],sig2x_list[k],sigdx_list[k],\
                             sig0z_list[k],sig1z_list[k],sig2z_list[k],sigdz_list[k],sig0_eqx_list[k],sig1_eqx_list[k],sig2_eqx_list[k],sigd_eqx_list[k],\
                             sig0_eqz_list[k],sig1_eqz_list[k],sig2_eqz_list[k],sigd_eqz_list[k],osmolist[k], sig_diff_list[k], cballist[k], ISlist[k], MBlist[k]]
            for a in range(len(Labels_sp)):
                v = a + (j*len(Distlist) + i)*len(Labels_sp)
                line = line +"%.6E  "
                labellist.append(concentrations[v])                    
            line = line % tuple(labellist)
            line = line +"\n" 
            g.write(line)
    g.close()

base.close()        
print 'time needed (s): ' + str(time.clock() - t0)
############################### THE END #######################################
