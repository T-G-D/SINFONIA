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

from bvp import solve_bvp
import numpy as num

############## other functions ######################################    
def calc_IS(Z, C):
    return 0.5 * num.sum(Z**2*10**C) 
############################### surface complexation model calculations #######
############################### all performed in simplex_point ################
class simplex_point:
    #class that contains all information and functions to perform surface complexation model calculation
    def __init__(self,K,A,T,Zcomp, Zspec,IS,X,cap, surf, dist, Activities = True):
        self.surf = surf
        self.dist = dist
        self.fact = -num.log10(num.exp(1))*1.60217662e-19/(1.38064852e-23 *298.15)
        self.fact2 = 1.60217662e-19*6.02214e23/ (self.surf[1]*self.surf[2])
        self.K = K
        self.T = T
        self.Zcomp = Zcomp
        self.Zspec = Zspec
        self.A = A[:]
        self.Atrans = num.transpose(self.A)
        self.cap = cap
        self.adj = X[1:]
        self.X = X[:]
        self.C = []
        self.G1 = []
        self.G2 = []
        self.Ttmp = []
        self.Etmp = []
        self.psi = []
        self.IS = IS
        self.dpsi_d = 0.0
        self.MB = 0.0
        self.EB = 0.0
        self.worst = 10000
        self.J = 0.0
        self.Activities = Activities
        

    def calc_EB(self, FULL_OUTPUT):
        #calculation of electrostatics
        self.psi = self.X[-2:]/self.fact
        self.sig = self.Ttmp[-2:]
        self.Etmp = num.zeros((2),float)
        self.Etmp[0] = (self.cap[0]*(self.psi[0]-self.psi[1]))/self.fact2
        self.Y[-2] = self.sig[0] - self.Etmp[0]
        if FULL_OUTPUT:
           psi, feld, sig, sig_diff_bvp, osmo = self.solve_PB(FULL_OUTPUT)
        else:
            self.solve_PB(FULL_OUTPUT)
        self.Etmp[1] = (- self.dpsi_d * 78.5 * 8.854187e-12 -(self.cap[0]*(self.psi[0]-self.psi[1])))/self.fact2
        self.Y[-1] = self.sig[1] - self.Etmp[1]
        if FULL_OUTPUT:
            return psi, feld, sig, sig_diff_bvp, osmo
        else:
            return
        
    def calc_IS(self):
        #calculation of ionic strength and Davies activity coefficients
        if self.Activities:
            self.IS = 0.5 * num.sum(self.Zspec**2*10**(self.C-self.G2))
            self.calc_act_coeff2()
            self.IS = 0.5 * num.sum(self.Zspec**2*10**(self.C-self.G2))
            self.calc_act_coeff1()
            self.calc_act_coeff2()
        else:
            self.IS = 0.5 * num.sum(self.Zspec**2*10**(self.C))
        return
    
    def calc_act_coeff1(self):
        #calculation of Davies activity coefficients for Solution Master species
        if self.Activities:
            self.G1 = num.ndarray((len(self.Zcomp)), float)
            for i in range(len(self.G1)):
                if i == 0:
                    self.G1[i] = 0
                elif i > 0:
                    self.G1[i] = -0.5092*self.Zcomp[i]**2 * ((self.IS**0.5 / (1 + self.IS**0.5))-0.3*self.IS)
        else:
            self.G1 = num.zeros((len(self.Zcomp)), float)
        return

    def calc_act_coeff2(self):
        #calculation of Davies activity coefficients for all Solution species
        if self.Activities:
            self.G2 = num.ndarray((len(self.Zspec)), float)
            for i in range(len(self.G2)):
                self.G2[i] = -0.5092*self.Zspec[i]**2 * ((self.IS**0.5 / (1 + self.IS**0.5)) -0.3*self.IS)
        else:
            self.G2 = num.zeros((len(self.Zspec)), float)
        return

    def calc_MB(self, EL):
        #calculation of chemical Mass balance
        self.X[1:] = self.adj[:]
        self.Y = num.zeros((len(self.X)),float)
        self.calc_act_coeff1()
        self.C = self.K + num.dot(self.A, (self.X+self.G1))
        self.calc_act_coeff2()
        self.Ttmp = num.dot(self.Atrans, 10**(self.C-self.G2))
        self.Y = self.Ttmp - self.T
        if EL:
            self.calc_EB(False)
        return
################# scipy PB solver ##############################################
    def solve_PB(self, FULL_OUTPUT = False):
        #function to solve the boundary value problem using Poisson Boltzmann for arbitrary electrolyte
        #some constants
        ew = 78.5 * 8.854187e-12               # permittivity of water (F/m)
        kbt = 1.38064852e-23 *298.15           # kb (J/K) * T in K
        e = 1.60217662e-19                     # electron charge in C
        A =6.02214e23 * 1000/ew                # a prefactor = Avogadro * 1000 /ew

        Q = self.Zspec * e
        C = 10**(self.C-self.G2)
        psi_d = self.psi[1]
        dbl = self.dist *1e-9                  # dbl = considered x-range

        tol_PB = 1e-6                      
        x = num.linspace(0,dbl,100)            # x-axis, discretizatioan of dbl, n steps

        #numerical solution of PB using scipy solve_bvp
        y0 = num.zeros((2, x.size), dtype = float)
        y0[0,0] = psi_d
        y0[1,0] = 0
        y0[0,-1]= psi_d
        y0[1,-1]= 0

        args = [Q,C,A,kbt,y0]

        def fun(x, y, args):
            Q = args[0]
            C = args[1]
            A = args[2]
            kbt = args[3]
            arg1 = num.zeros((x.size))
            for i in range(len(Q)):
                arg1 += Q[i]*C[i]*num.exp(-Q[i]*y[0]/kbt)
            arg1 = -A*arg1
            return num.vstack((y[1] , arg1))

        def bc(ya, yb, args):
            y0 = args[4]
            return num.array([ya[0]-y0[0,0] , yb[0]-y0[0,-1]])

        result = solve_bvp(fun, bc, x, y0, tol = tol_PB, max_nodes = 1000, args = args)
        self.dpsi_d = result.sol(0)[1]
        if FULL_OUTPUT == False:
            return
        elif FULL_OUTPUT == True:
            x_plot = num.linspace(0,dbl,1001)
            psi = result.sol(x_plot)[0]
            feld = result.sol(x_plot)[1] 
            sig = num.zeros((x_plot.size), dtype = float)
            for i in range(len(Q)):
                sig = sig + A*ew*Q[i]*C[i]* num.exp(-Q[i]*psi/kbt)
            sig_diff_bvp = num.trapz(sig, dx =x_plot[1])

            #calculation of osmotic pressure and ion_charge_balance throughout the pore
            ion_balance = num.zeros(len(psi), dtype = float)
            osmo = num.zeros(len(psi), dtype = float)
            for i in range(len(C)):
                osmo = osmo + C[i]*(num.exp(-Q[i]*psi/kbt)-1)   #for two parallel plates
                ion_balance = ion_balance + Q[i]/e*C[i]*num.exp(-Q[i]*psi/kbt)
            osmo = 6.02214e23 * 1000 * kbt * osmo -ew/2 * feld**2
            osmo = osmo[int(len(osmo)/2)]  #osmotic pressure is independent of distance to wall, therefore represented by one selected value close to the middle 
            return psi, feld, ion_balance, sig_diff_bvp, osmo
################################# NR main routine ###########################
def NR_fit(params, Bolz, point, EL):
    
    tr = params
    Bolz = Bolz
    #calculate residuals for the initial guess, Y 
    point.calc_MB(EL)
    point.calc_IS()
    point.calc_MB(EL)

    target_residual = point.T[1:]*tr
    target_residual[target_residual==0.0]= tr

    j = True
    count = 0
  
    while j:
        Y = point.Y[1:]
        #calculate Jacobian according to Westall 1980
        point.J = num.zeros((len(point.adj), len(point.adj)),float)
        for a in range(len(point.adj)):
            for b in range(len(point.adj)):
                for i in range(len(point.K)-1):
                    point.J[a][b] += point.A[i+1][a+1]*point.A[i+1][b+1]*10**point.C[i+1]/10**point.X[b+1]
        if EL:            
            a = len(point.adj)-1
            c = num.log10(num.exp(1))
            point.J[a-1][a-1] -= point.cap[0]/point.fact/point.fact2/10**point.adj[a-1]*c
            point.J[a-1][a] += point.cap[0]/point.fact/point.fact2/10**point.adj[a]*c
            point.J[a][a-1] += point.cap[0]/point.fact/point.fact2/10**point.adj[a-1]*c
            point.J[a][a] -= point.cap[0]/point.fact/point.fact2/10**point.adj[a]*c

            dpsid = point.dpsi_d
            point.adj[a] = point.adj[a] + num.log10(1+1e-4)
            point.psi = point.adj[-2:]/point.fact
            point.solve_PB(False)
            point.adj[a] = point.adj[a] - num.log10(1+1e-4)
            dpsid = (point.dpsi_d-dpsid)/(10**point.adj[a]*1e-4)
            point.J[a][a] += dpsid * 78.5 * 8.854187e-12/point.fact2
        
        Xtmp = 10**point.adj - num.dot(num.linalg.inv(point.J), Y)
           
        # sometimes when initial concentrations are low, Xtmp ranges above and below the target value, not converging.
        # To avoid this, we use the Bolzano approx: X(n+1) =(X(n+1)+X(n))/2 or, in case needed, its general equation:
        # X(n+1) =(X(n+1)+((2**i-1)X(n)))/(2**i)
        point.adj = num.log10((Xtmp+(2**Bolz-1)*10**point.adj)/(2**Bolz))
      
        point.calc_MB(EL)
        point.calc_IS()
        point.calc_MB(EL)

        count +=1
        
        Bool = True
        k = 0
        while Bool and k < len(Y):
            if num.fabs(Y[k]) > target_residual[k]:
                Bool = False
            k += 1

        if Bool == False:
            j = True
        else:
            j = False

    T_mod = point.T[1:]
    T_mod[T_mod==0.0]= 1.
    point.MB = num.max(num.fabs(point.Y[1:]/T_mod))       
    
    print str(count)+' iterations needed; final residual: '+str(point.MB)+'\n'

    return point

################################################################################
