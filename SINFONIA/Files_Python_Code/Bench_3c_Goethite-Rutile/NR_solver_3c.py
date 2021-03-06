from bvp import solve_bvp
import numpy as num

############## other functions ######################################    
def calc_IS(Z, C):
    return 0.5 * num.sum(Z**2*10**C) 
############################### surface complexation model calculations #######
############################### all performed in simplex_point ################
class simplex_point:
    #class used in NR_fit, that contains all information and functions to perform mass balance and electrostaic balance calculations
    def __init__(self,K,A,B,T,Zcomp, Zspec,IS,X,cap_xy,cap_z, surf_xy, surf_z, dist, eps_r, eps_0, Temp, kb, e, Activities = True):
        self.surf_xy = surf_xy
        self.surf_z = surf_z
        self.dist = dist
        self.fact = -num.log10(num.exp(1))*e/(kb*Temp)
        self.fact2_xy = e*6.02214086e23/ (self.surf_xy[1]*self.surf_xy[2])
        self.fact2_z = e*6.02214086e23/ (self.surf_z[1]*self.surf_z[2])
        self.K = K
        self.T = T
        self.Zcomp = Zcomp
        self.Zspec = Zspec
        self.A = A[:]
        self.B = B[:]
        self.Atrans = num.transpose(self.B)
        self.cap_xy = cap_xy
        self.cap_z = cap_z
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
        self.eps_r = eps_r
        self.eps_0 = eps_0
        self.Temp = Temp
        self.kb = kb
        self.e = e
        self.Activities = Activities
        

    def calc_EB(self, FULL_OUTPUT):
        #calculation of electrostatics
        self.psi = self.X[-8:]/self.fact
        self.sig = self.Ttmp[-8:]
        self.Etmp = num.zeros((8),float)
        self.Etmp[0] = (self.cap_xy[0]*(self.psi[0]-self.psi[1]))/self.fact2_xy
        self.Etmp[1] = (self.cap_xy[0]*(self.psi[1]-self.psi[0]) + self.cap_xy[1]*(self.psi[1]-self.psi[2]))/self.fact2_xy
        self.Etmp[2] = (self.cap_xy[1]*(self.psi[2]-self.psi[1]) + self.cap_xy[2]*(self.psi[2]-self.psi[3]))/self.fact2_xy
        self.Etmp[4] = (self.cap_z[0]*(self.psi[4]-self.psi[5]))/self.fact2_z
        self.Etmp[5] = (self.cap_z[0]*(self.psi[5]-self.psi[4]) + self.cap_z[1]*(self.psi[5]-self.psi[6]))/self.fact2_z
        self.Etmp[6] = (self.cap_z[1]*(self.psi[6]-self.psi[5]) + self.cap_z[2]*(self.psi[6]-self.psi[7]))/self.fact2_z
        if FULL_OUTPUT:
           psi, feld, sig, sig_diff_bvp, osmo = self.solve_PB(FULL_OUTPUT)
        else:
            self.solve_PB(FULL_OUTPUT)

        self.Etmp[3] = (self.cap_xy[2]*(self.psi[3]-self.psi[2]) - (self.dpsi_dxy * self.eps_r * self.eps_0))/self.fact2_xy 
        self.Etmp[7] = (self.cap_z[2]*(self.psi[7]-self.psi[6]) - (self.dpsi_dz * self.eps_r * self.eps_0))/self.fact2_z 
        self.Y[-8] = self.sig[0] - self.Etmp[0]
        self.Y[-7] = self.sig[1] - self.Etmp[1] 
        self.Y[-6] = self.sig[2] - self.Etmp[2]
        self.Y[-5] = self.sig[3] - self.Etmp[3]
        self.Y[-4] = self.sig[4] - self.Etmp[4]
        self.Y[-3] = self.sig[5] - self.Etmp[5] 
        self.Y[-2] = self.sig[6] - self.Etmp[6]
        self.Y[-1] = self.sig[7] - self.Etmp[7]
        
        if FULL_OUTPUT:
            return psi, feld, sig, sig_diff_bvp, osmo
        else:
            return
        
    def calc_IS(self):
        #calculation of ionic strength and Davies activity coefficients, valid for IS < 0.5M
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
        #calculation of Davies activity coefficients for components
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
        #calculation of chemical mass balance
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
        ew = self.eps_r * self.eps_0           # permittivity of water (F/m)
        kbt = self.kb *self.Temp               # kb (J/K) * T in K                  
        A =6.02214086e23 * 1000/ew             # a prefactor = Avogadro * 1000 /ew
        Q = self.Zspec * self.e                # for self.e the electron charge in C
        C = 10**(self.C-self.G2)
        psi_dxy = self.psi[-5]
        psi_dz = self.psi[-1]

        dbl = self.dist *1e-9                  # dbl = considered x-range

        tol_PB = 1e-6                      
        x = num.linspace(0,dbl,100)            # x-axis, discretizatioan of dbl, n steps

        #numerical solution of PB using scipy solve_bvp
        y0 = num.zeros((2, x.size), dtype = float)
        y0[0,0] = psi_dxy
        y0[1,0] = 0
        y0[0,-1]= psi_dz
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
        
        self.dpsi_dxy = result.sol(0)[1]
        self.dpsi_dz = -result.sol(dbl)[1]

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
                ion_balance = ion_balance + Q[i]/self.e*C[i]*num.exp(-Q[i]*psi/kbt)
            osmo = 6.02214086e23 * 1000 * kbt * osmo -ew/2 * feld**2
            osmo = osmo[int(len(osmo)/2)]  #osmotic pressure is independent of distance to wall, therefore represented by one selected value close to the middle 
            return psi, feld, ion_balance, sig_diff_bvp, osmo
################################# NR main routine ###########################
def NR_fit(params, Bolz, point, EL, eps_r, eps_0):
    #main routine optimizing the chemical equilibrium system
    tr = params
    Bolz = Bolz
    eps_r = eps_r
    eps_0 = eps_0
    
    #calculate residuals for the initial guess, Y 
    point.calc_MB(EL)
    point.calc_IS()
    point.calc_MB(EL)

    target_residual = point.T[1:]*tr
    target_residual[target_residual==0.0]= tr

    j = True
    count = 0
    conv_fact_xy = point.fact*point.fact2_xy/num.log10(num.exp(1))
    conv_fact_z = point.fact*point.fact2_z/num.log10(num.exp(1))
  
    while j:
        Y = point.Y[1:]
        #calculate Jacobian according to Westall 1980
        point.J = num.zeros((len(point.adj), len(point.adj)),float)
        for a in range(len(point.adj)):
            for b in range(len(point.adj)):
                for i in range(len(point.K)-1):
                    point.J[a][b] += point.B[i+1][a+1]*point.A[i+1][b+1]*10**point.C[i+1]/10**point.X[b+1]
        a = len(point.adj)-1

        if EL:
            #electrostatics for surface TiOH
            point.J[a-7][a-7] -= point.cap_xy[0]/conv_fact_xy/10**point.adj[a-7]
            point.J[a-7][a-6] += point.cap_xy[0]/conv_fact_xy/10**point.adj[a-6]
            point.J[a-6][a-7] += point.cap_xy[0]/conv_fact_xy/10**point.adj[a-7]            
            point.J[a-6][a-6] -= (point.cap_xy[0]+point.cap_xy[1])/conv_fact_xy/10**point.adj[a-6]
            point.J[a-6][a-5] += point.cap_xy[1]/conv_fact_xy/10**point.adj[a-5]
            point.J[a-5][a-6] += point.cap_xy[1]/conv_fact_xy/10**point.adj[a-6]
            point.J[a-5][a-5] -= (point.cap_xy[1]+point.cap_xy[2])/conv_fact_xy/10**point.adj[a-5]
            point.J[a-5][a-4] += point.cap_xy[2]/conv_fact_xy/10**point.adj[a-4]
            point.J[a-4][a-5] += point.cap_xy[2]/conv_fact_xy/10**point.adj[a-5]
            point.J[a-4][a-4] -= point.cap_xy[2]/conv_fact_xy/10**point.adj[a-4]

            #electrostatics for surface FeOH
            point.J[a-3][a-3] -= point.cap_z[0]/conv_fact_z/10**point.adj[a-3]
            point.J[a-3][a-2] += point.cap_z[0]/conv_fact_z/10**point.adj[a-2]
            point.J[a-2][a-3] += point.cap_z[0]/conv_fact_z/10**point.adj[a-3]            
            point.J[a-2][a-2] -= (point.cap_z[0]+point.cap_z[1])/conv_fact_z/10**point.adj[a-2]
            point.J[a-2][a-1] += point.cap_z[1]/conv_fact_z/10**point.adj[a-1]
            point.J[a-1][a-2] += point.cap_z[1]/conv_fact_z/10**point.adj[a-2]
            point.J[a-1][a-1] -= (point.cap_z[1]+point.cap_z[2])/conv_fact_z/10**point.adj[a-1]
            point.J[a-1][a] += point.cap_z[2]/conv_fact_z/10**point.adj[a]
            point.J[a][a-1] += point.cap_z[2]/conv_fact_z/10**point.adj[a-1]
            point.J[a][a] -= point.cap_z[2]/conv_fact_z/10**point.adj[a]
                
            #modificatin of Jacobian to accommodate charge regulation boundary conditions
            dpsid_xy1 = point.dpsi_dxy
            dpsid_z1 = point.dpsi_dz        
                
            point.adj[a-4] = point.adj[a-4] + num.log10(1+1e-4)
            point.psi = point.adj[-8:]/point.fact
            point.solve_PB(False)
            point.adj[a-4] = point.adj[a-4] - num.log10(1+1e-4)
            h1 = 10**point.adj[a-4]*1e-4
            slope_1 = (point.dpsi_dxy-dpsid_xy1)/h1
            slope_2 = (point.dpsi_dz-dpsid_z1)/h1
            point.J[a-4][a-4] += slope_1 * eps_r * eps_0/point.fact2_xy           
            point.J[a][a-4] += slope_2 * eps_r * eps_0/point.fact2_z

            point.adj[a] = point.adj[a] + num.log10(1+1e-4)
            point.psi = point.adj[-8:]/point.fact
            point.solve_PB(False)
            point.adj[a] = point.adj[a] - num.log10(1+1e-4)
            h2 = 10**point.adj[a]*1e-4
            slope_3 = (point.dpsi_dxy-dpsid_xy1)/h2
            slope_4 = (point.dpsi_dz-dpsid_z1)/h2
            point.J[a-4][a] += slope_3 * eps_r * eps_0/point.fact2_xy           
            point.J[a][a] += slope_4 * eps_r * eps_0/point.fact2_z
        
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
