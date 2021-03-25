function [T] = Benchmark_Rutile_Goethite (distance, pH, n_val)
% Model exported on Jun 15 2020, 09:11 by COMSOL 5.4.0.346.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');
model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 1);

model.component('comp1').mesh.create('mesh1');

model.component('comp1').physics.create('es', 'Electrostatics', 'geom1');
model.component('comp1').physics.create('ge', 'GlobalEquations', 'geom1');

model.study.create('std1');
model.study('std1').setGenConv(true);
model.study('std1').create('stat', 'Stationary');
model.study('std1').feature('stat').activate('es', true);
model.study('std1').feature('stat').activate('ge', true);
%% Parameter
model.param.set('L0', '0 [nm]');
model.param.set('Lf', '100 [nm]');
model.param.set('NA', 'N_A_const');
model.param.set('kb', 'k_B_const');
model.param.set('e', 'e_const');
model.param.set('F', 'F_const');
model.param.set('epsilon0', 'epsilon0_const');
model.param.set('R', 'R_const');
model.param.set('mi_0', '1 [mol/L]');
model.param.set('T', '298.15 [K]');
model.param.set('e_kbt', 'e/kbt');
model.param.set('kbt', 'kb*T');
model.param.set('teta', '(T-273.15)/100  [K]');
model.param.set('epsilon', '88.15-(41.4+(13.1-4.6*teta)*teta)*teta');
model.param.set('ew', 'epsilon0*epsilon');
model.param.set('pH', pH);
model.param.set('log_a_H', '-pH');
model.param.set('log_K_OH', '-14');
model.param.set('log_K_TiOH2', '5.8');
model.param.set('log_K_FeOH2','9.5');
model.param.set('z_H', '1');
model.param.set('z_OH', '-1');
model.param.set('z_C', '1');
model.param.set('z_A', '-1');
model.param.set('C_noextra', '0.001 [mol/L]');
model.param.set('A_noextra', '0.001 [mol/L]');
model.param.set('extra_C', '0 [mol/L]');
model.param.set('extra_A', '0 [mol/L]');
model.param.set('T_A', 'A_noextra + extra_A');
model.param.set('T_C', 'C_noextra + extra_C');

model.param.set('sites_a', '12.2 [1/nm^2]');
model.param.set('surf_a', '1 [m^2/g]');
model.param.set('mass_a', '1 [g/L]');
model.param.set('sites_b', '6.15 [1/nm^2]');
model.param.set('surf_b', '1 [m^2/g]');
model.param.set('mass_b', '1 [g/L]');
model.param.set('sa_F_a', '(surf_a*mass_a)/F');
model.param.set('sa_F_b', '(surf_b*mass_b)/F');
model.param.set('Ai', '0.5092');
model.param.set('F_RT', 'F/(R*T)');
model.param.set('C1a', '1.33 [F/m^2]');
model.param.set('C2a', '1e-4 [F/m^2]');
model.param.set('C3a', '1e-4 [F/m^2]');
model.param.set('C1b', '1.10 [F/m^2]');
model.param.set('C2b', '1e-4 [F/m^2]');
model.param.set('C3b', '1e-4 [F/m^2]');
model.param.set('T_surf_a', '(sites_a*surf_a*mass_a)/NA');
model.param.set('T_surf_b', '(sites_b*surf_b*mass_b)/NA');

%% geo

model.component('comp1').geom('geom1').create('i1', 'Interval');
model.component('comp1').geom('geom1').feature('i1').setIndex('coord', 'Lf', 1);
model.component('comp1').geom('geom1').feature('i1').setIndex('coord', 'L0', 0);
model.component('comp1').geom('geom1').run('fin');

%% Coupling
model.component('comp1').cpl.create('intop1', 'Integration');
model.component('comp1').cpl('intop1').set('axisym', true);
model.component('comp1').cpl.create('intop2', 'Integration');
model.component('comp1').cpl('intop2').set('axisym', true);
model.component('comp1').cpl('intop1').selection.geom('geom1', 0);
model.component('comp1').cpl('intop1').selection.set([1]);
model.component('comp1').cpl('intop2').selection.geom('geom1', 0);
model.component('comp1').cpl('intop2').selection.set([2]);

%% Variables
model.component('comp1').variable.create('var1');
model.component('comp1').variable('var1').set('log_coeff_H', '-Ai*(z_H^2)*((sqrt(IS)/(1+sqrt(IS)))-0.3*IS )');
model.component('comp1').variable('var1').set('log_coeff_OH', '-Ai*(z_OH^2)*((sqrt(IS)/(1+sqrt(IS)))-0.3*IS )');
model.component('comp1').variable('var1').set('log_coeff_C', '-Ai*(z_C^2)*((sqrt(IS)/(1+sqrt(IS)))-0.3*IS )');
model.component('comp1').variable('var1').set('log_coeff_A', '-Ai*(z_A^2)*((sqrt(IS)/(1+sqrt(IS)))-0.3*IS )');
% model.component('comp1').variable('var1').set('sigma0_a', 'C1*(Psi_0_a - Psi_1_a)');
model.component('comp1').variable('var1').set('sigma0_a', 'C1a*(Psi_0_a - intop1(V))');
model.component('comp1').variable('var1').set('T_sigma0_a', 'sa_F_a*sigma0_a');
model.component('comp1').variable('var1').set('sigma0_b', 'C1b*(Psi_0_b - intop2(V))');
model.component('comp1').variable('var1').set('T_sigma0_b', 'sa_F_b*sigma0_b');

%% Global equations
% ge1
model.component('comp1').physics('ge').feature('ge1').label('IS and log_as');
model.component('comp1').physics('ge').feature('ge1').setIndex('name', 'IS', 0, 0);
model.component('comp1').physics('ge').feature('ge1').setIndex('equation', 'IS - 0.5*((H/mi_0)*(z_H^2)+(OH/mi_0)*(z_OH^2)+(C/mi_0)*(z_C^2)+(A/mi_0)*(z_A^2))', 0, 0);
model.component('comp1').physics('ge').feature('ge1').setIndex('initialValueU', 'T_C [m^3/mol]', 0, 0);
model.component('comp1').physics('ge').feature('ge1').setIndex('initialValueUt', '0', 0, 0);
model.component('comp1').physics('ge').feature('ge1').setIndex('name', 'log_a_OH', 1, 0);
model.component('comp1').physics('ge').feature('ge1').setIndex('equation', 'log_a_OH - (log_coeff_OH + log10(OH/mi_0))', 1, 0);
model.component('comp1').physics('ge').feature('ge1').setIndex('initialValueU', '1', 1, 0);
model.component('comp1').physics('ge').feature('ge1').setIndex('initialValueUt', 0, 1, 0);
model.component('comp1').physics('ge').feature('ge1').setIndex('name', 'log_a_C', 2, 0);
model.component('comp1').physics('ge').feature('ge1').setIndex('equation', 'log_a_C - (log_coeff_C + log10(C/mi_0))', 2, 0);
model.component('comp1').physics('ge').feature('ge1').setIndex('initialValueU', '1', 2, 0);
model.component('comp1').physics('ge').feature('ge1').setIndex('initialValueUt', 0, 2, 0);
model.component('comp1').physics('ge').feature('ge1').setIndex('name', 'log_a_A', 3, 0);
model.component('comp1').physics('ge').feature('ge1').setIndex('equation', 'log_a_A - (log_coeff_A + log10(A/mi_0))', 3, 0);
model.component('comp1').physics('ge').feature('ge1').setIndex('initialValueU', '1', 3, 0);

% ge2
model.component('comp1').physics('ge').create('ge2', 'GlobalEquations', -1);
model.component('comp1').physics('ge').feature('ge2').label('H and mass action laws');
model.component('comp1').physics('ge').feature('ge2').set('CustomDependentVariableUnit', '1');
model.component('comp1').physics('ge').feature('ge2').set('DependentVariableQuantity', 'none');
model.component('comp1').physics('ge').feature('ge2').setIndex('CustomDependentVariableUnit', 'mol/L', 0, 0);
model.component('comp1').physics('ge').feature('ge2').setIndex('name', 'H', 0, 0);
model.component('comp1').physics('ge').feature('ge2').setIndex('equation', 'log_a_H - (log_coeff_H + log10(H/mi_0))', 0, 0);
model.component('comp1').physics('ge').feature('ge2').setIndex('initialValueU', '10^(-pH)', 0, 0);
model.component('comp1').physics('ge').feature('ge2').setIndex('initialValueUt', 0, 0, 0);

model.component('comp1').physics('ge').feature('ge2').setIndex('name', 'OH', 1, 0);
model.component('comp1').physics('ge').feature('ge2').setIndex('equation', 'log_a_OH+log_a_H-log_K_OH', 1, 0);
model.component('comp1').physics('ge').feature('ge2').setIndex('initialValueU', '10^(-(14-pH))', 1, 0);
model.component('comp1').physics('ge').feature('ge2').setIndex('initialValueUt', 0, 1, 0);

model.component('comp1').physics('ge').feature('ge2').setIndex('name', 'TiOH2_a', 2, 0);
model.component('comp1').physics('ge').feature('ge2').setIndex('equation', 'log10(TiOH2_a/mi_0) - log_a_H - log10(exp(-e_kbt*Psi_0_a))-log10(TiOH_a/mi_0)-log_K_TiOH2', 2, 0);
model.component('comp1').physics('ge').feature('ge2').setIndex('initialValueU', 1e-6, 2, 0);
model.component('comp1').physics('ge').feature('ge2').setIndex('initialValueUt', 0, 2, 0);


model.component('comp1').physics('ge').feature('ge2').setIndex('name', 'FeOH2_b', 3, 0);
model.component('comp1').physics('ge').feature('ge2').setIndex('equation', 'log10(FeOH2_b/mi_0) - log_a_H - log10(exp(-e_kbt*Psi_0_b))-log10(FeOH_b/mi_0)-log_K_FeOH2', 3, 0);
model.component('comp1').physics('ge').feature('ge2').setIndex('initialValueU', 1e-6, 3, 0);
model.component('comp1').physics('ge').feature('ge2').setIndex('initialValueUt', 0, 3, 0);

%g3
model.component('comp1').physics('ge').create('ge3', 'GlobalEquations', -1);
model.component('comp1').physics('ge').feature('ge3').label('mass balance');
model.component('comp1').physics('ge').feature('ge3').set('DependentVariableQuantity', 'none');
model.component('comp1').physics('ge').feature('ge3').setIndex('CustomDependentVariableUnit', 'mol/L', 0, 0);

model.component('comp1').physics('ge').feature('ge3').set('SourceTermQuantity', 'none');
model.component('comp1').physics('ge').feature('ge3').setIndex('CustomSourceTermUnit', 'mol/L', 0, 0);

model.component('comp1').physics('ge').feature('ge3').setIndex('name', 'C', 0, 0);
model.component('comp1').physics('ge').feature('ge3').setIndex('equation', 'C - T_C', 0, 0);
model.component('comp1').physics('ge').feature('ge3').setIndex('initialValueU', 'T_C', 0, 0);

model.component('comp1').physics('ge').feature('ge3').setIndex('name', 'A', 1, 0);
model.component('comp1').physics('ge').feature('ge3').setIndex('equation', 'A - T_A', 1, 0);
model.component('comp1').physics('ge').feature('ge3').setIndex('initialValueU', 'T_A', 1, 0);

model.component('comp1').physics('ge').feature('ge3').setIndex('name', 'TiOH_a', 2, 0);
model.component('comp1').physics('ge').feature('ge3').setIndex('equation', 'TiOH2_a + TiOH_a - T_surf_a', 2, 0);
model.component('comp1').physics('ge').feature('ge3').setIndex('initialValueU', 'T_surf_a', 2, 0);

model.component('comp1').physics('ge').feature('ge3').setIndex('name', 'FeOH_b', 3, 0);
model.component('comp1').physics('ge').feature('ge3').setIndex('equation', 'FeOH2_b + FeOH_b - T_surf_b', 3, 0);
model.component('comp1').physics('ge').feature('ge3').setIndex('initialValueU', 'T_surf_b', 3, 0);

%g4
model.component('comp1').physics('ge').create('ge4', 'GlobalEquations', -1);
model.component('comp1').physics('ge').feature('ge4').label('mass balance electrostatic');
model.component('comp1').physics('ge').feature('ge4').set('SourceTermQuantity', 'none');
model.component('comp1').physics('ge').feature('ge4').setIndex('CustomSourceTermUnit', 'mol/L', 0, 0);
model.component('comp1').physics('ge').feature('ge4').set('DependentVariableQuantity', 'none');
model.component('comp1').physics('ge').feature('ge4').setIndex('CustomDependentVariableUnit', 'V', 0, 0);

model.component('comp1').physics('ge').feature('ge4').setIndex('name', 'Psi_0_a', 0, 0);
model.component('comp1').physics('ge').feature('ge4').setIndex('equation', '-(TiOH_a/2) + (TiOH2_a/2)-T_sigma0_a', 0, 0);
model.component('comp1').physics('ge').feature('ge4').setIndex('initialValueU', -2.618154e-02, 0, 0);

model.component('comp1').physics('ge').feature('ge4').setIndex('name', 'Psi_0_b', 1, 0);
model.component('comp1').physics('ge').feature('ge4').setIndex('equation', '-(FeOH_b/2) + (FeOH2_b/2)-T_sigma0_b', 1, 0);
model.component('comp1').physics('ge').feature('ge4').setIndex('initialValueU', -2.618154e-02, 1, 0);

%% Electrostatic

% model.component('comp1').physics('es').create('pot1', 'ElectricPotential', 0);
% model.component('comp1').physics('es').feature('pot1').set('V0', 'Psi_3_a');
% model.component('comp1').physics('es').feature('pot1').selection.set([1]);
% model.component('comp1').physics('es').create('pot2', 'ElectricPotential', 0);
% model.component('comp1').physics('es').feature('pot2').selection.set([2]);
% model.component('comp1').physics('es').feature('pot2').set('V0', 'Psi_D_b');
% model.component('comp1').physics('es').feature('pot2').set('constraintType', 'unidirectionalConstraint');
% model.component('comp1').physics('es').feature('pot1').set('constraintType', 'unidirectionalConstraint');

model.component('comp1').physics('es').feature('ccn1').set('epsilonr_mat', 'userdef');
model.component('comp1').physics('es').feature('ccn1').set('epsilonr', {'epsilon' '0' '0' '0' 'epsilon' '0' '0' '0' 'epsilon'});
model.component('comp1').physics('es').create('scd1', 'SpaceChargeDensity', 1);
model.component('comp1').physics('es').feature('scd1').selection.set([1]);
model.component('comp1').physics('es').feature('scd1').set('rhoq', 'e*NA*(z_H*H*exp(-z_H*e_kbt*V) + z_OH*OH*exp(-z_OH*e_kbt*V)+ z_C*C*exp(-z_C*e_kbt*V) + z_A*A*exp(-z_A*e_kbt*V))');
model.component('comp1').physics('es').create('sfcd1', 'SurfaceChargeDensity', 0);
model.component('comp1').physics('es').feature('sfcd1').selection.all;
model.component('comp1').physics('es').feature('sfcd1').selection.set([1]);
model.component('comp1').physics('es').feature('sfcd1').set('rhoqs', 'sigma0_a');
model.component('comp1').physics('es').create('sfcd2', 'SurfaceChargeDensity', 0);
model.component('comp1').physics('es').feature('sfcd2').selection.all;
model.component('comp1').physics('es').feature('sfcd2').selection.set([2]);
model.component('comp1').physics('es').feature('sfcd2').set('rhoqs', 'sigma0_b');

%% mesh
model.component('comp1').mesh('mesh1').autoMeshSize(1);

%% Study
% model.study('std1').create('param', 'Parametric');
% model.study('std1').feature('param').setIndex('pname', 'pH', 0);
% model.study('std1').feature('param').setIndex('pname', 'extra_C', 1);
% model.study('std1').feature('param').setIndex('pname', 'extra_A', 2);
% model.study('std1').feature('param').setIndex('punit', 'mol/L', 1);
% model.study('std1').feature('param').setIndex('punit', 'mol/L', 2);
% model.study('std1').feature('param').setIndex('plistarr', 'range(2,1,12)', 0);
% %s1= '0.0, 0.0, 0.0, 0.0, 0.0, 3.594e-09, 1.040e-06, 1.022e-05, 1.020e-04';
% model.study('std1').feature('param').setIndex('plistarr', s1, 1);
% %s2 ='1.115e-02, 1.039e-03, 1.017e-04, 1.016e-05, 1.023e-06, 0.0, 0.0, 0.0, 0.0';
% model.study('std1').feature('param').setIndex('plistarr', s2, 2);
model.study('std1').create('param', 'Parametric');
model.study('std1').feature('param').set('paramselect', false);
model.study('std1').feature('param').setIndex('pname', 'Lf', 0);
model.study('std1').feature('param').setIndex('punit', 'nm', 0);
model.study('std1').feature('param').setIndex('plistarr', distance, 0);

model.sol.create('sol1');
model.sol('sol1').study('std1');

model.study('std1').feature('stat').set('notlistsolnum', 1);
model.study('std1').feature('stat').set('notsolnum', '1');
model.study('std1').feature('stat').set('listsolnum', 1);
model.study('std1').feature('stat').set('solnum', '1');

model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').feature('st1').set('study', 'std1');
model.sol('sol1').feature('st1').set('studystep', 'stat');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').feature('v1').set('control', 'stat');
model.sol('sol1').create('s1', 'Stationary');
model.sol('sol1').feature('s1').create('seDef', 'Segregated');
model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').feature('fc1').set('linsolver', 'dDef');
model.sol('sol1').feature('s1').feature.remove('fcDef');
model.sol('sol1').feature('s1').feature.remove('seDef');
model.sol('sol1').attach('std1');

model.batch.create('p1', 'Parametric');
model.batch('p1').study('std1');
model.batch('p1').create('so1', 'Solutionseq');
model.batch('p1').feature('so1').set('seq', 'sol1');
model.batch('p1').feature('so1').set('store', 'on');
model.batch('p1').feature('so1').set('clear', 'on');
model.batch('p1').feature('so1').set('psol', 'none');
model.batch('p1').set('pname', {'Lf'});
model.batch('p1').set('plistarr', {distance});
model.batch('p1').set('sweeptype', 'sparse');
model.batch('p1').set('probesel', 'all');
model.batch('p1').set('probes', {});
model.batch('p1').set('plot', 'off');
model.batch('p1').set('err', 'on');
model.batch('p1').attach('std1');
model.batch('p1').set('control', 'param');

model.sol.create('sol2');
model.sol('sol2').study('std1');
model.sol('sol2').label('Parametric Solutions 1');

model.batch('p1').feature('so1').set('psol', 'sol2');

model.sol('sol1').feature('s1').feature('fc1').set('dtech', 'hnlin');
model.sol('sol1').feature('s1').feature('fc1').set('rstepabs', 1);
model.sol('sol1').feature('s1').feature('fc1').set('rstep', 2);
model.sol('sol1').feature('s1').feature('fc1').set('minsteprecovery', '1e-4');
model.sol('sol1').feature('s1').feature('fc1').set('maxiter', 250);

model.batch('p1').run;



Rip=cell(n_val,1);
Lf = zeros(n_val,1);
IS = zeros(n_val,1);
log_a_H = zeros(n_val,1);
log_a_OH= zeros(n_val,1);
log_a_C= zeros(n_val,1); 
log_a_A= zeros(n_val,1);
H= zeros(n_val,1);
OH= zeros(n_val,1); 
C= zeros(n_val,1);
A= zeros(n_val,1);
TiOH_a= zeros(n_val,1); 
TiOH2_a= zeros(n_val,1); 
FeOH_b= zeros(n_val,1); 
FeOH2_b= zeros(n_val,1);
Psi_0_a= zeros(n_val,1);
sigma0_a= zeros(n_val,1);
Psi_0_b= zeros(n_val,1);
sigma0_b= zeros(n_val,1);
for i = 1:n_val
    Rip{i}=mpheval(model, {'V', 'es.Ex', 'ew','comp1.IS', 'log_a_H', 'comp1.log_a_OH', 'comp1.log_a_C', 'comp1.log_a_A', 'comp1.H', 'comp1.OH', 'comp1.C', 'comp1.A', 'comp1.TiOH_a', 'comp1.TiOH2_a', 'comp1.Psi_0_a', 'comp1.sigma0_a', 'Lf', 'comp1.FeOH_b', 'comp1.FeOH2_b', 'comp1.Psi_0_b', 'comp1.sigma0_b'}, 'dataset', 'dset2','outersolnum',i);
    Lf(i) = Rip{i}.d17(1);
    IS (i)= Rip{i}.d4(1);
    log_a_H(i) = Rip{i}.d5(1);
    log_a_OH(i)= Rip{i}.d6(1);
    log_a_C(i)= Rip{i}.d7(1);
    log_a_A(i)= Rip{i}.d8(1);
    H(i)= Rip{i}.d9(1);
    OH(i)= Rip{i}.d10(1);
    C(i)= Rip{i}.d11(1);
    A(i)= Rip{i}.d12(1);
    TiOH_a(i)= Rip{i}.d13(1);
    TiOH2_a(i)= Rip{i}.d14(1);
    Psi_0_a(i)= Rip{i}.d15(1);
    sigma0_a(i)= Rip{i}.d16(1);
    FeOH_b(i)= Rip{i}.d18(1);
    FeOH2_b(i)= Rip{i}.d19(1);
    Psi_0_b(i)= Rip{i}.d20(1);
    sigma0_b(i)= Rip{i}.d21(1);
    
    
end

vector_distance = Lf;
vector_distance_neg=vector_distance*(-1);
n_distances = length(Lf);
ew = Rip{1}.d3(1);
Kb = 1.3806E-23;% J/K
T = 298.15; %K
e = 1.6022E-19; % C
NA = 6.0221E23; % 1/mol
pH_range=pH*ones(n_distances,1);

R = 8.3145; % J/(mol*K)
F = 96485; %C/mol
Ai = (2*R*T)/F;   % 0.051385V
kappa_v=sqrt((ew.*R.*T)./(2.*IS.*1000.*(NA.*e).^2));  % V_info.d3(1,i) is supposed to be the ionic strength
% xrara=0.65;
% Eq = (exp(-xrara).*(exp(Psi_3_a./Ai)-1))./(exp(Psi_3_a./Ai)+1) ;
%Psi_zeta_1 = atanh(Eq).*(2*Ai);
zH=1;zOH=-1;zC=1;zA=-1;
Osmov = zeros(n_distances,1);
Osmo_integral = zeros(n_distances,1);
Felec = zeros(n_distances, 1);
for i=1:n_distances
    psi=Rip{i}.d1;
    Vx = Rip{i}.d2;

    SigmaD = ew*Vx;
    PsiD_a(i)=psi(1);PsiD_b(i)=psi(end);SigmaD_a(i)=SigmaD(1); SigmaD_b(i)=SigmaD(end);

    Osmo = Kb.*T.*NA.*(H(i).*(exp(-(zH.*e.*psi)./(Kb.*T))-1)+OH(i).*(exp(-(zOH.*e.*psi)./(Kb.*T))-1)+C(i).*(exp(-(zC.*e.*psi)./(Kb.*T))-1)+A(i).*(exp(-(zA.*e.*psi)./(Kb.*T))-1));

    Osmo = Osmo-(ew./2).*(Vx.^2);
        Osmov(i) = Osmo(round(length(Osmo)/2));
        if i==1
            Osmo_integral(i)=Osmov(i);
        else
            Osmo_integral(i)=trapz (vector_distance_neg(1:i),Osmov(1:i)');
            Felec(i)=2.*pi.*Osmo_integral(i);
        end

end
T = table(Lf, Osmov, Felec,  pH_range,kappa_v,IS, log_a_H, log_a_OH, log_a_C, log_a_A, H/1e3, OH/1e3, C/1e3, A/1e3, TiOH_a/1e3, TiOH2_a/1e3, FeOH_b/1e3, FeOH2_b/1e3,Psi_0_a, sigma0_a, Psi_0_b, sigma0_b, PsiD_a', PsiD_b',SigmaD_a', SigmaD_b');
filename =strcat('Data_Bench_RutileGoethite_','pH_',num2str(pH),'.xlsx');
T.Properties.VariableNames = {'distance', 'Osmo', 'Felectro','pH','Kappa', 'IS', 'log_a_H', 'log_a_OH', 'log_a_C', 'log_a_A', 'H', 'OH', 'C', 'A', 'TiOH_a', 'TiOH2_a','FeOH_b', 'FeOH2_b', 'Psi_0_a', 'sigma0_a','Psi_0_b', 'sigma0_b','PsiD_a', 'PsiD_b','SigmaD_a2', 'SigmaD_b2'};
writetable(T, filename, 'Sheet',1)
end