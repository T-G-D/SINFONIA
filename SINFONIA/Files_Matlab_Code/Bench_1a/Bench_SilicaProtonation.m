function [T] = Bench_SilicaProtonation()
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
model.param.set('pH', '6');
model.param.set('fix_I', '0.01 [mol/L]');
model.param.set('fix_log_lambda', '-0.05');
model.param.set('log_K_surf', '-6.81');
model.param.set('H', '10^-(pH)*mi_0');
model.param.set('OH', '(10^log_Kw)/(H)*mi_0*mi_0');
model.param.set('log_Kw', '-14');

model.param.set('Z_H', '1');
model.param.set('Z_OH', '-1');
model.param.set('Z_Na', '1');
model.param.set('Z_Cl', '-1');
model.param.set('Na_noextra', '0.01[mol/L]');
model.param.set('Cl_noextra', '0.01 [mol/L]');
model.param.set('extra_Na', '0 [mol/L]');
model.param.set('extra_Cl', '0 [mol/L]');
model.param.set('Cl', 'Cl_noextra + extra_Cl');
model.param.set('Na', 'Na_noextra + extra_Na');

model.param.set('sites_a', '4.6 [1/nm^2]');
model.param.set('surf_a', '1 [m^2/g]');
model.param.set('mass_a', '1 [g/L]');
model.param.set('sa_F_a', '(surf_a*mass_a)/F');
model.param.set('Ai', '0.5092');
model.param.set('F_RT', 'F/(R*T)');
model.param.set('C1a', '1e4 [F/m^2]');
model.param.set('T_surf', '(sites_a*surf_a*mass_a)/NA');
%% Geometry
model.component('comp1').geom('geom1').create('i1', 'Interval');
model.component('comp1').geom('geom1').feature('i1').setIndex('coord', 'L0', 0);
model.component('comp1').geom('geom1').feature('i1').setIndex('coord', 'Lf', 1);
model.component('comp1').geom('geom1').run('fin');

%% Couplings
model.component('comp1').cpl.create('intop1', 'Integration');
model.component('comp1').cpl('intop1').set('axisym', true);
model.component('comp1').cpl('intop1').selection.geom('geom1', 0);
model.component('comp1').cpl('intop1').selection.set([1]);

%% Variables
model.component('comp1').variable.create('var1');
model.component('comp1').variable('var1').set('sigma0_a', 'C1a*(Psi_0_a - intop1(V))');
model.component('comp1').variable('var1').set('T_sigma0_a', 'sa_F_a*sigma0_a');

%% Global equations
model.component('comp1').physics('ge').feature('ge1').label('mass acition laws');
model.component('comp1').physics('ge').feature('ge1').label('H and mass action laws');
model.component('comp1').physics('ge').feature('ge1').set('CustomDependentVariableUnit', '1');
model.component('comp1').physics('ge').feature('ge1').set('DependentVariableQuantity', 'none');
model.component('comp1').physics('ge').feature('ge1').setIndex('CustomDependentVariableUnit', 'mol/L', 0, 0);

model.component('comp1').physics('ge').feature('ge1').setIndex('name', 'SiO_a', 0, 0);
model.component('comp1').physics('ge').feature('ge1').setIndex('equation', 'log10(SiO_a/mi_0)+log10(H/mi_0)+log10(exp(-e_kbt*Psi_0_a))-log_K_surf-log10(SiOH_a/mi_0)', 0, 0);
model.component('comp1').physics('ge').feature('ge1').setIndex('initialValueU', '7.6374e-6 [mol/L]', 0, 0);
model.component('comp1').physics('ge').feature('ge1').setIndex('initialValueUt', 0, 0, 0);
model.component('comp1').physics('ge').feature('ge1').setIndex('description', '', 0, 0);

%g2
model.component('comp1').physics('ge').create('ge2', 'GlobalEquations', -1);
model.component('comp1').physics('ge').feature('ge2').label('mass balance');
model.component('comp1').physics('ge').feature('ge2').set('DependentVariableQuantity', 'none');
model.component('comp1').physics('ge').feature('ge2').setIndex('CustomDependentVariableUnit', 'mol/L', 0, 0);

model.component('comp1').physics('ge').feature('ge2').set('SourceTermQuantity', 'none');
model.component('comp1').physics('ge').feature('ge2').setIndex('CustomSourceTermUnit', 'mol/L', 0, 0);

model.component('comp1').physics('ge').feature('ge2').setIndex('name', 'SiOH_a', 0, 0);
model.component('comp1').physics('ge').feature('ge2').setIndex('equation', 'SiOH_a + SiO_a -T_surf', 0, 0);
model.component('comp1').physics('ge').feature('ge2').setIndex('initialValueU', '7.6374e-6 [mol/L]', 0, 0);
model.component('comp1').physics('ge').feature('ge2').setIndex('initialValueUt', 0, 0, 0);
model.component('comp1').physics('ge').feature('ge2').setIndex('description', '', 0, 0);

%g3

model.component('comp1').physics('ge').create('ge3', 'GlobalEquations', -1);
model.component('comp1').physics('ge').feature('ge3').label('mass balance electrostatic');
model.component('comp1').physics('ge').feature('ge3').set('SourceTermQuantity', 'none');
model.component('comp1').physics('ge').feature('ge3').setIndex('CustomSourceTermUnit', 'mol/L', 0, 0);
model.component('comp1').physics('ge').feature('ge3').set('DependentVariableQuantity', 'none');
model.component('comp1').physics('ge').feature('ge3').setIndex('CustomDependentVariableUnit', 'V', 0, 0);

model.component('comp1').physics('ge').feature('ge3').setIndex('name', 'Psi_0_a', 0, 0);
model.component('comp1').physics('ge').feature('ge3').setIndex('equation', '-SiO_a-T_sigma0_a', 0, 0);
model.component('comp1').physics('ge').feature('ge3').setIndex('initialValueU', -2.618154e-02, 0, 0);

%% Electrostatic

model.component('comp1').physics('es').feature('ccn1').set('epsilonr_mat', 'userdef');
model.component('comp1').physics('es').feature('ccn1').set('epsilonr', {'epsilon' '0' '0' '0' 'epsilon' '0' '0' '0' 'epsilon'});
model.component('comp1').physics('es').create('scd1', 'SpaceChargeDensity', 1);
model.component('comp1').physics('es').feature('scd1').selection.set([1]);
model.component('comp1').physics('es').feature('scd1').set('rhoq', 'e*NA*(H*Z_H*exp(-e_kbt*Z_H*V)+OH*Z_OH*exp(-e_kbt*Z_OH*V) +Na*Z_Na*exp(-e_kbt*Z_Na*V)+Cl*Z_Cl*exp(-e_kbt*Z_Cl*V))');
model.component('comp1').physics('es').create('sfcd1', 'SurfaceChargeDensity', 0);
model.component('comp1').physics('es').feature('sfcd1').selection.all;
model.component('comp1').physics('es').feature('sfcd1').selection.set([1]);
model.component('comp1').physics('es').feature('sfcd1').set('rhoqs', 'sigma0_a');
model.component('comp1').physics('es').create('sfcd2', 'SurfaceChargeDensity', 0);
model.component('comp1').physics('es').feature('sfcd2').selection.all;
model.component('comp1').physics('es').feature('sfcd2').selection.set([2]);
model.component('comp1').physics('es').feature('sfcd2').set('rhoqs', 'sigma0_a');

%% mesh
model.component('comp1').mesh('mesh1').autoMeshSize(1);
%% Study
model.study('std1').create('param', 'Parametric');
model.study('std1').feature('param').setIndex('pname', 'Lf', 0);
model.study('std1').feature('param').setIndex('plistarr', '28, 25, 23, 20, 18, 15, 13, 10, 7, 5, 4.5, 4, 3.5, 3, 2.5, 2, 1.75, 1.5, 1.25, 1.1, 1, 0.9, 0.8, 0.7, 0.6, 0.55, 0.5, 0.45, 0.4, 0.3, 0.2, 0.1', 0);
model.study('std1').feature('param').setIndex('punit', 'nm', 0);
model.study('std1').feature('param').set('pdistrib', true);

model.study('std1').create('param2', 'Parametric');
model.study('std1').feature.move('param', 0);  % move pos
model.study('std1').feature('param2').setIndex('pname', 'pH', 0);
model.study('std1').feature('param2').setIndex('plistarr', 'range(3,1,8)', 0);
model.study('std1').feature('param2').setIndex('pname', 'extra_Na', 1);
model.study('std1').feature('param2').setIndex('plistarr', '0., 0., 0., 0., 3.400e-07, 1.834e-06', 1);
model.study('std1').feature('param2').setIndex('punit', 'mol/L', 1);
model.study('std1').feature('param2').setIndex('pname', 'extra_Cl', 2);
model.study('std1').feature('param2').setIndex('plistarr', '1.114e-03, 1.109e-04, 1.104e-05, 9.500e-07, 0., 0.', 2);
model.study('std1').feature('param2').setIndex('punit', 'mol/L', 2);


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
model.sol('sol1').feature('s1').create('p1', 'Parametric');
model.sol('sol1').feature('s1').feature('p1').set('pname', {'pH' 'extra_Na' 'extra_Cl'});
model.sol('sol1').feature('s1').feature('p1').set('plistarr', {'range(3,1,8)' '0., 0., 0., 0., 3.400e-07, 1.834e-06' '1.114e-03, 1.109e-04, 1.104e-05, 9.500e-07, 0., 0.'});
model.sol('sol1').feature('s1').feature('p1').set('punit', {'' 'mol/L' 'mol/L'});
model.sol('sol1').feature('s1').feature('p1').set('sweeptype', 'sparse');
model.sol('sol1').feature('s1').feature('p1').set('preusesol', 'no');
model.sol('sol1').feature('s1').feature('p1').set('pcontinuationmode', 'no');
model.sol('sol1').feature('s1').feature('p1').set('plot', 'off');
model.sol('sol1').feature('s1').feature('p1').set('plotgroup', 'Default');
model.sol('sol1').feature('s1').feature('p1').set('probesel', 'all');
model.sol('sol1').feature('s1').feature('p1').set('probes', {});
model.sol('sol1').feature('s1').feature('p1').set('control', 'param2');
model.sol('sol1').feature('s1').set('control', 'stat');
model.sol('sol1').feature('s1').create('seDef', 'Segregated');
model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').feature('fc1').set('linsolver', 'dDef');
model.sol('sol1').feature('s1').feature.remove('fcDef');
model.sol('sol1').feature('s1').feature.remove('seDef');
model.sol('sol1').attach('std1');
model.sol('sol1').feature('s1').feature('fc1').set('dtech', 'hnlin');
model.sol('sol1').feature('s1').feature('fc1').set('rstep', 2);
model.sol('sol1').feature('s1').feature('fc1').set('minsteprecovery', '1e-4');
model.sol('sol1').feature('s1').feature('fc1').set('maxiter', 250);
model.sol('sol1').feature('s1').set('stol', '1e-12');

model.batch.create('p1', 'Parametric');
model.batch('p1').study('std1');
model.batch('p1').create('so1', 'Solutionseq');
model.batch('p1').feature('so1').set('seq', 'sol1');
model.batch('p1').feature('so1').set('store', 'on');
model.batch('p1').feature('so1').set('clear', 'on');
model.batch('p1').feature('so1').set('psol', 'none');
model.batch('p1').set('pname', {'Lf'});
model.batch('p1').set('plistarr', {'28, 25, 23, 20, 18, 15, 13, 10, 7, 5, 4.5, 4, 3.5, 3, 2.5, 2, 1.75, 1.5, 1.25, 1.1, 1, 0.9, 0.8, 0.7, 0.6, 0.55, 0.5, 0.45, 0.4, 0.3, 0.2, 0.1'});
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

model.batch('p1').run;



info = mphxmeshinfo(model);

V_info=mpheval(model, {'V', 'es.Ex', 'H', 'Na', 'Cl', 'OH', 'SiOH_a', 'SiO_a', 'comp1.Psi_0_a', 'comp1.sigma0_a'}, 'dataset', 'dset2','outersolnum','all');
[~,~,n_distances] = size(V_info.d1);
%% Force Van der wal vector
% P= -(H*r)/6*h^2
%para
vector_distance = [28, 25, 23, 20, 18, 15, 13, 10, 7, 5, 4.5, 4, 3.5, 3, 2.5, 2, 1.75, 1.5, 1.25, 1.1, 1, 0.9, 0.8, 0.7, 0.6, 0.55, 0.5, 0.45, 0.4, 0.3, 0.2, 0.1]*1e-9; %nm is the h


r = 1; %m
r2=10; %m
H = 1.2e-20; %J

Fvdw_sp = -(H.*r)./(6.*vector_distance.^2);
Fvdw_cc = -(H.*sqrt(r.*r2))./(6.*vector_distance.^2);
Fvdw_ss = (-H./(6.*vector_distance.^2))*((r.*r2)/(r+r2));

%% Calculating OSMO
H = [1, 0.1, 0.01, 0.001, 0.0001, 0.00001]; %mol/m3 1
OH = [1e-8,1e-7, 1e-6, 1e-5, 1e-4, 1e-3 ]; %mol/m3 -1
Cl = [11.114,10.1109, 10.01104, 10.00095, 10, 10]; %mol/m3 -1
Na = [10,10,10,10,10.00034, 10.001834]; %mol/m3 1

zH=1;zOH=-1;zNa=1;zCl=-1;
Kb = 1.3806E-23;% J/K
T = 298.15; %K
e = 1.6022E-19; % C
NA = 6.0221E23; % 1/mol
Osmov = zeros(n_distances,6);

e0 = 8.8542E-12;% F/m
er=78.45203739768931;
ew=er*e0;
R = 8.3145;
fix_I=0.01;
kappa = sqrt((ew*R*T)/(2*fix_I*1000*(NA*e)^2));
kappa_distance = (vector_distance.*kappa)./2;
Osmo_integral = zeros(n_distances,6);
Force_sp = zeros(n_distances, 6);
Force_ss = zeros(n_distances, 6);
Force_cc = zeros(n_distances, 6);
Felec = zeros(n_distances, 6);
vector_distance_neg=vector_distance*(-1);
% Boundaries
for i = 1:n_distances
    for j= 1:6
        psi = V_info.d1(j,:,i);
        Vx = V_info.d2(j,:,i);
        
        Osmo = Kb.*T.*NA.*(H(j).*(exp(-(zH.*e.*psi)./(Kb.*T))-1)+OH(j).*(exp(-(zOH.*e.*psi)./(Kb.*T))-1)+Na(j).*(exp(-(zNa.*e.*psi)./(Kb.*T))-1)+Cl(j).*(exp(-(zCl.*e.*psi)./(Kb.*T))-1));
        Osmo = Osmo-(ew./2).*(Vx.^2);
        Osmov(i,j) = Osmo(round(length(Osmo)/2));
        if i==1
            Osmo_integral(i,j)=Osmov(i,j);
        else
            Osmo_integral(i,j)=trapz (vector_distance_neg(1:i),Osmov(1:i,j)');
            Felec(i,j)=2.*pi.*Osmo_integral(i,j);
            Force_sp(i,j)= Felec(i,j)+Fvdw_sp(i);
            Force_ss(i,j)= Felec(i,j)+Fvdw_ss(i);
            Force_cc(i,j)= Felec(i,j)+Fvdw_cc(i);
        end
    end
end

%% Calculating Integral OSMO
pH_range=3:1:8;
% Table for exporting and comparing results
V_dis=[];K_dis=[];pHseq=[];PsiD_a=[];PsiD_b=[];sigD_a=[];sigD_b=[];Osmo_v=[];Osmo_Int =[];
Fvdw_list_sp =[]; Fvdw_list_ss =[]; Fvdw_list_cc =[];
Fsp=[]; Fss=[]; Fcc=[];
H_v=[]; Na_v=[]; Cl_v=[]; OH_v=[];
SiOH_a_V=[]; Psi0_a_V=[]; SiO_a_V=[]; Sigma0_a_V=[];
% 'V', 'es.Ex', 'H', 'Na', 'Cl', 'OH', 'SiOH_a', 'SiOH_b', 'SiO_a', 'SiO_b'

for i=1:length(pH_range)
V_dis = [V_dis, vector_distance];
K_dis = [K_dis, kappa_distance];
pHseq = [pHseq, pH_range(i)*ones(1,length(vector_distance))];
PsiD_a = [PsiD_a; squeeze(V_info.d1(i,1,:))];
PsiD_b = [PsiD_b; squeeze(V_info.d1(i,end,:))];
sigD_a = [sigD_a; ew*squeeze(V_info.d2(i,1,:))];
sigD_b = [sigD_b; ew*squeeze(V_info.d2(i,end,:))];
Osmo_v = [Osmo_v; Osmov(:,i)];
Osmo_Int = [Osmo_Int;Osmo_integral(:,i)];
Fvdw_list_sp =[Fvdw_list_sp, Fvdw_sp];
Fvdw_list_ss =[Fvdw_list_ss, Fvdw_cc]; 
Fvdw_list_cc =[Fvdw_list_cc, Fvdw_ss];
Fsp=[Fsp; Force_sp(:,i)]; 
Fss=[Fss; Force_ss(:,i)]; 
Fcc=[Fcc; Force_cc(:,i)];
H_v = [H_v; squeeze(V_info.d3(i,1,:))];
Na_v =[Na_v; squeeze(V_info.d4(i,1,:))]; 
Cl_v =[Cl_v; squeeze(V_info.d5(i,1,:))]; 
OH_v =[OH_v; squeeze(V_info.d6(i,1,:))];
SiOH_a_V=[SiOH_a_V; squeeze(V_info.d7(i,1,:))]; 
SiO_a_V=[SiO_a_V; squeeze(V_info.d8(i,1,:))]; 

Psi0_a_V= [Psi0_a_V; squeeze(V_info.d9(i,1,:))]; 
Sigma0_a_V= [Sigma0_a_V; squeeze(V_info.d10(i,1,:))]; 
end

T = table(V_dis', K_dis', pHseq', PsiD_a, PsiD_b, sigD_a, sigD_b, Osmo_v, Osmo_Int, Fvdw_list_sp', Fvdw_list_ss', Fvdw_list_cc', Fsp, Fss, Fcc, H_v, Na_v, Cl_v, OH_v, SiOH_a_V, SiO_a_V, Psi0_a_V,  Sigma0_a_V);
T.Properties.VariableNames = {'distance', 'kh_2', 'pH', 'PsiD_a', 'PsiD_b', 'sigD_a', 'sigD_b', 'Osmotic_pressure', 'Integral_Osmotic_pressure', 'Force_vdW_sphere_plane', 'Force_vdW_sphere_sphere', 'Force_vdW_cylinder_cylinder', 'Force_sphere_plane', 'Force_sphere_sphere', 'Force_cylinder_cylinder', 'H', 'Na', 'Cl', 'OH', 'SiOH_a', 'SiO_a', 'Psi0a', 'Sigma0a'};
T.Properties.VariableUnits={'m','-', '-', 'V', 'V', 'C/m2', 'C/m2', 'Pa', 'Pa*m', 'N', 'N', 'N', 'N', 'N', 'N', 'mol/m3', 'mol/m3', 'mol/m3', 'mol/m3', 'mol/m3', 'mol/m3', 'V', 'C/m2'};

filename ='Bench_SilicaProtonation.xlsx';
writetable(T, filename, 'Sheet',1);

end