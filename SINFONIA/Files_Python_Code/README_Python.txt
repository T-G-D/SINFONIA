The following folders contain the scripts for benchmarks 1a and 3c, written with Python 2.7.
To run the code, the following packages are required:

- Python (https://www.python.org/, version 2.7).
- Numpy (http://www.numpy.org/)
- Scipy (https://www.scipy.org/)
- matplotlib (http://matplotlib.org/)
- h5py (https://www.h5py.org/)

The main, running script is "NR_chem_xx.py", relying on the other two (bvp.py and NR_solver_xx.py) 
to perform all the calculations. Only modifications in "NR_chem_xx.py" scripts (i.e., # User input # 
and # Model definition # parts) need to be modified in order to simulate other chemical conditions.

NR_chem_1a.py is defined for a BSM. 
Two outputs can be selected:
- a global .dat file containing all calculations except diffuse layer potentials
- individual .dat files for every pH and distance computed including respective diffuse layer calculations

NR_chem_3c.py is defined for a FLM and all outputs are collected in a HDF5 library. 
The script is defined to calculate variable ionic strengths in one run too if necessary.
Specific print outs in the form of .dat files have to be selected individually from 
"Output_overall_prints_3c.py" and "NR_Force_3c.py". Calculating different geometries for charge regulation
force-distance curves and normalized forces (e.g., in N/m) is also possible when specified in the script.


####################################################################################################
DESCRIPTION OF THE REQUIRED INPUT
####################################################################################################
## Inside the NR_chem_xx.py script:

The user can modify the following characteristics of the chemical system within the sections: 
"USER INPUT"
- pH range 
- ionic strength
- pH corrections for main cations/anions
- surface separation distances
- in benchmark 1a, force-distance specifications are also defined in "NR_chem_1a.py"

"MODEL DEFINITIONS"
- involved master species / components
- component charges
- involved product species, solution + surface
- intrinsic K values
- surface properties such as capacitances, density sites, solid/liquid ratio and specific surface area.
- A/B matrices, defining chemical reactions in the system 

After this, the code can run alone, based on the following core routines:
"EXECUTION OF MODEL CALCULATIONS" 
The code will loop automatically to calculate all specified conditions according to the 
previously defined ionic strengths, pH values and surface distances. In this part of the script,
the initial vectors are constructed based on previous definitions. The actual problem solving is performed 
by creating an instance of the simplex_point class and and calling NR_fit (Newton-Raphson), the core routines of the SINFONIA code,
which are imported from the NR_solver_xx.py module. 
The NR_chem_xx.py script further retrieves the output values and prints a summary in the python shell. 
In benchmark 1a, the script produces automatically an .dat file with the output data. For the case of benchmark 3c, 
the script stores the data in the h5py library. In order to retrieve this data, the output .dat file can be 
obtained by using the "Output_overall_prints.py" script. In order to obtain information about force-distances in 
benchmark 3c, the user needs to run "NR_Force_3c.py".


## Inside the "NR_Force_3c.py" script:

The user can specify the geometry, hamaker constant and tip radius for force-distance calculations
in the "USER INPUT" section. The script will collect the data from the h5py library and perform
the corresponding calculations for every ionic strength, pH value and distances collected in the library.
Only calculations for charge regulation are specified in the script, but appropriate modifications can be
included in order to obtain other results such as forces based on constant charge and constant potential assumptions.
The results are collected in several lists where they are both (i) printed as an output figure, automatically stored as a .png,
and (ii) retrieved in an output .dat file, for every ionic strength (labelled accordingly).


## Inside the NR_solver_xx.py script:
"simplex_point"
- definition of the Gauss' law for charge regulation within the calculations of the electrostatic terms
- definition of activity corrections using Davies equation
- definition of mass balance conditions

"scipy PB solver"
- solver for the Poisson-Boltzmann equation based on the "bvp.py" code defined in python for solving ODE
- collection of electrostatic potential in the diffuse layer
- calculations of the osmotic pressure 

"NR main routine"
The mass balance and electrostatic term definitions are called in this part of the script to perform the iterations
of the Newton-Raphson method, based on the Jacobian matrix from Westall 1980 modified for charge regulation 
applications with dissimilar opposing surfaces.

####################################################################################################

