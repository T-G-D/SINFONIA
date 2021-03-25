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