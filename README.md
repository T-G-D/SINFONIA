# SINFONIA stands for Speciation and Interparticle Forces for Nanoscale Interactions.

SINFONIA is implemented in COMSOL and Python for chemical speciation and force-distance calculations 
of solid/liquid interfaces based on charge regulation. The provided examples represent two different Stern-Layer 
compartimentalizations regarding the Surface Complexation Model. This repository contains the associated scripts 
of both COMSOL and Python for two Benchmark cases:

- Benchmark 1a: protonation of two interacting silica surfaces in 0.01M NaCl and a wide pH range 
		(solution composition is corrected for acid base additions necessary to achieve a certain pH), 
                A Basic Stern Model (BSM) setting of the code is presented. It is however applied in a way such that it
		mimicks a Double Layer Model (DLM) approach.

- Benchmark 3c: protonation of interacting rutile and goethite surfaces in 0.001M NaCl at selected 
		pH values (without pH correction). A Four Layer Model (FLM) setting of the code is provided,
		though applied as a BSM.

These scripts serve as a guideline that can be used to help and inspire anybody, wanting to do similar 
types of simulations. Anybody using the provided scripts as a starting point for further calculations, 
please give proper credit to the original publication, soon to be specified.
