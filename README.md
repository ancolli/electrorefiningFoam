# v1.0.0

Cite as: A.N. Colli and J.M. Bisang.
DOI: 

# electrorefiningFoam
It is described how to simulate current distributions in electrochemical reactors under natural convection flows using OpenFOAM. The process of pre-processing, running, and post-processing a basic 2D case using the developed solver electroNaturalConvectionPimpleFoam and codedMixed boundary condition (BC) is demonstrated. The proposed strategy is fully functional in OpenFOAM version 7.

# Disclaimer
This offering is not approved or endorsed by OpeFOAM Foundation, producer and distributor of the OpenFOAM software via www.openfoam.org.

# Usage
In applications (A) you will find the scripts to compile the solver and post-processing utilities in order to solve current distribution in natural convection systems.
In tutorial (B) you will find an example of a 2D parallel-plate cell. 

# #  A) Applications
**1.**  Solver  
_A)_ Paste applications/utilities/solvers/electroNaturalConvectionPimpleFoam inside OpenFOAM user directory (Applications/Utilities/Solvers).  
_B)_ Open a terminal inside electroNaturalConvectionPimpleFoam.  
_C)_ Run wmake.  
**2.**  Total current density  
_A)_ Paste applications/utilities/postProcessing/wallFlux inside OpenFOAM user directory (Applications/Utilities/postProcessing).  
_B)_ Open a terminal inside wallFlux.  
_C)_ Run wmake.  
**3.**  Total averaged (over time) current density  
_A)_ Paste applications/utilities/postProcessing/potentialWallFluxAverage inside OpenFOAM user directory (Applications/Utilities/postProcessing).  
_B)_ Open a terminal inside potentialWallFluxAverage.  
_C)_ Run wmake.  


# #  B) Tutorial
**1-** Paste tutorial inside OpenFOAM user directory (Run/Tutorials).  
**2-** Enter to tutorial and open a Terminal.  
**3-** Modify transport properties (conductivity, kinetic parameters, equilibrium potentials, etc.) inside constant/transportProperties.  
**4-** Modify control properties (kind of electric control, under-over relaxation parameter, maxPot, tolerances) inside constant/controlProperties.   
**5-** Run ./Allrun.    

