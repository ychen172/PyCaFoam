# PyCaFOAM personal 2d CFD solver
Content:
- CTD: computational thermodynamics (conduction solver)
- IncStePISO: incompressible steady state piso solver (SD version is using steepest descent as oppose to direct method in solving Ax=b)

### Fully developed viscous channel flow
Setup for this case can be found in IncStPISO.

OpenFOAM solution:
![OF Solve](doc/OFUFullDev.png)

PyCaFOAM solution:
![PCU Solve](doc/PCFUFullDev.jpg)

Also the v velocity and pressure: 
![PCV Solve](doc/PCFVFullDev.jpg)
![PCP Solve](doc/PCFPFullDev.jpg)

### Jet shear layer with air
Because did not implement rhie chow interpolation, currently suffer checker board problem

However, the development of shear layer can still be vaguely seen

Here are the U velocity V velocity and pressure:
![MixU Solve](doc/UvelMixLayer.jpg)
![MixV Solve](doc/VvelMixlLayer.jpg)
![MixP Solve](doc/PMixLayer.jpg)

### Jet shear layer with air and non-uniform grid
As can be seen checker board can be alleviated using nonuniform grid
![MixUSm Solve](doc/USmooth.jpg)
![MixVSm Solve](doc/VSmooth.jpg)
![MixPSm Solve](doc/PSmooth.jpg)

### Problem suffered
- cannot accept reverse flow shows NaN big problem
- relies on non-uniform grid to kill the checker board problem
- advection is first order discretization, inconsistent with others
- mometum prediction is [U,V] solved together, higher cost then pressure correction, need to apply ADI

### Formulation
Structural grid is generated as follows
![Grid](doc/Grid.png)
The four colors represents the mapping between cell --> starting vertex --> X direction face --> Y direction face

Here is a incompressible steady state equation to solve
![MomEqu](doc/MomEqu.png)
Pressure and viscosity are both normalized by density
Pressure is also modeled as a source term in consistent with the pressure corrector equation (It does not matter if linear interpolation is used)

Formulation of advection term:
![Advection](doc/Advection.png)
First order upwind 
Simultaneously for the sumPhi matrix which will be used to implement div(mom equ) later when deriving pressure correction equ

Formulation of diffusion term:
![Diffusion](doc/Diffusion.png)
Linear interpolation used. Interpolation based on volume ratios

Formulation of pressure term:
![Pressure](doc/Pressure.png)
Linear

Formulation of discretized momentum equation:
![TotalMomDis](doc/TotalMomDis.png)

Formulation of PISO loop:
![PISO](doc/PISO.png)

