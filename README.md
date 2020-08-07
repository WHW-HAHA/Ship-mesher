# Ship mesher

## Introduction

This project was developed as part of my master thesis. Used to be Matlab version only, now the python version is underdevelopment.

The purpose: The scripts first reads typical ship model files, which each contains the scatter cloud of a typical ship frame(part beneath the water). Then convert this scatter cloud into a mesh file in vtk format. Associating with the wave modelling software, the vectors of wave static and dynamic pressure for each panel can be calculated to get the wave forces on the whole ship.   

## Typical results
* Typical scatter cloud
![pic1](https://user-images.githubusercontent.com/43483189/89647882-9ac1cf00-d8be-11ea-8067-52843b38ab4b.png)
* Typical ship mesh
![pic1](https://user-images.githubusercontent.com/43483189/89647960-b7f69d80-d8be-11ea-91a1-23c227bd333e.png)

* Typical visualizations in ParaView 
![pic1](https://user-images.githubusercontent.com/43483189/89648162-14f25380-d8bf-11ea-87fd-99e49da65c23.png)
![pic2](https://user-images.githubusercontent.com/43483189/89648174-1c196180-d8bf-11ea-9217-ac25e72bbeec.png)
