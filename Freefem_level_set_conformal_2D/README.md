# Freefem++ platform

Freefem++ platform for topology optimization via the level-set method.
At each iteration the mesh edges coincide with the underlying shape given by the level-set.
The objective function is a linear combination of the structure compliance and its volume.

## Installation (software needed)
1. Freefem++ (http://www.freefem.org/)
2. Python 2.7
3. MMG3D/2D (https://www.mmgtools.org/)
4. Advect & MshDist packages (https://github.com/loicNorgeot/ICStoolbox)

## Usage
1. Compile de C++ file "Delgado.cpp" with the command $ ff-c++ Delgado.cpp
2. Run $ Freefem++ Main.edp
3. Plot the convergence curve with $ python read_plot.py

## Credits

Created by Gabriel Delgado 

## Changelog

### V1.0