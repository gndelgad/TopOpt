# Freefem++ 3D platform

Freefem++ platform for 3D topology optimization via the level-set method.
The objective function is a linear combination of the structure compliance and its volume.

## Installation (software needed)
1. Freefem++ (http://www.freefem.org/)
2. MMG3D/2D (https://www.mmgtools.org/)
3. Python 2.7

## Usage
1. Compile de C++ file '"Delgado3D.cpp"' with the command '$ ff-c++ Delgado3D.cpp'
2. Generate the desired mesh of the design domain (rectangular parallelepiped with an isotropic structured mesh) 
via $ Freefem++ structured_mesh_constructor.edp. The resulting mesh is called '"IsoCube.mesh"'.
The size of the design domain can be modified at the end of the script.
3. Run the optimization algorithm with '$ Freefem++ Main.edp'
4. Visualize the shape at each iteration with xd3d. For that purpose do: '$ ./xd3d_visu/xd3d_new IsoCube.mesh'
5. Once the mesh opened, go in the xd3d GUI to "val", select ../LevelSet.psi and OK. Then press "Isosurface" and OK.
6. In order to generate a conformal mesh to the final shape, do '$ Freefem++ read_medit_MMG3D.edp'
7. Plot the convergence curve with $ python read_plot.py

## Comments
1. The '"Main.edp"' script defines: optimization loop + definitions of material properties + characteristics of the mesh
2. The above mentioned characteristics must coincide with the ones defined in '"structured_mesh_constructor.edp"'
3. The '"AuxiliarBasic.edp"' script defines: boundary conditions + physical equations + functional spaces + HJ equations.
4. No need to modify "shift3d.edp". This script drives the shift matrix construction for the finite difference part.
5. At the beginning of "Main.edp" script, the MUMPS library is needed to handle large matrices in parallel. 
In Linux use the line 'load "MUMPS"' meanwhile in Windows rather 'load "MUMPS_seq"'. 

## Credits

Created by Gabriel Delgado 

## Changelog

### V1.0