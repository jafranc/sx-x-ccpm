*C*urvature *C*apture in *P*orous *M*edia is a from image curvature measuring tool.
It uses a binarized mask to extract 3d-bounding surface meshes and measure local curvature on them.

## How to install

It can be installed from sources cloning it and submodules. It relies on *cmake* for makefile generation.
It requires *Eigen3*, *Boost-filesystem* and *libtiff* to be installed on the system.


```bash
	git clone https://gitlab.com/jacquesn7/ccpm
	cd ccpm
	git submodule init 
	git submodule update
	mkdir build-debug && cd build-debug/
	cmake ..
	make [-j nprocs]
```


## How to use

```bash

	path/to/bin/ccpm --image path/to/images.tiff --ard=15.,2,.5 -o /path/to/ccpm/generated/data/
```

`ard` arguments stands for angle, radius and distance that are input arguments from the surface mesh extraction algorithms.
The run will produce 2 files : a STL of the surface mesh and a CSV file in XYZ format with primary curvature k1 and secondary
curvature k2 that can be displayed in Paraview through `TableToPoint` filter.



