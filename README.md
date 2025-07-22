[*C*urvature *C*apture in *P*orous *M*edia](https://jafranc.github.io/sx-x-ccpm/) is a from image curvature measuring tool.
It uses a binarized mask to extract 3d-bounding surface meshes and measure local curvature on them.

## How to install

It can be installed from sources cloning it and submodules. It relies on *cmake* for makefile generation.
It requires *Eigen3*, *Boost* and *libtiff* to be installed on the system.


```bash
	git clone https://github.com/jafranc/sx-x-ccpm
	cd sx-x-ccpm
	git submodule init 
	git submodule update
	mkdir build-debug && cd build-debug/
	cmake ..
	make [-j nprocs]
```


## How to use

```bash
	path/to/bin/ccpm --image path/to/images.tiff [-c 4] -v 85,170,255 -o /path/to/ccpm/generated/data/
	path/to/bin/ccpm --image /path/to/ccpm/generated/data/images_xxx_cc_yy.tiff -n /path/to/ccpm/generated/data/isoVal_mapping.csv [--ard=15.,2,.5] -o /path/to/ccpm/generated/data/	
```

`ard` arguments stands for angle, radius and distance that are input arguments from the surface mesh extraction algorithms.
The run will produce 2 files : a STL/OFF of the surface mesh and a series of CSV files that can be displayed in Paraview 
through `TableToPoint` filter.
