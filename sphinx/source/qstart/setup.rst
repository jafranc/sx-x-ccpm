.. _setup:

How to install CCPM
=======================

Though this program should work on *Windows* based and *MacOs* based distro,
it is not tested in the CI/CD process nor developed on such distro and hence
considered not supported.

Linux-based
------------

These external dependencies are required for CCPM to work properly:

1. CGAL (as submodule)
2. CImg (as submodule)
3. Loguru (as submodule)
4. cxxopts (as submodule)

`CImg <https://cimg.eu/>`_ is used for processing binarized image, isolating connected component of gas labeled phase.
It is a well-known versatile header library for image processing. It offers simplicity and efficiency.

`CGAL <https://www.cgal.org/>`_ is the inria computational geometry toolbox, that is a reference in the meshing community.
Highly templatized, it offers state of the art remeshing and computational geometry implementation.

`loguru <https://github.com/emilk/loguru>`_ and `cxxopts <https://github.com/jarro2783/cxxopts>`_ are two handy library that are used for reasonable output and logging on one hand and
command line argument parsing and processing on the other hand.

Note that `ceres <http://ceres-solver.org/index.html>`_ non-linear optimization toolbox and `libtiff <https://libtiff.gitlab.io/libtiff/>`_ are required as dependencies for mesh processing and io operation. On debian distro, they are easily installable through,

.. code:: bash

    sudo apt install libtiff5-dev

`libboost-filesystem` is also used in order to make output directory and pathing. It can be install along all dependencies with

.. code:: bash

    sudo apt install libboost-all-dev

if not already present.

Code compilation
-----------------

Once the dependencies satisfied, a straight-forward `cmake <https://github.com/Kitware/CMake>`_ out-of-source install can be done

.. code:: bash

    mkdir build/ && cd build/
    cmake -DCMAKE_BUILD_TYPE=Release ..
    make -j 4


Code launch
------------

Once compiled the code is launch in several stages. Refer to command :ref:`command`
