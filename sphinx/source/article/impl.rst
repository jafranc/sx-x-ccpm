.. _impl:

How is it implemented
===============================

Implementation is leveraging CIMG's connected component and
CGAL's Polygon Mesh Processing (`PMP <https://doc.cgal.org/latest/Polygon_mesh_processing/index.html>`_) as well as CGAL's Surface Meshing (`SM <https://doc.cgal.org/latest/Surface_mesh/index.html>`_).

From image to 3D surface meshing
---------------------------------

From a segmented images, in which we have good confidence in the phase labeling. There exist :math:`n` liquid or gas phase
and let's assume :math:`1` solid phase. Depending on the main physical process (imbibition or drainage), one of these liquid or gas phase
is more connected than the other phases. Solid and the most connected phases can then be considered closure and only connected components of the
:math:`n-1` remaining phases can be processed.

This first step is then leveraging CIMG to read, isolate phase from the other phases and eventually each individual (large enough) connected components.connected.
Let us refered to them as :math:`B(i,c)` for phase :math:`i` anc component :math:`c`.
It uses `label() <https://cimg.eu/reference/structcimg__library_1_1CImg.html#a60b19453f1d63efd59dd99921cd6b513>`_ and io capability to produce an _.inr_(inria image format) output.

Surface Mesh generation
-------------------------

Then from each individual binarized connected components :math:`B(i,c)`, a 3D surface mesh is built, placing 50 initial points on the intersection of the surface with random
rays. This point cloud serve as initial stage for 3D surface triangulation algorith which is parameterizable with Angle-Radius-Distance (ARD) parameters.
They stand for the minimum angle in degree, the bound on radii of Delaunay balls and the bound on Hausdorff distance. Default is 30, 0.5 and 1 respectively.

.. note::

    The D (distance) is Hausdorff distance for inscribed triangle, it has no effect if it is larger than the square of radius.

.. note::

    The random process of placing initial points implies that two consecutive runs on the same initial image will results in different results.


Mesh quality enhancement and Gaussian curvature measurement
---------------------------------------------------------------------------
The method is leveraging preprocessing treatment from CGAL to enhance mesh quality around the area of interest. To do that, sharp edges are detected as exceeding 60 degree between two triangles.
Those points and edges are smoothed so that they share comparable values. Eventually Gaussian curvature is evaluated from interpolated local corrected curvatures.

