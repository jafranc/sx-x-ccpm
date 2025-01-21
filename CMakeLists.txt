cmake_minimum_required (VERSION 2.8.12)
project(ccpm)

set(exe "ccpm")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/build)

set(include "./include")
set(ccpm_include
        include/ )

include_directories(${include})

set( ccpm_src
    src/main.cpp )

find_package(OpenMP REQUIRED)

find_package( Threads REQUIRED )
# CImg for stl import
set(CIMG_DIR ./tpl/CImg)
include_directories(${CIMG_DIR})
add_compile_options("-Dcimg_display=0")

find_package( TIFF REQUIRED )
if(TIFF_FOUND)
  find_package(TIFF)
  include_directories(${TIFF_INCLUDE_DIR})
else()
  message(STATUS "These programs require the TIFF library, and will not be compiled.")
  #   message(STATUS "These programs require the X11 library, and will not be compiled.")
endif()


#cxxopts parser
set(CXXOPTS_DIR ./tpl/cxxopts/include)
include_directories(${CXXOPTS_DIR})
find_package(OpenMP)


# cgal
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
if(Eigen3_FOUND)
  include_directories(${EIGEN3_INCLUDE_DIRS})
endif()

# CGAL for all computational geom
set(CGAL_DIR ./tpl/cgal)
find_package(CGAL REQUIRED) #OPTIONAL_COMPONENTS Core ImageIO)
if(CGAL_FOUND)
  set(CGAL_ALL_INCLUDES
          ${CGAL_DIR}/SearchStructures/include
          ${CGAL_DIR}/Heat_method_3/include
          ${CGAL_DIR}/Ridges_3/include
          ${CGAL_DIR}/Periodic_3_triangulation_3/test/Periodic_3_triangulation_3/include
          ${CGAL_DIR}/Periodic_3_triangulation_3/include
          ${CGAL_DIR}/Jet_fitting_3/include
          ${CGAL_DIR}/Surface_mesh_shortest_path/test/Surface_mesh_shortest_path/include
          ${CGAL_DIR}/Surface_mesh_shortest_path/include
          ${CGAL_DIR}/HalfedgeDS/include
          ${CGAL_DIR}/Arrangement_on_surface_2/include
          ${CGAL_DIR}/Surface_mesh_topology/include
          ${CGAL_DIR}/Classification/include
          ${CGAL_DIR}/OpenNL/include
          ${CGAL_DIR}/Hyperbolic_triangulation_2/include
          ${CGAL_DIR}/Hyperbolic_triangulation_2/demo/Hyperbolic_triangulation_2/include
          ${CGAL_DIR}/Point_set_3/include
          ${CGAL_DIR}/Polygonal_surface_reconstruction/include
          ${CGAL_DIR}/Box_intersection_d/include
          ${CGAL_DIR}/Surface_mesh_deformation/include
          ${CGAL_DIR}/Testsuite/include
          ${CGAL_DIR}/Circular_kernel_3/test/Circular_kernel_3/include
          ${CGAL_DIR}/Circular_kernel_3/include
          ${CGAL_DIR}/Nef_2/test/Nef_2/include
          ${CGAL_DIR}/Nef_2/include
          ${CGAL_DIR}/Distance_3/include
          ${CGAL_DIR}/Surface_mesh_approximation/include
          ${CGAL_DIR}/Arithmetic_kernel/test/Arithmetic_kernel/include
          ${CGAL_DIR}/Arithmetic_kernel/include
          ${CGAL_DIR}/Convex_hull_d/test/Convex_hull_d/include
          ${CGAL_DIR}/Convex_hull_d/include
          ${CGAL_DIR}/Shape_detection/test/Shape_detection/include
          ${CGAL_DIR}/Shape_detection/include
          ${CGAL_DIR}/Convex_hull_3/test/Convex_hull_3/include
          ${CGAL_DIR}/Convex_hull_3/include
          ${CGAL_DIR}/LEDA/include
          ${CGAL_DIR}/Kernel_23/test/Kernel_23/include
          ${CGAL_DIR}/Kernel_23/include
          ${CGAL_DIR}/Surface_mesh_parameterization/include
          ${CGAL_DIR}/Minkowski_sum_3/include
          ${CGAL_DIR}/Optimal_transportation_reconstruction_2/include
          ${CGAL_DIR}/Visibility_2/test/Visibility_2/include
          ${CGAL_DIR}/Visibility_2/include
          ${CGAL_DIR}/Optimisation_basic/include
          ${CGAL_DIR}/Inscribed_areas/include
          ${CGAL_DIR}/Minkowski_sum_2/include
          ${CGAL_DIR}/Alpha_wrap_3/include
          ${CGAL_DIR}/Homogeneous_kernel/include
          ${CGAL_DIR}/Three/include
          ${CGAL_DIR}/Property_map/include
          ${CGAL_DIR}/Interval_support/include
          ${CGAL_DIR}/Nef_3/test/Nef_3/include
          ${CGAL_DIR}/Nef_3/include
          ${CGAL_DIR}/Random_numbers/include
          ${CGAL_DIR}/Periodic_4_hyperbolic_triangulation_2/include
          ${CGAL_DIR}/Periodic_4_hyperbolic_triangulation_2/demo/Periodic_4_hyperbolic_triangulation_2/include
          ${CGAL_DIR}/Envelope_2/include
          ${CGAL_DIR}/Surface_mesh_simplification/include
          ${CGAL_DIR}/Union_find/test/Union_find/include
          ${CGAL_DIR}/Union_find/include
          ${CGAL_DIR}/CGAL_Core/include
          ${CGAL_DIR}/Stream_lines_2/include
          ${CGAL_DIR}/Surface_sweep_2/include
          ${CGAL_DIR}/Circulator/include
          ${CGAL_DIR}/Point_set_2/include
          ${CGAL_DIR}/Orthtree/include
          ${CGAL_DIR}/Periodic_3_mesh_3/include
          ${CGAL_DIR}/Polyhedron/include
          ${CGAL_DIR}/Polyhedron/demo/Polyhedron/include
          ${CGAL_DIR}/Polyhedron/demo/Polyhedron/Plugins/Mesh_3/include
          ${CGAL_DIR}/Interpolation/test/Interpolation/include
          ${CGAL_DIR}/Interpolation/include
          ${CGAL_DIR}/Segment_Delaunay_graph_2/test/Segment_Delaunay_graph_2/include
          ${CGAL_DIR}/Segment_Delaunay_graph_2/include
          ${CGAL_DIR}/Spatial_searching/include
          ${CGAL_DIR}/Spatial_searching/benchmark/Spatial_searching/include
          ${CGAL_DIR}/CGAL_ipelets/include
          ${CGAL_DIR}/CGAL_ipelets/demo/CGAL_ipelets/include
          ${CGAL_DIR}/Hash_map/test/Hash_map/include
          ${CGAL_DIR}/Hash_map/include
          ${CGAL_DIR}/Surface_mesh/include
          ${CGAL_DIR}/STL_Extension/include
          ${CGAL_DIR}/Intersections_2/include
          ${CGAL_DIR}/Triangulation/test/Triangulation/include
          ${CGAL_DIR}/Triangulation/include
          ${CGAL_DIR}/Combinatorial_map/include
          ${CGAL_DIR}/Cone_spanners_2/include
          ${CGAL_DIR}/Profiling_tools/include
          ${CGAL_DIR}/Algebraic_kernel_for_spheres/test/Algebraic_kernel_for_spheres/include
          ${CGAL_DIR}/Algebraic_kernel_for_spheres/include
          ${CGAL_DIR}/GraphicsView/include
          ${CGAL_DIR}/GraphicsView/demo/L1_Voronoi_diagram_2/include
          ${CGAL_DIR}/GraphicsView/demo/Segment_Delaunay_graph_2/include
          ${CGAL_DIR}/GraphicsView/demo/Periodic_2_triangulation_2/include
          ${CGAL_DIR}/GraphicsView/demo/Segment_Delaunay_graph_Linf_2/include
          ${CGAL_DIR}/GraphicsView/demo/Triangulation_2/include
          ${CGAL_DIR}/Surface_mesher/include
          ${CGAL_DIR}/Polynomial/test/Polynomial/include
          ${CGAL_DIR}/Polynomial/include
          ${CGAL_DIR}/Modifier/include
          ${CGAL_DIR}/AABB_tree/include
          ${CGAL_DIR}/TDS_2/test/TDS_2/include
          ${CGAL_DIR}/TDS_2/include
          ${CGAL_DIR}/Periodic_2_triangulation_2/test/Periodic_2_triangulation_2/include
          ${CGAL_DIR}/Periodic_2_triangulation_2/include
          ${CGAL_DIR}/Stream_support/include
          ${CGAL_DIR}/Cartesian_kernel/include
          ${CGAL_DIR}/TDS_3/test/TDS_3/include
          ${CGAL_DIR}/TDS_3/include
          ${CGAL_DIR}/Polyline_simplification_2/include
          ${CGAL_DIR}/Polyline_simplification_2/demo/Polyline_simplification_2/include
          ${CGAL_DIR}/Surface_mesh_segmentation/include
          ${CGAL_DIR}/Mesh_3/include
          ${CGAL_DIR}/Generator/include
          ${CGAL_DIR}/Boolean_set_operations_2/include
          ${CGAL_DIR}/Skin_surface_3/include
          ${CGAL_DIR}/Straight_skeleton_2/test/Straight_skeleton_2/include
          ${CGAL_DIR}/Straight_skeleton_2/include
          ${CGAL_DIR}/Nef_S2/test/Nef_S2/include
          ${CGAL_DIR}/Nef_S2/include
          ${CGAL_DIR}/Installation/include
          ${CGAL_DIR}/Voronoi_diagram_2/test/Voronoi_diagram_2/include
          ${CGAL_DIR}/Voronoi_diagram_2/include
          ${CGAL_DIR}/Convex_hull_2/test/Convex_hull_2/include
          ${CGAL_DIR}/Convex_hull_2/include
          ${CGAL_DIR}/Triangulation_3/test/Triangulation_3/include
          ${CGAL_DIR}/Triangulation_3/include
          ${CGAL_DIR}/Kernel_d/test/Kernel_d/include
          ${CGAL_DIR}/Kernel_d/include
          ${CGAL_DIR}/Intersections_3/include
          ${CGAL_DIR}/Barycentric_coordinates_2/include
          ${CGAL_DIR}/Point_set_processing_3/include
          ${CGAL_DIR}/Convex_decomposition_3/include
          ${CGAL_DIR}/Principal_component_analysis_LGPL/include
          ${CGAL_DIR}/Interval_skip_list/include
          ${CGAL_DIR}/SMDS_3/include
          ${CGAL_DIR}/Tetrahedral_remeshing/include
          ${CGAL_DIR}/Snap_rounding_2/include
          ${CGAL_DIR}/NewKernel_d/include
          ${CGAL_DIR}/Set_movable_separability_2/include
          ${CGAL_DIR}/Segment_Delaunay_graph_Linf_2/test/Segment_Delaunay_graph_Linf_2/include
          ${CGAL_DIR}/Segment_Delaunay_graph_Linf_2/include
          ${CGAL_DIR}/Mesher_level/include
          ${CGAL_DIR}/Partition_2/include
          ${CGAL_DIR}/QP_solver/include
          ${CGAL_DIR}/Filtered_kernel/include
          ${CGAL_DIR}/Triangulation_2/test/Triangulation_2/include
          ${CGAL_DIR}/Triangulation_2/include
          ${CGAL_DIR}/Modular_arithmetic/include
          ${CGAL_DIR}/Surface_mesh_skeletonization/include
          ${CGAL_DIR}/Polytope_distance_d/include
          ${CGAL_DIR}/Algebraic_foundations/include
          ${CGAL_DIR}/Optimal_bounding_box/include
          ${CGAL_DIR}/Advancing_front_surface_reconstruction/include
          ${CGAL_DIR}/Generalized_map/include
          ${CGAL_DIR}/Principal_component_analysis/include
          ${CGAL_DIR}/Bounding_volumes/include
          ${CGAL_DIR}/Polygon/include
          ${CGAL_DIR}/Solver_interface/include
          ${CGAL_DIR}/Alpha_shapes_3/test/Alpha_shapes_3/include
          ${CGAL_DIR}/Alpha_shapes_3/include
          ${CGAL_DIR}/Spatial_sorting/include
          ${CGAL_DIR}/Algebraic_kernel_for_circles/test/Algebraic_kernel_for_circles/include
          ${CGAL_DIR}/Algebraic_kernel_for_circles/include
          ${CGAL_DIR}/Envelope_3/include
          ${CGAL_DIR}/Triangulation_on_sphere_2/include
          ${CGAL_DIR}/Algebraic_kernel_d/test/Algebraic_kernel_d/include
          ${CGAL_DIR}/Algebraic_kernel_d/include
          ${CGAL_DIR}/Linear_cell_complex/include
          ${CGAL_DIR}/Alpha_shapes_2/include
          ${CGAL_DIR}/BGL/include
          ${CGAL_DIR}/Number_types/test/Number_types/include
          ${CGAL_DIR}/Number_types/include
          ${CGAL_DIR}/Straight_skeleton_extrusion_2/include
          ${CGAL_DIR}/Poisson_surface_reconstruction_3/include
          ${CGAL_DIR}/Shape_regularization/test/Shape_regularization/include
          ${CGAL_DIR}/Shape_regularization/include
          ${CGAL_DIR}/Matrix_search/include
          ${CGAL_DIR}/Circular_kernel_2/test/Circular_kernel_2/include
          ${CGAL_DIR}/Circular_kernel_2/include
          ${CGAL_DIR}/Scale_space_reconstruction_3/include
          ${CGAL_DIR}/Polygon_mesh_processing/include
          ${CGAL_DIR}/Subdivision_method_3/include
          ${CGAL_DIR}/Distance_2/include
          ${CGAL_DIR}/Mesh_2/include
          ${CGAL_DIR}/Apollonius_graph_2/test/Apollonius_graph_2/include
          ${CGAL_DIR}/Apollonius_graph_2/include
          ${CGAL_DIR}/Weights/test/Weights/include
          ${CGAL_DIR}/Weights/include
          ${CGAL_DIR}/CGAL_ImageIO/include
          )

  include_directories(${CGAL_ALL_INCLUDES})
endif ()

find_package(Boost COMPONENTS filesystem REQUIRED)

add_executable(${exe} ${ccpm_src})
target_link_libraries(${exe}  ${TIFF_LIBRARIES} Eigen3::Eigen OpenMP::OpenMP_CXX Threads::Threads CGAL::CGAL ${Boost_FILESYSTEM_LIBRARY})

#tests
#add_subdirectory (tests)