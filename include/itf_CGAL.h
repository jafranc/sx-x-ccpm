/*
 * itf_CGAL.h
 *
 *  Created on: 20 déc. 2018
 *      Author: jfranc
 */

#ifndef ITF_CGAL_H_
#define ITF_CGAL_H_

#ifdef Success // avoid Eigen(via CGAL) X11 Success clash
#undef Success
#endif

#include <functional>

#include <CGAL/IO/STL.h>

/* CGAL header */
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/boost/graph/properties_Surface_mesh.h>

#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Surface_mesh_default_criteria_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Gray_level_image_3.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/Monge_via_jet_fitting.h>
#include <CGAL/wlop_simplify_and_regularize_point_set.h>
#include <CGAL/bilateral_smooth_point_set.h>
#include <CGAL/Polygon_mesh_slicer.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Polygon_mesh_processing/interpolated_corrected_curvatures.h>

//mesh refine and fair
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/refine.h>
#include <CGAL/Polygon_mesh_processing/fair.h>
//mesh smoothing
#include <CGAL/Polygon_mesh_processing/angle_and_area_smoothing.h>
#include <CGAL/Polygon_mesh_processing/detect_features.h>

//for least square plane fitting
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_Bezier_curve_traits_2.h>
#include <CGAL/Arrangement_2.h>

#include "../tpl/cgal/STL_Extension/include/CGAL/array.h"
#include "Poly_rings.h"

//#include <CGAL/natural_neighbor_coordinates_3.h>
//#include <CGAL/surface_neighbor_coordinates_3.h>
//#include <CGAL/interpolation_functions.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>

namespace PMP = CGAL::Polygon_mesh_processing;
//TODO refactor itf around Porjector Slicer I/O and Refiner classes
//TODO longer Bezier fitting impl
namespace ccpm {


/***
 * Interface template class which
 *  - input either an image (*.inr) or an suface (*.off)
 * 	- output image MCF cgal formated files (*.cgal)
 *
 * 	with possible access to intermediate (*.off) if needed
 */

    class itf_to_CGAL {

    public:

        //for c2t3
        typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
        typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;
        //for input
        typedef Tr::Geom_traits GT;
        typedef CGAL::Gray_level_image_3<GT::FT, GT::Point_3> Gray_level_image;
        typedef CGAL::Implicit_surface_3<GT, Gray_level_image> Surface_3;

        // for output (MCF)
        typedef CGAL::Simple_cartesian<double> Kernel;
//        typedef CGAL::Exact_predicates_inexact_constructions_kernel EpickKernel;
        typedef Kernel::Point_3 Point;
        typedef Kernel::Point_2 Point_2;
        typedef Kernel::Vector_3 Vector;
        typedef CGAL::Surface_mesh<Point> Triangle_mesh;
        //for interpolant
        typedef std::map<Point, double, Kernel::Less_xyz_3> CGAL_map;
//        typedef CGAL::Data_access<CGAL_map> CGAL_mapgetter;
        typedef std::vector<std::pair<Point, double> > PDVec;

        typedef boost::graph_traits<Triangle_mesh>::vertex_descriptor vertex_descriptor;
        typedef boost::graph_traits<Triangle_mesh>::face_descriptor face_descriptor;

        //curvature meas.
        typedef CGAL::Monge_via_jet_fitting<Kernel> Monge_via_jet_fitting;
        typedef Monge_via_jet_fitting::Monge_form Monge_form;
        typedef Triangle_mesh::Property_map<vertex_descriptor, int> Vertex_PM_type;
        typedef Triangle_mesh::Property_map <vertex_descriptor, std::vector<double>> Curvature_PM_type;

        //Gaussian, mean curvature another way
        typedef  Triangle_mesh::Property_map<vertex_descriptor, double > CurvatureMaps_type;

        //for local refinement
        //TODO look at what was there and what/why we extracted from
//    typedef T_PolyhedralSurf_rings<Triangle_mesh, Vertex_PM_type > Poly_rings
        typedef Poly_rings <Triangle_mesh, Vertex_PM_type> Poly_rings_def;


        //Kd Tree for shortest diatance
 class Image_point_property_map
{
  const std::vector<Point>& points;
public:
  typedef Point value_type;
  typedef const value_type& reference;
  typedef std::size_t key_type;
  typedef boost::lvalue_property_map_tag category;

  Image_point_property_map(const std::vector<Point>& pts):points(pts){}

  reference operator[](key_type k) const {return points[k];}

  friend reference get(const Image_point_property_map& ppmap,key_type i)
  {return ppmap[i];}
};

typedef CGAL::Search_traits_3<Kernel>                                        Traits_base;
typedef CGAL::Search_traits_adapter<std::size_t,Image_point_property_map,Traits_base> Traits;


typedef CGAL::Orthogonal_k_neighbor_search<Traits>                      K_neighbor_search;
typedef K_neighbor_search::Tree                                         Tree;
typedef Tree::Splitter                                                  Splitter;
typedef K_neighbor_search::Distance                                     Distance;



        struct refinement_parameters {

            bool isotropic = false;
            bool faired = true;
            // default parameter values and global variables
            unsigned int d_fitting = 2;
            unsigned int d_monge = 2;
            unsigned int nb_rings = 3;//seek min # of rings to get the required #pts
            unsigned int nb_points_to_use = 10;//
            bool verbose = (dbg_lvl > 2);
            unsigned int min_nb_points = (d_fitting + 1) * (d_fitting + 2) / 2;
        } ref_params;

        //for plane least suqare fitting
//        typedef Kernel::Line_2 Line;
        //Bezier
        typedef CGAL::CORE_algebraic_number_traits Nt_traits;
        typedef Nt_traits::Rational NT;
        typedef Nt_traits::Rational Rational;
        typedef Nt_traits::Algebraic Algebraic;
        typedef CGAL::Cartesian<Rational> Rat_kernel;
        typedef CGAL::Cartesian<Algebraic> Alg_kernel;
//        typedef Rat_kernel::Point_2 Rat_point;
        typedef CGAL::Arr_Bezier_curve_traits_2<Rat_kernel, Alg_kernel, Nt_traits>
                BTraits;
//        typedef Traits::X_monotone_curve_2 Bezier_x_monotone_curve;
//        typedef Traits::Curve_2 Bezier_curve;
        typedef CGAL::Arrangement_2<BTraits> Arrangement;

    private:
        static constexpr const double curvature_ratio = 0.1;
        static constexpr const double epsilon_box = 2.;
        static constexpr const int dbg_lvl = 0;

        CGAL_map m_F;


    public:

        itf_to_CGAL() : c2t3_(tr_), processed_sm_(false),
                        processed_refined_(false), do_refined_(false) {
            fname_ = "";
            //min_edge_length_ = 0.0;

            // surf mesh.
            // lower bound on minimum angle in degrees
            m_angle_ = 30.;
            // bound on radii of surface Delaunay balls
            m_radius_ = 0.5;
            // bound on Hausdorff distance does not play any role if bigger than
            // the square of the Uniform_size_criterion
            m_distance_ = 1.0;
        }

        //setters
        void set_input(const char *fname) {
            fname_ = fname;
            input_.read(fname_);
        }

        void set_ARD(float angle, float radius, float distance) {
            m_angle_ = angle;
            m_radius_ = radius;
            m_distance_ = distance;
            if (m_distance_ > m_radius_ * m_radius_)
                std::cout << "[Warning]  Distance has no impact if dist>radius^2 \n";
        }

        void set_isotropic_ref() {
            ref_params.isotropic = true;
        }

        void set_refined() {
            do_refined_ = true;
        }

        void set_faired() {
            ref_params.faired = true;
        }//TODO work on memory interface between CImg and CGAL

        const CGAL::Image_3 &get_input() {
            return input_;
        }

        const Triangle_mesh &get_tmesh() const {
            return tmesh_;
        }

        void input_neighFile(const std::string &fnamein) {
            //load csv
            std::ifstream csvfile(fnamein);
            int x, y, z, v;
            char delim;
            std::string header;
            getline(csvfile, header);
            while (csvfile.good()) {
                csvfile >> x >> delim >> y >> delim >> z >> delim >> v;
                m_F[Point{x, y, z}] = v;
            }

        }


        //TODO improve constness
        const C2t3 &get_surfmesh() {
            //do some processing
            process_mesh();
            return c2t3_;
        }

        const itf_to_CGAL &save_surf_off(const char *fname) {
            process_mesh();
            std::cerr << "\t [completed]\n";
            std::ofstream out(fname);
            CGAL::output_surface_facets_to_off(out, c2t3_);

            std::cout << "\n Final number of points: " << tr_.number_of_vertices() << "\n";

            return *this;
        }


        const itf_to_CGAL &save_surf_stl(const char *fname) {

            //prep for

            //prep for STL
            string v = "stl";
            auto sname = std::string(fname);
            if (sname.find(v) != std::string::npos) {
                size_t lastindex = sname.find_last_of(".");
                string prefix = sname.substr(0, lastindex);

                string postfix = sname.substr(lastindex + 1);
                string off_fname = prefix + (".off");
//                if (postfix == "inr") {//DO treat image to produce off
                save_surf_off(off_fname.c_str());
//                }
                read_off(off_fname.c_str());
                auto cpm = process_refined(false);
                //TODO refactor / clean up
//                read_off("refined.off");

                std::ofstream out(fname, std::ios::binary), csvfile(prefix + (".csv"));
                CGAL::IO::write_STL(out, tmesh_);


//                cpm = process_refined(false);

                PDVec coords;
                std::vector<Point> pts;
                for (const auto pd: m_F)
                    pts.push_back(pd.first);
                //for n-values
                Image_point_property_map ppmap(pts);
                Tree tree(boost::counting_iterator<std::size_t>(0),
                    boost::counting_iterator<std::size_t>(pts.size()),
                    Splitter(),
                    Traits(ppmap));
                Distance tr_dist(ppmap);

                auto [mcurv, gcurv] = process_curvature();

                csvfile << "x, y, z, k1, k2, km, kG, n, e" << std::endl;
                for (auto vd: boost::make_iterator_range(boost::vertices(tmesh_))) {
                    const auto pt = get(CGAL::vertex_point, tmesh_, vd);
                    const auto mf = get(cpm, vd);

                    //try #1
                    uintmax_t nv = tmesh_.num_vertices();
                    uintmax_t ne = tmesh_.num_edges();
                    uintmax_t nf = tmesh_.num_faces();
                    intmax_t euler = nv - ne + nf;
//try #2
//          LCC_3::Dart_const_descriptor root = m_lcc.dart_descriptor(CGAL::get_default_random().get_int(0,
//                                                            static_cast<int>(m_lcc.number_of_darts())));


                    K_neighbor_search search(tree, pt/*query*/,25/*nb nearest neigh*/ ,0,true,tr_dist);
                        csvfile << pt[0] << "," << pt[1] << "," << pt[2] << "," << mf[0] << "," << mf[1] << ","
                        << get(mcurv,vd) << "," << get(gcurv,vd) << ", ";

                        //get max over the 3-closest points
                    double avg = 0.;
                    for(K_neighbor_search::iterator it = search.begin(); it != search.end(); it++) {
                       avg = std::max(avg,m_F.at(pts[it->first]));
                    }

                    csvfile << avg << "," << euler << std::endl;

                }
                std::cout << "Final number of points: " << tr_.number_of_vertices() << "\n";
            } else
                throw std::exception();
            return *this;
        }

        const itf_to_CGAL &read_off(const char *fname) {
            std::ifstream input(fname);
            input >> tmesh_;

            std::cout << "Mesh Info : number of vertices " << tmesh_.number_of_vertices() << std::endl;
            std::cout << "Mesh Info : number of edges " << tmesh_.number_of_edges() << std::endl;
            std::cout << "Mesh Info : number of half-edges " << tmesh_.number_of_halfedges() << std::endl;
            std::cout << "Mesh Info : number of faces " << tmesh_.number_of_faces() << std::endl;

            return *this;
        }

    protected:
        ~itf_to_CGAL() {};


    private:
        Tr tr_; // intermediate result triangulation
        //if surface mesh deduced from image. C2t3 is intermediate data structure to Triangle_mesh
        C2t3 c2t3_;
        //this is the input for refinement, slicing, projection and MCF
        Triangle_mesh tmesh_;

        CGAL::Image_3 input_;
        std::vector<Kernel::Vector_3> normal;


    private:
        //double min_edge_length_;
        bool processed_sm_, processed_refined_;

        float m_angle_, m_radius_, m_distance_;
        bool do_refined_ = false;

        const char *fname_;

        void process() {
            process_mesh();
            if (do_refined_)
                process_refined(true);
        }

        //black box function
        void process_mesh() {
            if (!processed_sm_) {
                Gray_level_image tmp(input_, 2.9f);

                // the sphere actually bounds the whole image.
                GT::Point_3 bounding_sphere_center(input_.image()->xdim/2,input_.image()->ydim/2, input_.image()->zdim/2);
                auto max_dim = std::max( std::max(input_.image()->xdim, input_.image()->ydim) , input_.image()->zdim );
                GT::FT bounding_sphere_squared_radius = max_dim * max_dim  * 2.;
                GT::Sphere_3 bounding_sphere(bounding_sphere_center,
                                             bounding_sphere_squared_radius);

                // definition of the surface, with 10^-5 as relative precision
                Surface_3 surface(tmp, bounding_sphere, 1e-5);

                // defining meshing criteria   // -- control size mesh parameters
                CGAL::Surface_mesh_default_criteria_3<Tr> criteria(m_angle_,
                                                                   m_radius_,
                                                                   m_distance_);


                std::cout << "Using parameters:\n \t angle: "
                          << m_angle_ << "\n\t radius: "
                          << m_radius_ << "\n\t distance: "
                          << m_distance_ << std::endl;

                // meshing surface, with the "manifold without boundary" algorithm
                //the output mesh is guaranteed to be a manifold surface without boundary
                CGAL::make_surface_mesh(c2t3_, surface, criteria, CGAL::Manifold_tag(),50);
                //the output mesh is guaranteed to be manifold but may have boundaries.
                //CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());
                //the output mesh is guaranteed to be manifold but may have boundaries.
//                CGAL::make_surface_mesh(c2t3_, surface, criteria, CGAL::Manifold_with_boundary_tag(),350);



                ref_params.faired = true;
                if (ref_params.faired) {


                    std::vector<vertex_descriptor> in_points;
                    BOOST_FOREACH(vertex_descriptor vd, tmesh_.vertices()) {
                                    in_points.push_back(vd);
                                }

                    bool success = CGAL::Polygon_mesh_processing::fair(tmesh_, in_points);

                }

                processed_sm_ = true;
            }
        }

        std::pair<CurvatureMaps_type, CurvatureMaps_type> process_curvature()
        {
            CurvatureMaps_type mcurv, gcurv;
            bool created;

            boost::tie(mcurv, created) = tmesh_.add_property_map<vertex_descriptor, double >("v:mcm", 0);
            assert(created);
            boost::tie(gcurv, created) = tmesh_.add_property_map<vertex_descriptor, double >("v:gcm", 0);
            assert(created);

              PMP::interpolated_corrected_curvatures(tmesh_,
                CGAL::parameters::vertex_mean_curvature_map(mcurv)
                     .vertex_Gaussian_curvature_map(gcurv)
                     .ball_radius(0.5));

            return std::make_pair(mcurv,gcurv);
        }

        Curvature_PM_type process_refined(bool do_actually_refine) {

            //CUP
            Curvature_PM_type cpm;
            bool created;
            std::vector<double> delta;
            boost::tie(cpm, created) = tmesh_.add_property_map<vertex_descriptor, std::vector<double>>("v:cpm", {0, 0});
//            assert(created);

            if( true ) //ref_paramas.smoothing
            {
                typedef boost::property_map<Triangle_mesh, CGAL::edge_is_feature_t>::type EIFMap;
                EIFMap eif = get(CGAL::edge_is_feature, tmesh_);
                PMP::detect_sharp_edges(tmesh_, 60, eif);

                int sharp_counter = 0;
                for(auto e : edges(tmesh_))
                    if(get(eif, e))
                        ++sharp_counter;

                std::cout << "\n" << sharp_counter << " sharp edges" << std::endl;

                int nb_iterations = 30;
                PMP::angle_and_area_smoothing(tmesh_, CGAL::parameters::number_of_iterations(nb_iterations)
                        .use_safety_constraints(false) // authorize all moves
                        .edge_is_constrained_map(eif));


            if( true ) //ref_params.faired
            {
                std::vector<Point> in_points;  //container for data points
                BOOST_FOREACH(vertex_descriptor vd, tmesh_.vertices()) {
                                in_points.push_back(get(CGAL::vertex_point, tmesh_, vd));
                            }
                Kernel::Iso_cuboid_3 c3 = CGAL::bounding_box(in_points.begin(), in_points.end());
                in_points.clear();
                // second arg is epsilon detection
                std::set<face_descriptor> inside_faces = inside_face(c3, epsilon_box, tmesh_);


                std::vector<vertex_descriptor> point_from_faces;
//                BOOST_FOREACH(face_descriptor fd, inside_faces)
//                                BOOST_FOREACH(vertex_descriptor vd, vertices_around_face(tmesh_.halfedge(fd),
//                                                                                         tmesh_))
//                                                                                         point_from_faces.push_back(
//                                                        vd);

                for(auto ed : boost::make_iterator_range(edges(tmesh_)) ) {
                                if (get(eif, ed)) {
                                    point_from_faces.push_back(boost::source(ed,tmesh_));
                                    point_from_faces.push_back(boost::target(ed,tmesh_));
                                }
                            }

                bool success = CGAL::Polygon_mesh_processing::fair(tmesh_, point_from_faces);
                }
            }

            //DEBUG
            if (dbg_lvl > 1) {
                std::ofstream faired_off("faired.off");
                faired_off.precision(17);
                faired_off << tmesh_;
                faired_off.close();
            }

            if (!processed_refined_) {

                std::cout << " Analyzing refinement \n";

                std::set<face_descriptor> selected_faces;

                unsigned int nb_vertices_considered = 0;
                Vertex_PM_type vpm;
                boost::tie(vpm, created) = tmesh_.add_property_map<vertex_descriptor, int>("v:vpm", -1);
                assert(created);


                double mean_delta = 0.0;

                //filter inside BBox
                std::vector<Point> in_points;  //container for data points
                BOOST_FOREACH(vertex_descriptor vd, tmesh_.vertices()) {
                                in_points.push_back(get(CGAL::vertex_point, tmesh_, vd));
                            }
                Kernel::Iso_cuboid_3 c3 = CGAL::bounding_box(in_points.begin(), in_points.end());
                in_points.clear();
                // second arg is epsilon detection
                std::set<face_descriptor> inside_faces = inside_face(c3, epsilon_box, tmesh_);


                //gathering faces
                BOOST_FOREACH(face_descriptor f, inside_faces)
                                BOOST_FOREACH(vertex_descriptor vd, vertices_around_face(tmesh_.halfedge(f), tmesh_)) {
                                                in_points.clear();
                                                //gather points around the vertex using rings
                                                gather_fitting_points(tmesh_, vd, in_points, vpm);
                                                //skip if the nb of points is to small
                                                if (in_points.size() < ref_params.min_nb_points) {
                                                    std::cerr << "not enough pts for fitting this vertex"
                                                              << in_points.size() << std::endl;
                                                    vd++;
                                                }

                                                Monge_form monge_form;
                                                Monge_via_jet_fitting monge_fit;
                                                monge_form = monge_fit(in_points.begin(), in_points.end(),
                                                                       ref_params.d_fitting,
                                                                       ref_params.d_monge);
                                                //verbose txt output
                                                if (ref_params.verbose) {
                                                    std::vector<Point>::iterator itbp = in_points.begin(), itep = in_points.end();
//                                                    std::cout << "in_points list : " << std::endl;
//                                                    for (; itbp != itep; itbp++) std::cout << *itbp << std::endl;
                                                    std::cout << "--- vertex " << ++nb_vertices_considered
                                                              << " : " << get(CGAL::vertex_point, tmesh_, vd)
                                                              << std::endl
                                                              << "number of points used : " << in_points.size()
                                                              << std::endl
                                                              //<< monge_form << std::endl
                                                              << "coeff: " <<
                                                              monge_form.coefficients()[0] << " "
                                                              << monge_form.coefficients()[1]
                                                              << std::endl;
                                                }
                                                delta.push_back(
                                                        monge_form.coefficients()[0] * monge_form.coefficients()[0] +
                                                        monge_form.coefficients()[1] * monge_form.coefficients()[1]);

                                                cpm[vd] = monge_form.coefficients();//CUP


                                            }

                mean_delta = 1.0 / delta.size() * std::accumulate(delta.begin(), delta.end(), 0.0);
                std::cout << "refinement -- mean delta :" << mean_delta << std::endl;

                //now put face in bucket
                if (do_actually_refine) {
                    BOOST_FOREACH(face_descriptor f, inside_faces)
                                    BOOST_FOREACH(vertex_descriptor vd,
                                                  vertices_around_face(tmesh_.halfedge(f), tmesh_)) {
                                                    in_points.clear();
                                                    //gather points around the vertex using rings
                                                    gather_fitting_points(tmesh_, vd, in_points, vpm);
                                                    //skip if the nb of points is to small
                                                    if (in_points.size() < ref_params.min_nb_points) {
                                                        std::cerr
                                                                << "not enough pts for fitting this vertex"
                                                                << in_points.size()
                                                                << std::endl;
                                                        vd++;
                                                    }

                                                    Monge_form monge_form;
                                                    Monge_via_jet_fitting monge_fit;
                                                    monge_form = monge_fit(in_points.begin(), in_points.end(),
                                                                           ref_params.d_fitting, ref_params.d_monge);
                                                    //verbose txt output
                                                    if (ref_params.verbose) {
                                                        std::vector<Point>::iterator itbp = in_points.begin(), itep = in_points.end();
                                                        std::cout << "in_points list : " << std::endl;
//                                                    for (; itbp != itep; itbp++) std::cout << *itbp << std::endl;
                                                        std::cout << "--- vertex " << ++nb_vertices_considered
                                                                  << " : " << get(CGAL::vertex_point, tmesh_, vd)
                                                                  << std::endl
                                                                  << "number of points used : " << in_points.size()
                                                                  << std::endl
                                                                  << " monge cdt : " << monge_fit.condition_number()
                                                                  << std::endl
                                                                  //<< monge_form << std::endl
                                                                  << "coeff: " <<
                                                                  monge_form.coefficients()[0] << " "
                                                                  << monge_form.coefficients()[1]
                                                                  << std::endl;
                                                    }
                                                    if (monge_form.coefficients()[0] * monge_form.coefficients()[0] +
                                                        monge_form.coefficients()[1] * monge_form.coefficients()[1] >
                                                        curvature_ratio * mean_delta) {
                                                        selected_faces.insert(f);
                                                        BOOST_FOREACH(face_descriptor fi,
                                                                      faces_around_face(tmesh_.halfedge(f),
                                                                                        tmesh_)) selected_faces.insert(
                                                                                fi);
                                                    }
                                                }

                    //refining faces
                    if (ref_params.isotropic) {
                        PMP::isotropic_remeshing(
                                selected_faces,
                                1.,
                                tmesh_);
                    }
                    else
                    {
                        std::vector<face_descriptor> new_faces;
                        std::vector<vertex_descriptor> new_vertices;

                        PMP::refine(tmesh_,
                                    selected_faces,
                                    std::back_inserter(new_faces),
                                    std::back_inserter(new_vertices),
                                    CGAL::Polygon_mesh_processing::parameters::density_control_factor(5.));

                    }

                    //DEBUG
                    if (dbg_lvl > 1) {
                        std::ofstream refined_off("refined.off");
                        refined_off.precision(17);
                        refined_off << tmesh_;
                        refined_off.close();
                    }



                }// end of do_actually_refined

                processed_refined_ = true;
            }

            return cpm;

        }


        void gather_fitting_points(const Triangle_mesh &m, vertex_descriptor vd,
                                   std::vector<Point> &in_points,
                                   Vertex_PM_type &vpm) {
//container to collect vertices of v on the PolyhedralSurf
            std::vector<vertex_descriptor> gathered;
            //initialize
            in_points.clear();
            //OPTION -p nb_points_to_use, with nb_points_to_use != 0. Collect
            //enough rings and discard some points of the last collected ring to
            //get the exact "nb_points_to_use"
            if (ref_params.nb_points_to_use != 0) {
                Poly_rings_def::collect_enough_rings(m, vd, ref_params.nb_points_to_use, gathered, vpm);
                if (gathered.size() > ref_params.nb_points_to_use) gathered.resize(ref_params.nb_points_to_use);
            } else { // nb_points_to_use=0, this is the default and the option -p is not considered;
                // then option -a nb_rings is checked. If nb_rings=0, collect
                // enough rings to get the min_nb_points required for the fitting
                // else collect the nb_rings required
                if (ref_params.nb_rings == 0)
                    Poly_rings_def::collect_enough_rings(m, vd, ref_params.min_nb_points, gathered, vpm);
                else Poly_rings_def::collect_i_rings(m, vd, ref_params.nb_rings, gathered, vpm);
            }
            //store the gathered points
            BOOST_FOREACH(vertex_descriptor vdi, gathered)in_points.push_back(get(CGAL::vertex_point, m, vdi));

        }

        std::set<face_descriptor> inside_face(const Kernel::Iso_cuboid_3 &c, double eps, const Triangle_mesh &m) {

            std::set<face_descriptor> selected;
            double xmin, xmax, ymin, ymax, zmin, zmax;
            xmin = c.xmin() + eps;
            xmax = c.xmax() - eps;
            ymin = c.ymin() + eps;
            ymax = c.ymax() - eps;
            zmin = c.zmin() + eps;
            zmax = c.zmax() - eps;

            std::cout << " Bounding box " << c << std::endl;

            BOOST_FOREACH(face_descriptor fd, m.faces()) {
                            BOOST_FOREACH(vertex_descriptor vd, vertices_around_face(m.halfedge(fd), m)) {
                                            Point p = get(CGAL::vertex_point, m, vd);
                                            if (p[0] > xmin && p[0] < xmax &&
                                                p[1] > ymin && p[1] < ymax &&
                                                p[2] > zmin && p[2] < zmax)
                                                selected.insert(fd);
                                        }
                        }

            return selected;

        }

    public:
        //conversion function
        const vector <Kernel::Vector_3> &get_normal() const {
            return normal;
        }
    };

}//end namespace mema

#endif /* ITF_CGAL_H_ */
