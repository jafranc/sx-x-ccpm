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
#include <boost/filesystem.hpp>

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

#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>

//CC
#include <boost/graph/connected_components.hpp>
#include <boost/graph/incremental_components.hpp>
#include <boost/filesystem/path.hpp>

#include "../tpl/loguru/loguru.hpp"

namespace PMP = CGAL::Polygon_mesh_processing;

namespace ccpm
{

/***
 * Interface template class which
 *  - input either an image (*.inr) or an suface (*.off)
 *  - output image MCF cgal formated files (*.cgal)
 *
 *  with possible access to intermediate (*.off) if needed
 */

class itf_to_CGAL
{

public:

  //for c2t3
  typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
  typedef CGAL::Complex_2_in_triangulation_3< Tr > C2t3;
  //for input
  typedef Tr::Geom_traits GT;
  typedef CGAL::Gray_level_image_3< GT::FT, GT::Point_3 > Gray_level_image;
  typedef CGAL::Implicit_surface_3< GT, Gray_level_image > Surface_3;

  // for output (MCF)
  typedef CGAL::Simple_cartesian< double > Kernel;
  typedef Kernel::Point_3 Point;
  typedef Kernel::Point_2 Point_2;
  typedef Kernel::Vector_3 Vector;
  typedef CGAL::Surface_mesh< Point > Triangle_mesh;
  //for interpolant
  typedef std::map< Point, double, Kernel::Less_xyz_3 > CGAL_map;
  typedef std::vector< std::pair< Point, double > > PDVec;

  typedef boost::graph_traits< Triangle_mesh >::vertex_descriptor vertex_descriptor;
  typedef boost::graph_traits< Triangle_mesh >::vertices_size_type vertex_index;
  typedef boost::graph_traits< Triangle_mesh >::edge_descriptor edge_descriptor;
  typedef boost::graph_traits< Triangle_mesh >::face_descriptor face_descriptor;

  //curvature meas.
  typedef Triangle_mesh::Property_map< vertex_descriptor, int > Vertex_PM_type;
  typedef Triangle_mesh::Property_map< vertex_descriptor, std::vector< double > > Curvature_PM_type;

  //Gaussian, mean curvature another way
  typedef Triangle_mesh::Property_map< vertex_descriptor, double > CurvatureMaps_type;

  //Kd Tree for shortest diatance
  class Image_point_property_map
  {
public:
    const std::vector< Point > & points;

    typedef Point value_type;
    typedef const value_type & reference;
    typedef std::size_t key_type;
    typedef boost::lvalue_property_map_tag category;

    Image_point_property_map(): points( std::vector< Point >()) {}

    Image_point_property_map( const std::vector< Point > & pts ): points( pts ) {}

    reference operator[]( key_type k ) const { return points[k]; }

    friend reference get( const Image_point_property_map & ppmap, key_type i ) { return ppmap[i]; }
  };

  typedef CGAL::Search_traits_3< Kernel > Traits_base;
  typedef CGAL::Search_traits_adapter< std::size_t, Image_point_property_map, Traits_base > Traits;


  typedef CGAL::Orthogonal_k_neighbor_search< Traits > K_neighbor_search;
  typedef K_neighbor_search::Tree Tree;
  typedef Tree::Splitter Splitter;
  typedef K_neighbor_search::Distance Distance;

  //for plane least suqare fitting
  //Bezier
  typedef CGAL::CORE_algebraic_number_traits Nt_traits;
  typedef Nt_traits::Rational NT;
  typedef Nt_traits::Rational Rational;
  typedef Nt_traits::Algebraic Algebraic;
  typedef CGAL::Cartesian< Rational > Rat_kernel;
  typedef CGAL::Cartesian< Algebraic > Alg_kernel;
  typedef CGAL::Arr_Bezier_curve_traits_2< Rat_kernel, Alg_kernel, Nt_traits >
    BTraits;
  typedef CGAL::Arrangement_2< BTraits > Arrangement;

private:
  double TRIPLE_VALUE;

  CGAL_map m_F;
  Triangle_mesh::Property_map< vertex_descriptor, double > triple_map_;


public:

  itf_to_CGAL(): c2t3_( tr_ ), processed_sm_( false ), do_refined_( false )
  {
    fname_ = "";
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
  int set_input( const char *fname )
  {
    fname_ = fname;
    return (input_.read( fname_ )) ? 0 : 1;
  }

  void set_cutoff( double cval )
  {
    TRIPLE_VALUE = cval;
  }

  int set_ARD( float angle, float radius, float distance )
  {
    m_angle_ = angle;
    m_radius_ = radius;
    m_distance_ = distance;
    if( m_distance_ > m_radius_ * m_radius_ )
      LOG_S( ERROR ) << "[Warning]  Distance has no impact if dist>radius^2 \n";

    return 0;
  }

  void set_refined()
  {
    do_refined_ = true;
  }

  const CGAL::Image_3 & get_input()
  {
    if( input_.image() == 0 )
      LOG_S( ERROR ) << " Input undefined or ill-defined";
    return input_;
  }

  int input_neighFile( const std::string & fnamein )
  {
    //load csv
    std::ifstream csvfile( fnamein );
    int x, y, z, v;
    char delim;
    std::string header;
    getline( csvfile, header );
    while( csvfile.good())
    {
      csvfile >> x >> delim >> y >> delim >> z >> delim >> v;
      m_F[Point{x, y, z}] = v;
    }

    return 0;

  }

  const itf_to_CGAL & save_surf_off( const char *fname )
  {
    process_mesh();
    LOG_S( INFO ) << "\t [completed]";
    std::ofstream out( fname );
    CGAL::output_surface_facets_to_off( out, c2t3_ );

    LOG_S( INFO ) << " Final number of points: " << tr_.number_of_vertices();

    return *this;
  }


  int save_surf_stl( const char *fname, std::vector< double > & angles )
  {

    //prep for
    //prep for STL
    std::string v = "stl";
    auto fullpath = boost::filesystem::path( fname );
    auto sname = fullpath.leaf().string();
    if( sname.find( v ) != std::string::npos )
    {
      size_t lastindex = sname.find_last_of( "." );
      std::string prefix = sname.substr( 0, lastindex );

      std::string postfix = sname.substr( lastindex + 1 );
      std::string off_fname = prefix + (".off");
      save_surf_off( off_fname.c_str());
      read_off( off_fname.c_str());
      auto cpm = process_refined( false );
      //TODO refactor / clean up

      std::ofstream out( (fullpath.branch_path() / (prefix + (".stl"))).c_str(), std::ios::binary ), csvfile( (fullpath.branch_path() / ( prefix + (".csv")).c_str()) );
      LOG_S( INFO ) << "writing " <<(fullpath.branch_path() / (prefix + (".stl")) );
      LOG_S( INFO ) << "writing " <<(fullpath.branch_path() / (prefix + (".csv")) );
      CGAL::IO::write_STL( out, tmesh_ );

      //for n-values
      auto [mcurv, gcurv] = process_curvature();

      PDVec coords;
      std::vector< Point > pts;
      for( const auto pd: m_F )
        pts.push_back( pd.first );
      //produce angles
      build_tree( pts );
      double integral_constant = get_totalIntegral( gcurv );
      LOG_S( INFO ) << "integral_constant :: " << integral_constant;

      LOG_S( INFO ) << " angle :: label, numel, euler \n [iv, area]";
      auto vvec = process_triple_contact();

      //here is the loop on contact liens
      LOG_S( INFO ) << "Start processing CC";
      for( const auto &[k, v]: vvec.first )
      {
        const auto & map = vvec.second[k];           //is the set of faces
        auto euler = v;
        if( map.size() > 10 )
        {
          auto basicStringName = std::string( prefix + "-csvfile_" + std::to_string( k ) + ".csv" );
          LOG_S( INFO ) << "writing " << basicStringName;
          std::ofstream basicOfstream(
            basicStringName.c_str());
          basicOfstream << " # x, y ,z, angle \n";
          angles.push_back( get_angle( map, gcurv, [integral_constant]( double t ) {
            return acos((integral_constant - t) / 2. / M_PI - 1. ) * 180 / M_PI;
          }, basicOfstream ));
          LOG_S( INFO ) << k << ", " << map.size()
                        << "," << euler << ", "
                        << angles.back();
        }
        else
          LOG_S( INFO ) << "Discard CC # " << k;

      }


      csvfile << "x, y, z, km, kG, n, e" << std::endl;
      for( auto vd: boost::make_iterator_range( boost::vertices( tmesh_ )))
      {
        const auto pt = get( CGAL::vertex_point, tmesh_, vd );
        const auto mf = get( cpm, vd );

        //try #1
        uintmax_t nv = tmesh_.num_vertices();
        uintmax_t ne = tmesh_.num_edges();
        uintmax_t nf = tmesh_.num_faces();
        intmax_t euler = nv - ne + nf;

        csvfile << pt[0] << "," << pt[1] << "," << pt[2] << "," << mf[0] << "," << mf[1] << ","
                << get( mcurv, vd ) << "," << get( gcurv, vd ) << ", ";
        csvfile << tags_projection( pt ) << "," << euler << std::endl;

      }
      LOG_S( INFO ) << "Final number of points: " << tr_.number_of_vertices();
    }
    else
      throw std::exception();

    return 0;
  }

  const itf_to_CGAL & read_off( const char *fname )
  {
    std::ifstream input( fname );
    input >> tmesh_;

    LOG_S( INFO ) << "Mesh Info : number of vertices " << tmesh_.number_of_vertices();
    LOG_S( INFO ) << "Mesh Info : number of edges " << tmesh_.number_of_edges();
    LOG_S( INFO ) << "Mesh Info : number of half-edges " << tmesh_.number_of_halfedges();
    LOG_S( INFO ) << "Mesh Info : number of faces " << tmesh_.number_of_faces();

    return *this;
  }

protected:
  ~itf_to_CGAL() {};


private:
  Tr tr_;       // intermediate result triangulation
  //if surface mesh deduced from image. C2t3 is intermediate data structure to Triangle_mesh
  C2t3 c2t3_;
  //this is the input for refinement, slicing, projection and MCF
  Triangle_mesh tmesh_;

  CGAL::Image_3 input_;
  std::vector< Kernel::Vector_3 > normal;

  //double interpolation
  std::unique_ptr< Tree > ptree_;
  std::unique_ptr< Distance > p_tr_dist_;

  typedef boost::component_index< int > Components;
  //double min_edge_length_;
  bool processed_sm_;

  float m_angle_, m_radius_, m_distance_;
  bool do_refined_ = false;

  const char *fname_;

  void process()
  {
    process_mesh();
    if( do_refined_ )
      process_refined( true );
  }

  //black box function
  void process_mesh()
  {
    if( !processed_sm_ )
    {
      const auto & input = get_input();
      Gray_level_image tmp( input, 2.9f );

      // the sphere actually bounds the whole image.
      GT::Point_3 bounding_sphere_center( input.image()->xdim / 2, input.image()->ydim / 2,
                                          input.image()->zdim / 2 );
      auto max_dim = std::max( std::max( input.image()->xdim, input.image()->ydim ), input.image()->zdim );
      GT::FT bounding_sphere_squared_radius = max_dim * max_dim * 2.;
      GT::Sphere_3 bounding_sphere( bounding_sphere_center,
                                    bounding_sphere_squared_radius );

      // definition of the surface, with 10^-5 as relative precision
      Surface_3 surface( tmp, bounding_sphere, 1e-5 );

      // defining meshing criteria   // -- control size mesh parameters
      CGAL::Surface_mesh_default_criteria_3< Tr > criteria( m_angle_,
                                                            m_radius_,
                                                            m_distance_ );


      LOG_S( INFO ) << "Using parameters:\n \t angle: "
                    << m_angle_ << "\n\t radius: "
                    << m_radius_ << "\n\t distance: "
                    << m_distance_;

      // meshing surface, with the "manifold without boundary" algorithm
      //the output mesh is guaranteed to be a manifold surface without boundary
      CGAL::make_surface_mesh( c2t3_, surface, criteria, CGAL::Manifold_tag(), 50 );
      processed_sm_ = true;
    }
  }

  std::pair< std::unordered_map< int, int >,
             std::unordered_map< int, std::set< face_descriptor > > > process_triple_contact()
  {
    //TODO build a connected component list of faces whose is in contact with one vertices of triple line
    build_triple_map();
    std::unordered_map< int, int > eulers;
    std::unordered_map< int, std::set< face_descriptor > > vvec;
    std::vector< vertex_descriptor > cc;

    //try #2 -- build filtered graph on triple value
    Is_triple filter( triple_map_, tmesh_, TRIPLE_VALUE );
    boost::filtered_graph< Triangle_mesh, Is_triple, Is_triple > filtered_mesh( tmesh_, filter, filter );
    // compute cc
    std::vector< vertex_index > rank( num_vertices( filtered_mesh ));
    std::vector< vertex_descriptor > parent( num_vertices( filtered_mesh ), vertex_descriptor{0} );
    BOOST_FOREACH( vertex_descriptor vd, vertices( filtered_mesh ))
    {
      if( Triangle_mesh::null_vertex() != vd )
        parent[vd] = vd;

    }
    boost::disjoint_sets ds( &rank[0], &parent[0] );
    boost::incremental_components( filtered_mesh, ds );
    Components comp( parent.begin(), parent.end());
    vertex_index numCC = 0;
    BOOST_FOREACH( vertex_index cd, comp )
    {
      numCC = std::max( numCC, cd );
    }
    LOG_S( INFO ) << " number of CC : " << numCC;
    BOOST_FOREACH( vertex_index cd, comp )
    {
      CGAL::Face_around_target_iterator< Triangle_mesh > fbegin, fend;
      BOOST_FOREACH( vertex_index cc, comp[cd] )
      {
        auto vd = vertex_descriptor{cc};
        if( filter( vd ))
        {
          for( boost::tie( fbegin, fend ) = faces_around_target(
                 tmesh_.halfedge( vd ), tmesh_ );
               fbegin != fend;
               ++fbegin )
          {
            vvec[cd].insert( *fbegin );

          }
        }
      }
    }

    LOG_S( INFO ) << "End inserting CC";


    for( int i = 1; i <= numCC; ++i )
    {

      uintmax_t nv = 0, ne = 0;
      uintmax_t nf = vvec[i].size();


      if( nf > 35 )
      {
        std::tie( nv, ne ) = filtered_by_component( vvec[i] );
      }
      intmax_t euler = nv - ne + nf;
      eulers[i] = euler;          //tell you is ribbon is cut or perforated

    }

    LOG_S( INFO ) << "End Euler counting";

    return std::make_pair< std::unordered_map< int, int >,
                           std::unordered_map< int, std::set< face_descriptor > > >( std::move( eulers ),
                                                                                     std::move( vvec ));
  }

  std::pair< int, int > filtered_by_component( const std::set< face_descriptor > & face_set )
  {

    std::set< edge_descriptor > edges;
    std::set< vertex_descriptor > verts;
    for( const auto & fd: face_set )
    {
      //edges around faces
      for( auto [eb, ee] = edges_around_face( tmesh_.halfedge( fd ), tmesh_ ); eb != ee; ++eb )
      {
        edges.insert( *eb );
        verts.insert( boost::source( *eb, tmesh_ ));
        verts.insert( boost::target( *eb, tmesh_ ));
      }


    }

    return std::make_pair( verts.size(), edges.size());


  }

  void build_triple_map()
  {
    bool created;
    boost::tie( triple_map_, created ) = tmesh_.add_property_map< vertex_descriptor, double >( "v:triple", 0 );

    BOOST_FOREACH( vertex_descriptor vd, tmesh_.vertices())
    {
      get( triple_map_, vd ) = tags_projection( get( CGAL::vertex_point, tmesh_, vd ));
    }
  }

  //filter struct
  struct Is_triple
  {

    Is_triple() {};

    Is_triple( const Triangle_mesh::Property_map< vertex_descriptor, double > & mapin, const Triangle_mesh & tm,
               double tv )
      : m_tag_values( mapin ), m_tm( tm ), m_tv( tv ) {};
    Triangle_mesh::Property_map< vertex_descriptor, double > m_tag_values;
    Triangle_mesh m_tm;
    double m_tv;

    bool operator()( const vertex_descriptor & vd ) const
    {
      bool is_triple = get( m_tag_values, vd ) >= m_tv;
      return is_triple;
    }

    bool operator()( const edge_descriptor & ed ) const
    {

      auto vds = source( ed, m_tm ), vdt = target( ed, m_tm );
      return this->operator()( vds ) && this->operator()( vdt );
    }

  };

  double get_totalIntegral( const CurvatureMaps_type & gcurv )
  {
    double integrate;
    BOOST_FOREACH( face_descriptor fd, tmesh_.faces())
    {
      auto area = PMP::face_area( fd, tmesh_ );
      auto KG = 0.;
      auto count = 0;
      CGAL::Vertex_around_face_iterator< Triangle_mesh > vbegin, vend;
      for( boost::tie( vbegin, vend ) = vertices_around_face( tmesh_.halfedge( fd ), tmesh_ );
           vbegin != vend;
           ++vbegin, ++count )
      {
        KG += get( gcurv, *vbegin );
      }

      KG /= count;
      integrate += std::fabs( area ) * KG;
    }

    return integrate;
  }

  double get_angle( const std::set< face_descriptor > & cci, const CurvatureMaps_type & gcurv,
                    const std::function< double(double) > & fn, std::ofstream & csvfile )
  {
    double integrate = 0., total_area = 0.;

    //integrate
    for( const auto & fd: cci )
    {
      auto area = PMP::face_area( fd, tmesh_ );
      //value of gaussian -- interpolate properly or not
      auto KG = 0.;
      auto count = 0;
      CGAL::Vertex_around_face_iterator< Triangle_mesh > vbegin, vend;
      for( boost::tie( vbegin, vend ) = vertices_around_face( tmesh_.halfedge( fd ), tmesh_ );
           vbegin != vend;
           ++vbegin, ++count )
      {
        KG += get( gcurv, *vbegin );
      }


      KG /= count;
      integrate += std::fabs( area ) * KG;
      total_area += std::fabs( area );

    }

    LOG_S( INFO ) << "[" << integrate << "," << total_area << "]";
    for( const auto & fd: cci )
    {

      CGAL::Vertex_around_face_iterator< Triangle_mesh > vbegin, vend;
      for( boost::tie( vbegin, vend ) = vertices_around_face( tmesh_.halfedge( fd ), tmesh_ );
           vbegin != vend;
           ++vbegin )
      {
        const auto pt = get( CGAL::vertex_point, tmesh_, *vbegin );
        csvfile << pt[0] << ", " << pt[1] << "," << pt[2] << "," << fn( integrate ) << std::endl;
      }
    }


    return fn( integrate );
  }

  void build_tree( const std::vector< Point > & pts )
  {

    LOG_S( INFO ) << "Building neighboring tree";
    Image_point_property_map ppmap( pts );
    ptree_ = std::make_unique< Tree >( boost::counting_iterator< std::size_t >( 0 ),
                                       boost::counting_iterator< std::size_t >( pts.size()),
                                       Splitter(),
                                       Traits( ppmap ));
    p_tr_dist_ = std::make_unique< Distance >( ppmap );
  }

  double tags_projection( const Point & pt )
  {

    //for n-values
    K_neighbor_search search( *ptree_, pt /*query*/, 10 /*nb nearest neigh*/, 0, true, *p_tr_dist_,
                              true /*ordered*/ );
    auto const & pts = p_tr_dist_->point_property_map().points;
    //get max over the 3-closest points
    double avg = 0.;
    for( K_neighbor_search::iterator it = search.begin(); it != search.end(); it++ )
    {
      if( m_F.at( pts[it->first] ) == 2 )        //hit solid
      {
        avg = 2;
        break;
      }
      else
        avg = std::max( avg, m_F.at( pts[it->first] ));
    }
    return avg;
  }

  std::pair< CurvatureMaps_type, CurvatureMaps_type > process_curvature()
  {
    LOG_S( INFO ) << " Processing curvature";
    CurvatureMaps_type mcurv, gcurv;
    bool created;

    boost::tie( mcurv, created ) = tmesh_.add_property_map< vertex_descriptor, double >( "v:mcm", 0 );
    assert( created );
    boost::tie( gcurv, created ) = tmesh_.add_property_map< vertex_descriptor, double >( "v:gcm", 0 );
    assert( created );

    PMP::interpolated_corrected_curvatures( tmesh_,
                                            CGAL::parameters::vertex_mean_curvature_map( mcurv )
                                              .vertex_Gaussian_curvature_map( gcurv )
                                              .ball_radius( 0.5 ));

    return std::make_pair( mcurv, gcurv );
  }

  Curvature_PM_type process_refined( bool do_actually_refine )
  {

    //CUP
    Curvature_PM_type cpm;
    bool created;
    std::vector< double > delta;
    boost::tie( cpm, created ) = tmesh_.add_property_map< vertex_descriptor, std::vector< double > >( "v:cpm",
                                                                                                      {0, 0} );

    return cpm;

  }
};

}//end namespace

#endif /* ITF_CGAL_H_ */
