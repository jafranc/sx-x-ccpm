//////////////////////////////////////////////////////////////////////////////////
//
// Created by jfranc on Mon 06 2024.
//////////////////////////////////////////////////////////////////////////////////

#include <iostream>
//argparser
#include <cxxopts.hpp>
#include <boost/filesystem.hpp>

#define LOGURU_WITH_STREAMS 1

#include "../tpl/loguru/loguru.cpp"

#include "itf_policies.h"
#include "itf_CImg.h"
#include "itf_CGAL.h"

std::tuple< std::string, std::string, std::string > split_name( const std::string & fname )
{
  size_t premindex = fname.find_last_of( "/" );
  size_t lastindex = fname.find_last_of( "." );
  return std::make_tuple( fname.substr( 0, premindex ), fname.substr( premindex, lastindex - premindex ),
                          fname.substr( lastindex + 1, fname.size()));
}


int main( int argc, char * *argv )
{
  //loguru
  loguru::init( argc, argv );
  loguru::add_file( "ccpm.log", loguru::Append, loguru::Verbosity_MAX );

  //block options
  cxxopts::Options options( "ccpm", "Curvature capture in porous media" );
  options.set_width( 70 )
    .set_tab_expansion()
    .allow_unrecognised_options()
    .add_options()
    ( "d,debug", "Enable various level of verbosity and debug",
    cxxopts::value< int >()->default_value( "0" ))
  //option for image pre-treatment generating images on hte values input
    ( "i,image", "image file name", cxxopts::value< std::string >())
    ( "u,isoval", "isovalues at which isolate phases", cxxopts::value< std::vector< int > >())
  // - mesh related
    ( "ard", "surface mesh definition",
    cxxopts::value< std::vector< double > >()->default_value( "30.,2.,1." ))
    ( "c,ncomponent", "number of components", cxxopts::value< int >()->default_value( "4" ))
  //options for curvature isolations
    ( "n,ncoords", "neighborhood coord file name", cxxopts::value< std::string >())
    ( "x,cutoff", "isovalues above which cutoff the boundary", cxxopts::value< double >())
  //droplets image generators options
    ( "t,gen-drop", "generate droplets with multiple angles", cxxopts::value< int >())
    ( "a,test-angle", "testing-ange", cxxopts::value< std::vector< double > >())           //left in non-production
    ( "r,test-rot", "testing-rot", cxxopts::value< std::vector< double > >())           //left in non-production
  //misc / usual
    ( "o,output", "output repository", cxxopts::value< std::string >())
    ( "h,help", "Print usage" );


  auto options_parsed = options.parse( argc, argv );
  //create the repo
  //help print
  if( options_parsed["help"].count())
  {
    std::cout << options.help() << std::endl;
    exit( 0 );
  }
  //workflow options
  else if( options_parsed.count( "output" ) > 0 )
  {

    boost::filesystem::path output_dir = boost::filesystem::path( options_parsed["output"].as< std::string >());
    boost::filesystem::create_directory( output_dir );
    //req. interfaces
    ccpm::interface< ccpm::itf_to_CImg< u_int8_t, u_int16_t > > cimg_interface;
    //droplets genertions
    if( options_parsed.count( "gen-drop" ) > 0 && options_parsed.count( "test-angle" ) > 0 &&
        options_parsed.count( "test-rot" ) > 0 )
    {
      cimg_interface.get_superposed( output_dir.string(), options_parsed["gen-drop"].as< int >(),
                                     options_parsed["test-angle"].as< std::vector< double > >(),
                                     options_parsed["test-rot"].as< std::vector< double > >(),
                                     ccpm::interface< ccpm::itf_to_CImg< u_int8_t, u_int16_t > >::get_contactSphere );

      return 0;
    }
    else
    {

      std::string prefix, base;
      auto names = split_name( options_parsed["image"].as< std::string >().c_str());
      size_t lastindex = options_parsed["image"].as< std::string >().find_last_of( "." );
      auto ext = std::get< 2 >( names );

      ccpm::interface< ccpm::itf_to_CGAL > cgal_interface;

      if( ext == "tiff" || ext == "tif" )
      {
        LOG_S( INFO ) << "image file " << options_parsed["image"].as< std::string >() << std::endl;
        cimg_interface.set_input( options_parsed["image"].as< std::string >().c_str());
        //if isoval, we generate label lists and isolate each component and their CC
        if( options_parsed.count( "isoval" ) > 0 )
        {
          cimg_interface.get_mapping( output_dir.string() + ("/isoVal_"));
          cimg_interface.to_isoValue( output_dir.string() + ("/isoVal_"),
                                      options_parsed["isoval"].as< std::vector< int > >(),
                                      options_parsed["ncomponent"].as< int >());
          cimg_interface.to_cc_images( output_dir.string() + ("/isoVal_"),
                                       options_parsed["isoval"].as< std::vector< int > >());
          return 0;
        }

        cimg_interface.get_output().save_inr( "from-image.inr" );
        cgal_interface.set_input( "from-image.inr" );       // inr format is the input of SurfaceMesh
      }
      else if( ext == "off" )
      {
        cgal_interface.set_input( options_parsed["image"].as< std::string >().c_str());
      }


      cimg_interface.release();
      cgal_interface.input_neighFile( options_parsed["ncoords"].as< std::string >());
      cgal_interface.set_refined();
      cgal_interface.set_ARD(
        options_parsed["ard"].as< std::vector< double > >()[0],
        options_parsed["ard"].as< std::vector< double > >()[1],
        options_parsed["ard"].as< std::vector< double > >()[2] );
      cgal_interface.set_cutoff( options_parsed["cutoff"].as< double >());
      std::vector< double > angles;
      cgal_interface.save_surf_stl( boost::filesystem::path( output_dir / (std::get< 1 >( names ) + (".stl"))).c_str(),
                                    angles );

    }
  }
  else
    throw std::invalid_argument( "Output directory required" );


  return 0;
}
