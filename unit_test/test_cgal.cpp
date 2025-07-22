//////////////////////////////////////////////////////////////////////////////////
//
// Created by jfranc on Wed 05 2025.
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
#include <iostream>

#define LOGURU_WITH_STREAMS 1

#include "itf_policies.h"
#include "itf_CGAL.h"
#include <gtest/gtest.h>

#include "../tpl/loguru/loguru.cpp"

TEST( TestCGAL, test_all )
{
  ccpm::interface< ccpm::itf_to_CGAL > cgal_interface;
  int ret = 0;
  EXPECT_EQ( cgal_interface.set_input( "/opt/output/from-image.inr" ), 0 ) << " Error reading input ";
  ret = cgal_interface.input_neighFile( std::string( "/opt/output/isoVal_mapping.csv" ).c_str());
  EXPECT_EQ( ret, 0 ) << "Import of neighboring data failed" << std::endl;
  //
  cgal_interface.set_refined();
  cgal_interface.set_ARD( 30., 4., 1. );
  cgal_interface.set_cutoff( 5. );
  std::vector< double > angles, refval{26, 36};
  ret = cgal_interface.save_surf_stl( std::string( "/opt/output/isoVal_3_cc_1.stl" ).c_str(), angles );
  EXPECT_EQ( ret, 0 ) << "Curvature meas and integration failed " << std::endl;
  std::sort( std::begin( angles ), std::end( angles ));
  for( auto iangle = std::cbegin( angles ), iref = std::cbegin( refval );
       iangle != std::cend( angles ) && iref != std::cend( refval ); iangle++, iref++ )
    EXPECT_NEAR( *iangle, *iref, 10 ) << *iangle << " is not equal to refval " << *iref;
}

int main( int argc, char *argv[] )
{
  testing::InitGoogleTest( &argc, argv );
  return RUN_ALL_TESTS();
}
