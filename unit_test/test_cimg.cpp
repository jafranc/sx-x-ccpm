//////////////////////////////////////////////////////////////////////////////////
//
// Created by jfranc on Tue 05 2025.
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
#include <iostream>

#define LOGURU_WITH_STREAMS 1

#include "itf_policies.h"
#include "itf_CImg.h"

#include "../tpl/loguru/loguru.cpp"

#include <gtest/gtest.h>

TEST( TestCimg, test_generator )
{

  ccpm::interface< ccpm::itf_to_CImg< u_int8_t, u_int16_t > > cimg_interface;
  //
  double expected_vol = cimg_interface.get_superposed( std::string( "/opt/output" ),
                                                       128 /*size*/,
                                                       std::vector< double >{0.9, 0.8},
                                                       std::vector< double >{180},
                                                       decltype(cimg_interface)::get_contactSphere
                                                       ).sum();

  EXPECT_DOUBLE_EQ( expected_vol, 4241680 ) << " Expected Sphere Volume is not satisfied \n \t theoretical value: " << 4241680
                                            << "\n\t calculated value: " << expected_vol << std::endl;

  //
  expected_vol = cimg_interface.get_superposed( std::string( "/opt/output" ),
                                                128 /*size*/,
                                                std::vector< double >{0.9},
                                                std::vector< double >{},
                                                decltype(cimg_interface)::get_cylinder
                                                ).sum();

  EXPECT_DOUBLE_EQ( expected_vol, 2880852 ) << " Expected Cylinder Volume is not satisfied \n \t theoretical value: " << 2880852
                                            << "\n\t calculated value: " << expected_vol << std::endl;
  //
  expected_vol = cimg_interface.get_superposed( std::string( "/opt/output" ),
                                                64 /*size*/,
                                                std::vector< double >{0.9},
                                                std::vector< double >{},
                                                decltype(cimg_interface)::get_torus
                                                ).sum();

  EXPECT_DOUBLE_EQ( expected_vol, 468284 ) << " Expected Cylinder Volume is not satisfied \n \t theoretical value: " << 468284
                                           << "\n\t calculated value: " << expected_vol << std::endl;

}

TEST( TestCimg, test_decompose )
{
  ccpm::interface< ccpm::itf_to_CImg< u_int8_t, u_int16_t > > cimg_interface;
  cimg_interface.set_input(
    std::string( "/opt/output/contactSphere128_0.900000_180.000000_0.800000_.tiff" ).c_str());
  int ret = 0;
  //if isoval, we generate label lists and isolate each component and their CC
  {
    ret = cimg_interface.get_mapping( std::string( "/opt/output/isoVal_" ).c_str());
    EXPECT_EQ( ret, 0 ) << "Mapping Failed " << std::endl;
    ret = cimg_interface.to_isoValue( std::string( "/opt/output/isoVal_" ).c_str(),
                                      std::vector< int >{1, 2, 3}, 3 );
    EXPECT_EQ( ret, 0 ) << "Isovalue Failed " << std::endl;
    ret = cimg_interface.to_cc_images( std::string( "/opt/output/isoVal_" ).c_str(),
                                       std::vector< int >{1, 2, 3} );
    EXPECT_EQ( ret, 0 ) << "CC Failed " << std::endl;
  }


}

TEST( TestCimg, test_INRize )
{

  ccpm::interface< ccpm::itf_to_CImg< u_int8_t, u_int16_t > > cimg_interface;
  cimg_interface.set_input(
    std::string( "/opt/output/isoVal_3_cc_1.tiff" ).c_str());

  //part of CIMG testing guaranteed -- no need to test
  cimg_interface.get_output().save_inr( "/opt/output/from-image.inr" );
}

int main( int argc, char *argv[] )
{
  testing::InitGoogleTest( &argc, argv );
  return RUN_ALL_TESTS();
}
