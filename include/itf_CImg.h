/*
 * itf_CImg.h
 *
 *  Created on: 20 déc. 2018
 *      Author: jfranc
 */

#ifndef ITF_CIMG_H_
#define ITF_CIMG_H_

/* CImg header */
#define cimg_use_tiff
#define cimg_use_openmp

#include <CImg.h>
#include <set>
#include <cassert>
#include <fstream>
#include <iostream>
#include <numeric>
#include <execution>

#include "../tpl/loguru/loguru.hpp"
/***
 * Interface template class which
 *  - input an image (cimg format)
 *  -  output a posttreated (cimg format) image via CIMG
 *
 *  TIPS: to work with Grayscale image mesher from CGAL
 *        output INR (*.inr)
 */
namespace ccpm
{

template< class T, class V >
class itf_to_CImg
{

public:
  itf_to_CImg(): processed_( false )
  {
    fname_ = "";
    LOG_S( INFO ) << "using itf to CImg template default ctor instantiation \n";
  }

  void set_input( const char *fname )
  {
    fname_ = fname;
    input_ = cimg_library::CImg< T >( fname_ );
  }

  int get_mapping( const std::string & prefix )
  {
    std::ofstream csv( prefix + ("mapping.csv"));
    csv << "x,y,z,n\n";
    auto grad = input_.get_gradient();
    LOG_S( INFO ) << "\n Mapping  info : " << grad.size() << " " << grad.max() << " " << grad.min() << std::endl;

    cimg_forXYZC( input_, x, y, z, c )
    {
      const uint64_t pos = input_.offset( x, y, z, c );
      std::set< T > neigh;
      int s = 0;

      if( grad[0][pos] + grad[1][pos] + grad[2][pos] != 0.0 )
      {
        csv << x << "," << y << "," << z;

        if( input_.offset( x - 1, y, z, c ) < input_.size())
          neigh.insert( input_[input_.offset( x - 1, y, z, c )] );
        if( input_.offset( x + 1, y, z, c ) < input_.size())
          neigh.insert( input_[input_.offset( x + 1, y, z, c )] );
        if( input_.offset( x, y - 1, z, c ) < input_.size())
          neigh.insert( input_[input_.offset( x, y - 1, z, c )] );
        if( input_.offset( x, y + 1, z, c ) < input_.size())
          neigh.insert( input_[input_.offset( x, y + 1, z, c )] );
        if( input_.offset( x, y, z - 1, c ) < input_.size())
          neigh.insert( input_[input_.offset( x, y, z - 1, c )] );
        if( input_.offset( x, y, z + 1, c ) < input_.size())
          neigh.insert( input_[input_.offset( x, y, z + 1, c )] );

        for( auto n: neigh )
          s += (int) n;
        csv << "," << (int) s;
        csv << "\n";

      }
    }

    return 0;

  }

  int to_isoValue( const std::string & prefix, const std::vector< int > & isoval, int N /*kept component*/ )
  {
    for( auto iso: isoval )
    {
      int delta = 1;
      cimg_library::CImg< V > img = (+input_);
      //discard non isoval
#pragma omp parallel for
      for( V *pv = img.data(); pv != img.data() + img.size(); ++pv )
      {
        (*pv > iso - delta && *pv < iso + delta) ? *pv = 1 : *pv = 0;
      }

      //discard some components
      img.erode( 5 );
      img.dilate( 5 );
      img.label();
      std::vector< std::size_t > size_list;

      //for each ccs count
#pragma omp parallel for
      for( int v1 = 0; v1 < img.max(); ++v1 )
      {
        int c = 0;
        std::for_each( img.data(), img.data() + img.size(), [v1, &c]( const V & vi ) { (vi == v1) ? ++c : 0; } );
        size_list.push_back( c );
#pragma omp critical
        LOG_S( INFO ) << "For iso " << iso << " push  " << c << " for  " << v1;
      }
      std::vector< V > index_list( size_list.size());
      std::iota( index_list.begin(), index_list.end(), 0 );
      //make the largest cc first
      std::sort( index_list.begin(), index_list.end(), [&size_list]( const V & v1, const V & v2 ) {
        return size_list[v1] < size_list[v2];
      } );

      //discard smaller cc s
#pragma omp parallel for
      for( auto *ptr = img.data(); ptr != (img.data() + img.size()); ++ptr )
      {
        *(ptr) = (std::find( std::next( index_list.begin(), N + 1 ), index_list.end(), *(ptr)) !=
                  index_list.end()) ? 0 : *(ptr);
      }
      img.save_tiff((prefix + std::to_string( iso ) + ("_cc.tiff")).c_str());
    }

    return 0;
  }

  int to_cc_images( const std::string & prefix, const std::vector< int > & isoval )
  {
    for( auto iso: isoval )
    {
      cimg_library::CImg< V > img;
      img.load_tiff((prefix + std::to_string( iso ) + ("_cc.tiff")).c_str());
      int count = 0;

#pragma omp parallel for
      for( int i = 1; i <= img.max(); ++i )
      {
        auto copy = (+img);
        int c = 0;
        std::for_each( copy.data(), copy.data() + copy.size(), [i, &c]( V & v ) { (v == i) ? ++c : v = 0; } );
#pragma omp critical
        if( c > 125 )
        {
          copy.save_tiff(
            (prefix + std::to_string( iso ) + ("_cc_") + std::to_string( i ) + (".tiff")).c_str());
          count++;
        }
      }

      LOG_S( INFO ) << "Splitting " << img.max() << " components ["<< iso <<"] into " << count << " images of more than 125 px";
    }


    return 0;

  }

  cimg_library::CImg< V >
  get_superposed( const std::string & output_dir, int n, const std::vector< double > & al,
                  const std::vector< double > & rot,
                  cimg_library::CImg< V >(&func)( int, double ))
  {
    assert( rot.size() == al.size() - 1 );
    std::string identity;
    cimg_library::CImg< V > img( n, n, n, 1, 0.0 );
    img += func( n, al[0] );
    for( int i = 0; i < rot.size(); ++i )
    {

      max_add( img, func( n, al[i + 1] ).get_rotate( rot[i], 1 /*nearest*/, 0 /*dirichlet*/ ), 3 );
      identity += std::to_string( rot[i] ) + "_" + std::to_string( al[i + 1] ) + "_";
    }

    //TODO change
    std::string cname = output_dir + "/contactSphere" + std::to_string( n ) + "_" + std::to_string( al[0] ) + "_";
    cname += identity;
    cname += ".tiff";
    img.save_tiff( cname.c_str());

    return img;

  }

  static cimg_library::CImg< V >
  get_cylinder( int n, double al )
  {
    cimg_library::CImg< V > img( n, n, n, 1, 0.0 );
    int r = static_cast< int >(double(30.0 / 128.0) * n), xc = img.width() / 2, yc = img.height() / 2, zc =
      img.depth() / 2;

    cimg_forXYZ( img, x, y, z )
    {
      if( y >= yc + al/2 *img.height() )                    //|| (y <= yc - al/2*img.height()) )
        img( x, y, z ) = 2.0;
      else if(((x - xc) * (x - xc) + (z - zc) * (z - zc) < r * r)
              && ( y < yc + al/2 * img.height()) )
//                        && ( y < yc + al/2 * img.height()|| y > yc - al/2 * img.height() ))
        img( x, y, z ) = 3.0;
      else
        img( x, y, z ) = 1.0;
    }


    return img;
  }

  static cimg_library::CImg< V >
  get_torus( int n, double al )
  {
    auto imgs = get_contactSphere( n, al ), imgc = get_cylinder( n, al );
    cimg_forXYZ( imgs, x, y, z )
    {

      if( imgc( x, y, z ) == 3.0 /*gas*/ && imgs( x, y, z ) == 3.0 )
      {
        imgs( x, y, z ) = 2.0;        /*solid*/
      }

    }

    return imgs;

  }


  void max_add( cimg_library::CImg< V > & img, const cimg_library::CImg< V > & other, int mask_value )
  {
    assert( img.size() == other.size());
    cimg_forXYZ( img, x, y, z )
    {
      if( img( x, y, z ) == mask_value && img( x, y, z ) != other( x, y, z ))
        img( x, y, z ) = other( x, y, z );
      else if( other( x, y, z ) == mask_value && img( x, y, z ) != other( x, y, z ))
        img( x, y, z ) = img( x, y, z );
      else
        img( x, y, z ) = std::max( img( x, y, z ), other( x, y, z ));
    }

  }

  static cimg_library::CImg< V > get_contactSphere( int n, double al )
  {
    cimg_library::CImg< V > img( n, n, n, 1, 0.0 );
    int r = static_cast< int >(double(60.0 / 128.0) * n), xc = img.width() / 2, yc = img.height() / 2, zc =
      img.depth() / 2;
    cimg_forXYZ( img, x, y, z )
    {
      if(((x - xc) * (x - xc) + (y - yc) * (y - yc) + (z - zc) * (z - zc) < r * r)
         && y < yc + al * r )
        img( x, y, z ) = 3.0;
      else if( y >= yc + al * r )
        img( x, y, z ) = 2.0;
      else
        img( x, y, z ) = 1.0;
    }


    return img;

  }


  const cimg_library::CImg< V > & get_output()
  {
    //do some processing
    process();
    return output_;
  }

  void release()
  {
    input_.clear();
  }

protected:
  ~itf_to_CImg() {};

private:
  cimg_library::CImg< T > input_;
  cimg_library::CImg< V > output_;

  const char *fname_;
  bool processed_;

  const char *rename()
  {
    std::string name = fname_;
    std::string v = "tiff";
    std::string newName;
    if( name.find( v ) != std::string::npos )
    {
      size_t lastindex = name.find_last_of( "." );
      newName = name.substr( 0, lastindex ) + (".inr");
    }
    return newName.c_str();
  }

  //black box function
  void process()
  {
    if( !processed_ )
    {
      LOG_S( INFO ) << " input [stats] (white/size)" << input_.get_normalize( 0, 1 ).sum() << " / " << input_.size()
                    << std::endl;
      output_ = input_.get_normalize( 0, 255 );
      LOG_S( INFO ) << " output [stats] (white/size)" << output_.sum() << " / " << output_.size() << std::endl;
      processed_ = true;
    }
  }

};

}


#endif /* ITF_CIMG_H_ */
