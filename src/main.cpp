//////////////////////////////////////////////////////////////////////////////////
//
// Created by jfranc on Mon 06 2024.
//////////////////////////////////////////////////////////////////////////////////

#include <iostream>
//argparser
#include <cxxopts.hpp>
#include <boost/filesystem.hpp>

#include "itf_policies.h"
#include "itf_CImg.h"
#include "itf_CGAL.h"

std::tuple<string,string,string> split_name(const string& fname)
{
         size_t premindex = fname.find_last_of("/");
         size_t lastindex = fname.find_last_of(".");
         return std::make_tuple(fname.substr(0,premindex), fname.substr(premindex, lastindex-premindex), fname.substr(lastindex,fname.size()));
}


int main(int argc, const char **argv) {
    //block options
    cxxopts::Options options("ccpm", "Curvature capture in porous media");
    options.set_width(70)
            .set_tab_expansion()
            .allow_unrecognised_options()
            .add_options()
                    ("d,debug", "Enable various level of verbosity and debug",
                     cxxopts::value<int>()->default_value("0"))
                    ("i,image", "image file name", cxxopts::value<std::string>())
                    ("t,test", "testing", cxxopts::value<std::string>())
                    ("ard", "surface mesh definition",
                     cxxopts::value<std::vector<double>>()->default_value("30.,2.,1."))
                    ("o,output", "output repository", cxxopts::value<std::string>())
                    ("h,help", "Print usage");


    auto options_parsed = options.parse(argc, argv);
    //create the repo
    //help print
    if (options_parsed["help"].count()) {
        std::cout << options.help() << std::endl;
        exit(0);
    }
    else if(options_parsed.count("output")>0) {
        boost::filesystem::path output_dir = boost::filesystem::path(options_parsed["output"].as<std::string>());
        boost::filesystem::create_directory(output_dir);

        //req. interfaces
        ccpm::interface<ccpm::itf_to_CImg<float, float>> Itf;
        string prefix, base,  ext;
        std::cout << "image file " << options_parsed["image"].as<std::string>() << std::endl;
        Itf.set_input(options_parsed["image"].as<std::string>().c_str());
        auto names = split_name(options_parsed["image"].as<std::string>().c_str());

        if (options_parsed.count("test")>0){

            Itf.to_mlOtsu(3, output_dir.string() + ("/otsu_"));
            return 0;
        }


        ccpm::interface<ccpm::itf_to_CGAL> Itf2;


        Itf.get_output().save_inr("from-image.inr");
        Itf2.set_input("from-image.inr");
        Itf2.set_refined();
        Itf2.set_ARD(
                options_parsed["ard"].as<std::vector<double>>()[0],
                options_parsed["ard"].as<std::vector<double>>()[1],
                options_parsed["ard"].as<std::vector<double>>()[2]);

        Itf2.save_surf_stl(boost::filesystem::path(output_dir / ( std::get<1>(names) + (".stl"))).c_str());
    }
    else
        throw std::invalid_argument("Output directory required");



    return 0;
}