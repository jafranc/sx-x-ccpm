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
#include<CImg.h>

/***
 * Interface template class which
 *  - input an image (cimg format)
 * 	-  output a posttreated (cimg format) image via CIMG
 *
 * 	TIPS: to work with Grayscale image mesher from CGAL
 * 	      output INR (*.inr)
 */
namespace ccpm{

    template<class T, class V>
    class itf_to_CImg {

    public:
        itf_to_CImg() : processed_(false) {
            fname_ = "";
            std::cout << "using itf to CImg template default ctor instantiation \n";
        }

        void set_input(const char *fname) {
            fname_ = fname;
            input_ = cimg_library::CImg<T>(fname_);
        }

        const cimg_library::CImg <T> &get_input() const {
            return input_;
        }

        const cimg_library::CImg <V> &get_output() {
            //do some processing
            process();
            return output_;
        }

    protected:
        ~itf_to_CImg() {};

    private:
        cimg_library::CImg <T> input_;
        cimg_library::CImg <V> output_;
        const char *fname_;
        bool processed_;

        const char *rename() {
            std::string name = fname_;
            std::string v = "tiff";
            std::string newName;
            if (name.find(v) != std::string::npos) {
                size_t lastindex = name.find_last_of(".");
                newName = name.substr(0, lastindex) + (".inr");
            }
            return newName.c_str();
        }

        //black box function
        void process() {
            if (!processed_) {
                output_ = input_.normalize(0, 255);
                output_.blur(1);
                processed_ = true;
            }
        }

    };

}


#endif /* ITF_CIMG_H_ */
