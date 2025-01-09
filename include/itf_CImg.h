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


        void to_mlOtsu(int nclasses, const std::string& prefix)
        {
           input_.blur_median(2);
           float m, M = input_.max_min(m);
           int nbins = 100;
           cimg_library::CImg<uint32_t> hist = input_.get_histogram(nbins, m, M),curr_ids(nclasses-1), thres_ids(nclasses-1);
           cimg_library::CImg<float> zeroth_moment(hist.size()), first_moment(hist.size());

           zeroth_moment(0) = first_moment(0) = hist(0);
           for(int i=1; i<hist.size(); ++i) {
               zeroth_moment(i) = zeroth_moment(i - 1) + hist(i);
               first_moment(i) = first_moment(i - 1) + i*hist(i);
           }

            std::function<float(uint32_t i, uint32_t j)> _get_var_btw_classes =
                    [&zeroth_moment, &first_moment](uint32_t i, uint32_t j)
            {
               if(i==0) {
                   if (zeroth_moment(i) > 0)
                       return std::pow(first_moment(j), 2) / zeroth_moment(j);
               }
               else{

                   float zeroth_moment_ij = zeroth_moment(j) - zeroth_moment(i-1);
                   if(zeroth_moment_ij>0) {
                       float first_moment_ij = first_moment(j) - first_moment(i - 1);
                       return std::pow(first_moment_ij, 2) / zeroth_moment_ij;
                   }

               }

                return 0.;
            };// end helper _get_var_btw_classes


            std::function<float(int,int,int,int,float&)> _get_thres_idx;
            _get_thres_idx = [&_get_var_btw_classes, &_get_thres_idx,&zeroth_moment,&first_moment,&curr_ids,&thres_ids]
                    (int hist_idx, int thres_idx, int nbins, int thres_count, float& sigma_max) {
                if (thres_idx < thres_count) {
                    for(int idx=hist_idx; idx < nbins - thres_count + thres_idx; ++idx) {
                        curr_ids(thres_idx) = idx;
                        sigma_max = _get_thres_idx(idx+1,thres_idx+1,
                                                   nbins,thres_count,
                                                   sigma_max);
                    }

                } else {
                    float sigma = (_get_var_btw_classes(0, curr_ids(0))) +
                                  (_get_var_btw_classes(curr_ids[thres_count - 1] + 1, nbins - 1));

                    for (int idx = 0; idx < thres_count - 1; ++idx) {
                        sigma += _get_var_btw_classes(curr_ids[idx] + 1, curr_ids[idx + 1]);
                    }

                    if (sigma > sigma_max) {
                        sigma_max = sigma;
                        thres_ids = curr_ids;
                    }
                }

                return sigma_max;
            };//end of _get_thres_idx helper

            float sigma_max = 0.;
            _get_thres_idx(0, 0, nbins, nclasses-1, sigma_max);

            // look into thres indices
            for (int i = 0; i < nclasses; ++i) {
                float min_as_val = (i==0) ? 0 : (M-m)*(2*thres_ids(i-1)+1)/(float)(2*nbins),
                        max_as_val = (i==nclasses-1) ? M : (M-m)*(2*thres_ids(i)+1)/(float)(2*nbins);
                std::cout <<"[" << min_as_val << "-" << max_as_val << "] writing " << (prefix + std::to_string(i) + (".png")) << std::endl;
                (input_.get_threshold(min_as_val) &
                (cimg_library::CImg<T>(input_,false).fill(M)-input_).get_threshold(M-max_as_val))
                    .get_erode(0).get_dilate(0).get_normalize(0,255).save_png( (prefix + std::to_string(i) + (".png")).c_str() );
            }


            };


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
                output_.blur(5);
                processed_ = true;
            }
        }

    };

}


#endif /* ITF_CIMG_H_ */
