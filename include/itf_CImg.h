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

#include<CImg.h>
#include <set>
#include <fstream>
#include <iostream>
#include <numeric>
#include <execution>
/***
 * Interface template class which
 *  - input an image (cimg format)
 * 	-  output a posttreated (cimg format) image via CIMG
 *
 * 	TIPS: to work with Grayscale image mesher from CGAL
 * 	      output INR (*.inr)
 */
namespace ccpm {

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

        void get_mapping(const std::string &prefix) {
//             cimg_library::CImg<V> img = (+input_);
            std::ofstream csv(prefix + ("mapping.csv"));
            csv << "x,y,z,n\n";
            auto grad = input_.get_gradient();
            std::cout << "\n Mapping  info : " << grad.size() << " " << grad.max() << " " << grad.min() << std::endl;

            cimg_forXYZC(input_, x, y, z, c) {
                            const uint64_t pos = input_.offset(x, y, z, c);
                            std::set<T> neigh;
                            int s = 0;

                            if (grad[0][pos] + grad[1][pos] + grad[2][pos] != 0.0) {
                                csv << x << "," << y << "," << z;

                                if (input_.offset(x - 1, y, z, c) < input_.size())
                                    neigh.insert(input_[input_.offset(x - 1, y, z, c)]);
                                if (input_.offset(x + 1, y, z, c) < input_.size())
                                    neigh.insert(input_[input_.offset(x + 1, y, z, c)]);
                                if (input_.offset(x, y - 1, z, c) < input_.size())
                                    neigh.insert(input_[input_.offset(x, y - 1, z, c)]);
                                if (input_.offset(x, y + 1, z, c) < input_.size())
                                    neigh.insert(input_[input_.offset(x, y + 1, z, c)]);
                                if (input_.offset(x, y, z - 1, c) < input_.size())
                                    neigh.insert(input_[input_.offset(x, y, z - 1, c)]);
                                if (input_.offset(x, y, z + 1, c) < input_.size())
                                    neigh.insert(input_[input_.offset(x, y, z + 1, c)]);

                                for (auto n: neigh)
                                    s += (int) n;
                                csv << "," << (int) s;
                                csv << "\n";

                            }
                        }

        }

        void to_isoValue(const std::string &prefix, const std::vector<int> &isoval, int N /*kept component*/) {
            for (auto iso: isoval) {
                int delta = 1;
                cimg_library::CImg<V> img = (+input_);
#pragma omp parallel for
                for (V *pv = img.data(); pv != img.data() + img.size(); ++pv) {
                    (*pv > iso - delta && *pv < iso + delta) ? *pv = 1 : *pv = 0;
                }
//                std::for_each(img.data(), img.data() + img.size(),[iso,delta](T&v){ (v>iso-delta && v<iso+delta) ? v=1 : v=0;});
//                img.erode(3).dilate(3).label(true);

                std::vector<std::size_t> size_list;

                for (int v1 = 0; v1 < img.max(); ++v1) {
                    int c = 0;
                    std::for_each(img.data(), img.data() + img.size(), [v1, &c](const V &vi) { (vi == v1) ? ++c : 0; });
                    size_list.push_back(c);
                }
                std::vector<V> index_list(size_list.size());
                std::iota(index_list.begin(), index_list.end(), 0);
                std::sort(index_list.begin(), index_list.end(), [&size_list](const V &v1, const V &v2) {
                    return size_list[v1] > size_list[v2];
                });


////               #pragma omp parallel for
//                std::for_each( std::execution::par_unseq, img.data(), img.data() + img.size(), [&index_list,N](V& v)
//                    { v = (std::find(std::next(index_list.begin(),N+1),index_list.end(),v)!=index_list.end()) ? 0 : v;
//                        if(std::find(std::next(index_list.begin(),N+1),index_list.end(),v)!=index_list.end())
//                            std::cerr << "Discard " << v << "\n";
//                    });

#pragma omp parallel for
                for (auto *ptr = img.data(); ptr != (img.data() + img.size()); ++ptr) {
//                    if(std::find(std::next(index_list.begin(),N+1),index_list.end(),*(ptr))!=index_list.end())
//                        std::cerr << "Discard " << *ptr << "\n";
                    *(ptr) = (std::find(std::next(index_list.begin(), N + 1), index_list.end(), *(ptr)) !=
                              index_list.end()) ? 0 : *(ptr);
                }
                img.save_tiff((prefix + std::to_string(iso) + ("_cc.tiff")).c_str());
            }
        }

        void to_cc_images(const std::string &prefix, const std::vector<int> &isoval) {
            for (auto iso: isoval) {
                cimg_library::CImg<V> img;
                img.load_tiff((prefix + std::to_string(iso) + ("_cc.tiff")).c_str());

#pragma omp parallel for
                for (int i = 1; i <= img.max(); ++i) {
                    auto copy = (+img);
                    int c = 0;
                    std::for_each(copy.data(), copy.data() + copy.size(), [i, &c](V &v) { (v == i) ? ++c : v = 0; });
#pragma omp critical
                    if (c > 125)
                        copy.save_tiff(
                                (prefix + std::to_string(iso) + ("_cc_") + std::to_string(i) + (".tiff")).c_str());
                }
            }
        }

        void to_mlOtsu(int nclasses, const std::string &prefix) {
//           input_.blur_median(2);
            float m, M = input_.max_min(m);
            int nbins = 100;
            cimg_library::CImg<uint32_t> hist = input_.get_histogram(nbins, m, M), curr_ids(nclasses - 1), thres_ids(
                    nclasses - 1);
            cimg_library::CImg<float> zeroth_moment(hist.size()), first_moment(hist.size());

            zeroth_moment(0) = first_moment(0) = hist(0);
            for (int i = 1; i < hist.size(); ++i) {
                zeroth_moment(i) = zeroth_moment(i - 1) + hist(i);
                first_moment(i) = first_moment(i - 1) + i * hist(i);
            }

            std::function<float(uint32_t i, uint32_t j)> _get_var_btw_classes =
                    [&zeroth_moment, &first_moment](uint32_t i, uint32_t j) {
                        if (i == 0) {
                            if (zeroth_moment(i) > 0)
                                return std::pow(first_moment(j), 2) / zeroth_moment(j);
                        } else {

                            float zeroth_moment_ij = zeroth_moment(j) - zeroth_moment(i - 1);
                            if (zeroth_moment_ij > 0) {
                                float first_moment_ij = first_moment(j) - first_moment(i - 1);
                                return std::pow(first_moment_ij, 2) / zeroth_moment_ij;
                            }

                        }

                        return 0.;
                    };// end helper _get_var_btw_classes


            std::function<float(int, int, int, int, float &)> _get_thres_idx;
            _get_thres_idx = [&_get_var_btw_classes, &_get_thres_idx, &zeroth_moment, &first_moment, &curr_ids, &thres_ids]
                    (int hist_idx, int thres_idx, int nbins, int thres_count, float &sigma_max) {
                if (thres_idx < thres_count) {
                    for (int idx = hist_idx; idx < nbins - thres_count + thres_idx; ++idx) {
                        curr_ids(thres_idx) = idx;
                        sigma_max = _get_thres_idx(idx + 1, thres_idx + 1,
                                                   nbins, thres_count,
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
            _get_thres_idx(0, 0, nbins, nclasses - 1, sigma_max);

            // look into thres indices
            for (int i = 0; i < nclasses; ++i) {
                float min_as_val = (i == 0) ? 0 : (M - m) * (2 * thres_ids(i - 1) + 1) / (float) (2 * nbins),
                        max_as_val = (i == nclasses - 1) ? M : (M - m) * (2 * thres_ids(i) + 1) / (float) (2 * nbins);
                std::cout << "[" << min_as_val << "-" << max_as_val << "] writing "
                          << (prefix + std::to_string(i) + (".tiff")) << std::endl;
                cimg_library::CImg<T> img = (input_.get_threshold(min_as_val) &
                                             (cimg_library::CImg<T>(input_, false).fill(M) - input_).get_threshold(
                                                     M - max_as_val))
                        .get_erode(2).get_dilate(2).get_normalize(0, 255);

                img.save_tiff((prefix + std::to_string(i) + (".tiff")).c_str());
                if (i > 0) { // no need to label solid
                    img.label();
                    img.save_tiff((prefix + std::to_string(i) + ("_cc.tiff")).c_str());
                    for (int j = 0; j < img.max(); ++j) {

                        auto copy = (+img);
                        int c = 0;
                        std::for_each(copy.data(), copy.data() + copy.size(),
                                      [&c, j](float &v) { if (v == j) { ++c, v = 1; } else { v = 0; }});
                        std::cout << " cc : " << j << " count " << c << std::endl;
                        if (c > 15) {
                            copy.normalize(0, 255).save_tiff(
                                    (prefix + std::to_string(i) + ("_") + std::to_string(j) + ("_cc.tiff")).c_str());
                        }

                    }
                }

            }
        };

        const cimg_library::CImg <V> &
        get_contactSphere(int n, const std::vector<double> &al, const std::vector<double> &rot) {
            assert(rot.size() == al.size() - 1);
            std::string identity;
            cimg_library::CImg<V> img(n, n, n, 1, 0.0);
            img += get_contactSphere(n, al[0]);
            for (int i = 0; i < rot.size(); ++i) {

                max_add(img, get_contactSphere(n, al[i + 1]).get_rotate(rot[i], 1/*nearest*/, 0/*dirichlet*/), 3);
                identity += std::to_string(rot[i]) + "_" + std::to_string(al[i + 1]) + "_";
            }

            std::string cname = "/tmp/contactSphere" + std::to_string(n) + "_" + std::to_string(al[0]) + "_";
            cname += identity;
            cname += ".tiff";
            img.save_tiff(cname.c_str());

        }

        void max_add(cimg_library::CImg <V> &img, const cimg_library::CImg <V> &other, int mask_value) {
            assert(img.size() == other.size());
            cimg_forXYZ(img, x, y, z) {
                        if (img(x, y, z) == mask_value && img(x, y, z) != other(x, y, z))
                            img(x, y, z) = other(x, y, z);
                        else if (other(x, y, z) == mask_value && img(x, y, z) != other(x, y, z))
                            img(x, y, z) = img(x, y, z);
                        else
                            img(x, y, z) = std::max(img(x, y, z), other(x, y, z));
                    }

        }

        cimg_library::CImg <V> get_contactSphere(int n, double al) {
            cimg_library::CImg<V> img(n, n, n, 1, 0.0);
            int r = static_cast<int>(double(60.0 / 128.0) * n), xc = img.width() / 2, yc = img.height() / 2, zc =
                    img.depth() / 2;
            cimg_forXYZ(img, x, y, z) {
                        if (((x - xc) * (x - xc) + (y - yc) * (y - yc) + (z - zc) * (z - zc) < r * r)
                            && y < yc + al * r)
                            img(x, y, z) = 3.0;
                        else if (y >= yc + al * r)
                            img(x, y, z) = 2.0;
                        else
                            img(x, y, z) = 1.0;
                    }


            return img;

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
                std::cerr << " input [stats] (white/size)" << input_.get_normalize(0, 1).sum() << " / " << input_.size()
                          << std::endl;
                output_ = input_.get_normalize(0, 255);
//                output_.blur(1);
                std::cerr << " output [stats] (white/size)" << output_.sum() << " / " << output_.size() << std::endl;
                processed_ = true;
            }
        }

    };

}


#endif /* ITF_CIMG_H_ */
