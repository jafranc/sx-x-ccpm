//////////////////////////////////////////////////////////////////////////////////
//
// Created by jfranc on Thu 03 2025.
//////////////////////////////////////////////////////////////////////////////////

#ifndef CCPM_UTILITIES_HPP
#define CCPM_UTILITIES_HPP


#include <vector>

namespace ccpm {
    namespace utility {


/*        struct UnionFind {
            // A vector of parents.
            std::vector<int> parent_;
            std::vector<uint32_t> rank_;

            // Constructs a union-find data structure with n elements.
            UnionFind(size_t n) {
                // Clear the parent vector.
                parent_.clear();
                rank_.clear();

                // Resize the parent vector to n elements.
                parent_.resize(n);
                rank_.resize(n);

                // Fill the parent vector with -1.
                std::fill(parent_.begin(), parent_.end(), -1);
                std::fill(rank_.begin(), rank_.end(), 0);
            }

            std::tuple<int,int> get_max_rank()
            {
                auto mit = std::max_element(rank_.begin(),rank_.end());
                return std::make_tuple(*mit,std::distance(rank_.begin(),mit));
            }

            std::vector< std::vector<int> > get_ordered_component(int n=-1)
            {          //order all
                std::vector<int> idx(rank_.size()); std::iota(idx.begin(),idx.end(),0);
                std::sort(idx.begin(),idx.end(),[this](int r1, int r2){ return (rank_[r2]<rank_[r1]); });

                std::vector< std::vector<int> > out((n>0) ? n : rank_.size());
                for (int i = 0; i < out.size(); ++i)
                    for(int ii = 0; ii< parent_.size(); ++ii)
                    {
                        if(set_find(parent_[ii])==idx[i])
                            out[i].push_back(ii);
                    }

                return out;

                }

            // Creates a new set with the given element.
            size_t make_set(size_t k) {
                // Set the parent of the given element to itself.
                parent_[k] = k;
                rank_[k] = 0;

                // Return the given element.
                return k;
            }

            // Unites the sets containing the given elements.
            void set_union(size_t i, size_t j) {
                // Set the parent of the element i to the element j.
                if (parent_[i] != parent_[j]) {
                    if (rank_[i] < rank_[j]) {
                        parent_[i] = j;
                        for(auto& p : parent_)
                            if (p == i) p = j;
                        rank_[j]+=rank_[i]+1;
                    } else {
                        parent_[j] = i;
                        for(auto& p : parent_)
                            if (p == i) p = j;
                        rank_[i]+=rank_[j]+1;
                    }
                }
            }

            // Finds the representative of the set containing the given element.
            size_t set_find(size_t k) {
                // If the parent of the given element is itself, then the given element is the representative.
                if (parent_[k] == k) {
                    return k;
                }

                // Otherwise, recursively find the representative of the set containing the parent of the given element.
                return set_find(parent_[k]);
            }
        };*/

//        class UnionFind {
//
//        public:
//            UnionFind() = delete;
//            UnionFind(int id):parent_(this),rank_(0),id_(id){};
//
//            UnionFind(const UnionFind &) = default;
//            ~UnionFind(){};
//
//            bool operator<(const UnionFind &other) const {
//                return rank_ < other.rank_;
//            }
//
//            bool operator>(const UnionFind &other) const {
//                return rank_ > other.rank_;
//            }
//
////            bool operator++() {
////                rank_++;
////            }
//
//            static std::shared_ptr<UnionFind> find(std::shared_ptr<UnionFind> u) {
//
//                if (u->parent_.get() != u.get())
//                    u->setParent(find(u->parent_));
//                return u->parent_;
//
//            }
//
//            inline int getRank() const { return rank_; }
//
//            static void makeUnion(std::shared_ptr<UnionFind> u0, std::shared_ptr<UnionFind> u1) {
//
//                u0 = find(u0);
//                u1 = find(u1);
//
//                if (u0 == u1) {
//                    return;
//                }
////                else if (u0 > u1) {
////                    u1->setParent(u0);
////                    return u0;
////                }
//                else if (u1 > u0) {
//                    u0->setParent(u1);
//                    u1->rank_+=1;
//                    return;
//                } else {
//                    u1->setParent(u0);
//                    u0->rank_+=1;
//                    return;
//                }
//
//            }
//
////            static std::shared_ptr<UnionFind>  makeUnion(std::vector<std::shared_ptr<UnionFind>> &sets) {
////                std::shared_ptr<UnionFind> n = nullptr;
////
////                if (sets.empty())
////                    return n;
////                else if (sets.size() == 1)
////                    return sets[0];
////                else {
////                    for (int i = 0; i < sets.size() - 1; ++i) {
////                        n = makeUnion(sets[i], sets[i + 1]);
////                    }
////
////                    return n;
////                }
////            }
//
//            void setParent(const std::shared_ptr<UnionFind> u) {
////                parent_.reset(u.get());
//                parent_ = u;
//            }
//
//
//        private:
//            std::shared_ptr<UnionFind> parent_;
//            int rank_;
//            const int id_;
//        };
//
    } // utility
} // ccpm

#endif //CCPM_UTILITIES_HPP
