/*
 * Poly_rings.h
 *
 *  Created on: 14 févr. 2019
 *      Author: jfranc
 */

#ifndef POLY_RINGS_H_
#define POLY_RINGS_H_

#include <cassert>

using namespace std;

template<class SurfaceMesh, class VertexPropertyMap>
class Poly_rings {
protected:
//  //Polyhedron
//  typedef typename ::Vertex Vertex;
//  typedef typename ::Halfedge Halfedge;
//  typedef typename ::Facet Facet;
//  typedef typename ::Halfedge_around_vertex_circulator  Halfedge_around_vertex_circulator;
//  typedef typename ::Vertex_iterator Vertex_iterator;

    typedef typename boost::graph_traits<SurfaceMesh>::vertex_descriptor vertex_descriptor;
    //vertex property map
//   typedef typename boost::property_traits<VertexPropertyMap>::value_type vpm_value_type;
//   typedef typename boost::property_traits<VertexPropertyMap>::key_type vpm_key_type;

public:
    //vertex indices are initialised to -1
    static void reset_ring_indices(std::vector<vertex_descriptor> &vces,
                                   VertexPropertyMap &vpm);

    //i >= 1; from a start vertex on the current i-1 ring, push non-visited neighbors
    //of start in the nextRing and set indices to i. Also add these vertices in all.
    static void push_neighbours_of(const SurfaceMesh &m,
                                   vertex_descriptor start,
                                   int ith,
                                   std::vector<vertex_descriptor> &nextRing,
                                   std::vector<vertex_descriptor> &all,
                                   VertexPropertyMap &vpm);

    //i >= 1, from a currentRing i-1, collect all neighbors, set indices
    //to i and store them in nextRing and all.
    static void collect_ith_ring(const SurfaceMesh &m, int ith,
                                 std::vector<vertex_descriptor> &currentRing,
                                 std::vector<vertex_descriptor> &nextRing,
                                 std::vector<vertex_descriptor> &all,
                                 VertexPropertyMap &vpm);

public:

    static void collect_i_rings(const SurfaceMesh &m, vertex_descriptor vd,
                                int ring_i,
                                std::vector<vertex_descriptor> &all,
                                VertexPropertyMap &vpm);

    static void collect_enough_rings(const SurfaceMesh &m,
                                     vertex_descriptor vd,
                                     unsigned int min_nb,
                                     std::vector<vertex_descriptor> &all,
                                     VertexPropertyMap &vpm);


};


template<class SurfaceMesh, class VertexPropertyMap>
void Poly_rings<SurfaceMesh, VertexPropertyMap>::
reset_ring_indices(std::vector<vertex_descriptor> &vces, VertexPropertyMap &vpm) {
    BOOST_FOREACH(vertex_descriptor vd, vces) vpm[vd] = -1;
}

template<class SurfaceMesh, class VertexPropertyMap>
void Poly_rings<SurfaceMesh, VertexPropertyMap>::
push_neighbours_of(const SurfaceMesh &m,
                   vertex_descriptor start,
                   int ith,
                   std::vector<vertex_descriptor> &nextRing,
                   std::vector<vertex_descriptor> &all,
                   VertexPropertyMap &vpm) {

    BOOST_FOREACH(vertex_descriptor vd, vertices_around_target(start, m)) if (vpm[vd] !=
                                                                              -1) { continue; }//if visited: next
                    else {
                        vpm[vd] = ith;
                        nextRing.push_back(vd);
                        all.push_back(vd);
                    }
}

template<class SurfaceMesh, class VertexPropertyMap>
void Poly_rings<SurfaceMesh, VertexPropertyMap>::
collect_ith_ring(const SurfaceMesh &m, int ith,
                 std::vector<vertex_descriptor> &currentRing,
                 std::vector<vertex_descriptor> &nextRing,
                 std::vector<vertex_descriptor> &all,
                 VertexPropertyMap &vpm) {

    BOOST_FOREACH(vertex_descriptor vd, currentRing) push_neighbours_of(m, vd, ith, nextRing, all, vpm);

}

template<class SurfaceMesh, class VertexPropertyMap>
void Poly_rings<SurfaceMesh, VertexPropertyMap>::
collect_i_rings(const SurfaceMesh &m, vertex_descriptor vd,
                int ring_i,
                std::vector<vertex_descriptor> &all,
                VertexPropertyMap &vpm) {
    std::vector<vertex_descriptor> current_ring, next_ring;

    assert(ring_i >= 1);
    //initialize
    vpm[vd] = 0;
    current_ring.push_back(vd);
    all.push_back(vd);

    for (int i = 1; i <= ring_i; i++) {
        collect_ith_ring(m, i, current_ring, next_ring, all, vpm);
        //next round must be launched from p_nextRing...
        current_ring.clear();
        std::swap(current_ring, next_ring);
    }
    //clean up
    reset_ring_indices(all, vpm);

}

template<class SurfaceMesh, class VertexPropertyMap>
void Poly_rings<SurfaceMesh, VertexPropertyMap>::
collect_enough_rings(const SurfaceMesh &m,
                     vertex_descriptor vd,
                     unsigned int min_nb,
                     std::vector<vertex_descriptor> &all,
                     VertexPropertyMap &vpm) {
    //declare
    std::vector<vertex_descriptor> current_ring, next_ring;
    //initialize
    vpm[vd] = 0;
    current_ring.push_back(vd);
    all.push_back(vd);

    int i = 1;
    while ((all.size() < min_nb) && (current_ring.size() != 0)) {
        collect_ith_ring(m, i, current_ring, next_ring, all, vpm);
        //next round must be launched from p_nextRing...
        current_ring.clear();
        std::swap(current_ring, next_ring);
        i++;
    }
    //clean up
    reset_ring_indices(all, vpm);
}

#endif /* POLY_RINGS_H_ */
