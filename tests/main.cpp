//////////////////////////////////////////////////////////////////////////////////
//
// Created by jfranc on Fri 04 2025.
//////////////////////////////////////////////////////////////////////////////////
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/foreach.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/incremental_components.hpp>
#include <boost/property_map/dynamic_property_map.hpp>
#include <boost/graph/subgraph.hpp>

#include <iostream>

struct node_t {
    typedef boost::vertex_property_tag kind;
};


template<int NP, int NC>
class Network {

public:

//typedef
//boost::subgraph<
//boost::adjacency_list<
//        boost::vecS            // edge list
//        , boost::vecS            // vertex list
//        , boost::undirectedS     // directedness
//        , boost::property<node_t, int> // property associated with vertices
//        , boost::property<boost::edge_index_t, int>
//        >>
//        Graph;
    using vecS = boost::vecS;
    typedef boost::subgraph<boost::adjacency_list<vecS, vecS, boost::undirectedS, boost::property<node_t, int>, boost::property<boost::edge_index_t, int>>> Graph;

    typedef boost::graph_traits<Graph>::edge_descriptor edge_descriptor;
    typedef boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
    typedef boost::graph_traits<Graph>::vertices_size_type VertexIndex;
//typedef boost::property_map<Graph, boost::vertex_index_t>::type PMap;
    typedef std::vector<int> PMap;
    typedef boost::component_index<VertexIndex> Components;

    struct Is_overflown {

        Is_overflown() {};

        Is_overflown(const PMap &map, const Graph &tm) : m_map(map), m_tm(tm) {};
        PMap m_map;
        Graph m_tm;

        // -
        bool operator()(const vertex_descriptor &vd) const {
            return m_map[vd] > 1;
        }

        bool operator()(const edge_descriptor &ed) const {

            auto vds = source(ed, m_tm), vdt = target(ed, m_tm);
            return (m_map[vds] > 1) &&
                   (m_map[vdt] > 1);
        }

    };


    static auto ds_component(Graph &cg1) {

        std::vector<VertexIndex> rank(boost::num_vertices(cg1));
        std::vector<vertex_descriptor> parent(boost::num_vertices(cg1));

//        decltype(rank) a = "bananana";

        boost::disjoint_sets ds(&rank[0], &parent[0]);
        boost::initialize_incremental_components(cg1, ds);
      boost::incremental_components(cg1, ds);

        return std::make_pair(ds, parent);
    }

};


int main() {
    //set graph
    Network<1,1>::Graph c(11);

    boost::add_edge(0, 1, c);
    boost::add_edge(0, 2, c);
    boost::add_edge(0, 3, c);

    boost::add_edge(6, 9, c);
    boost::add_edge(7, 9, c);
    boost::add_edge(8, 9, c);

    boost::add_edge(4, 1, c);
    boost::add_edge(4, 2, c);
    boost::add_edge(5, 2, c);
    boost::add_edge(5, 3, c);

    boost::add_edge(4, 6, c);
    boost::add_edge(4, 7, c);
    boost::add_edge(5, 7, c);
    boost::add_edge(5, 8, c);

    boost::add_edge(5, 10, c);
    //set indicator
    Network<1,1>::PMap capacity = {10, 10, 10, 10, 0, 0, 10, 10, 10, 10, 10};
    // -

    //filter it
    Network<1,1>::Is_overflown filter(capacity, c);
    boost::filtered_graph<Network<1,1>::Graph, Network<1,1>::Is_overflown, Network<1,1>::Is_overflown> fg(c, filter, filter);

    //print it
    for (auto vd: boost::make_iterator_range(boost::vertices((fg))))
        std::cout << vd << ", ";

//connected
/*    std::vector<Network::VertexIndex> rank(11);
    std::vector<Network::vertex_descriptor> parent(11);
    boost::disjoint_sets ds(&rank[0], &parent[0]);
//    int nc = boost::connected_components(fg,&comp[0]);
    boost::initialize_incremental_components(fg, ds);
    boost::incremental_components(fg, ds);*/

    Network<1,1> n;
    auto [ds1, parent] = n.ds_component(c);
//    ds1 = "banana";
    Network<1,1>::Components comp(parent.begin(), parent.end());

    std::cout << "component \n";
    typedef Network<1,1>::VertexIndex VI;
    BOOST_FOREACH( VI cd , comp) {
        std::cout << "\n " << cd << ": \n";
        BOOST_FOREACH(VI cc, comp[cd])
            std::cout << cc << ", ";
    }

    //Is it better bad as long as we only browse filtered vertices ?

    return 0;
}