#ifndef AUDI_INVERT_MAP_HPP
#define AUDI_INVERT_MAP_HPP

#include <boost/numeric/ublas/matrix.hpp>

#include <audi.hpp>

namespace audi
{

std::vector<gdual_d> invert_map(const std::vector<gdual_d> map_in)
{
    // We check that the map contains N elements having symbol size of N
    auto map_degree = map_in[0].degree();
    std::vector<std::vector<gdual_d>> map_terms;
    for (const auto &map : map_in) {
        std::vector<gdual_d> map_d;
        for (decltype(map_degree) i = 0; i < map_degree; ++i) {
            map_d.push_back();
        }
    }
}

} // end of namespace audi

#endif
