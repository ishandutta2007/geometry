// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2021 Tinko Bartels, Berlin, Germany.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_STRATEGY_CARTESIAN_SIDE_ROUNDED_INPUT_HPP
#define BOOST_GEOMETRY_STRATEGY_CARTESIAN_SIDE_ROUNDED_INPUT_HPP

#include <boost/geometry/core/config.hpp>

#include <boost/geometry/strategies/side.hpp>


namespace boost { namespace geometry
{

namespace strategy { namespace side
{


template <typename CalculationType, int Coeff1 = 5, int Coeff2 = 32>
struct side_rounded_input
{
    using cs_tag = cartesian_tag;
    template <typename P1, typename P2, typename P>
    static inline int apply(P1 const& p1, P2 const& p2, P const& p)
    {

        using CT = std::conditional_t
            <
                std::is_void<CalculationType>::value,
                double,
                CalculationType
            >;

        using namespace boost::geometry;
        CT _1 = get<0>(p1);
        CT _2 = get<1>(p1);
        CT _3 = get<0>(p2);
        CT _4 = get<1>(p2);
        CT _5 = get<0>(p);
        CT _6 = get<1>(p);

        auto const det = (_1 - _5) * (_4 - _6) - (_2 - _6) * (_3 - _5);
        CT const two = 2;
        CT const eps = std::numeric_limits<CT>::epsilon() / two;
        CT const err_bound = (Coeff1 * eps + Coeff2 * eps * eps) *
            (  (std::abs(_1) + std::abs(_5)) * (std::abs(_4) + std::abs(_6))
             + (std::abs(_3) + std::abs(_5)) * (std::abs(_2) + std::abs(_6)));
        return (det > err_bound) - (det < -err_bound);
    }
};

}} // namespace strategy::side

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_STRATEGY_CARTESIAN_SIDE_ROUNDED_INPUT_HPP
