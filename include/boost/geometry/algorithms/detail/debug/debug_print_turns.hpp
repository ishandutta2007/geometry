// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2014-2015, Oracle and/or its affiliates.

// Contributed and/or modified by Menelaos Karavelas, on behalf of Oracle

// Licensed under the Boost Software License version 1.0.
// http://www.boost.org/users/license.html

#ifndef BOOST_GEOMETRY_ALGORITHMS_DETAIL_DEBUG_PRINT_TURNS_HPP
#define BOOST_GEOMETRY_ALGORITHMS_DETAIL_DEBUG_PRINT_TURNS_HPP

#include <boost/geometry/core/config.hpp>

#if defined(BOOST_GEOMETRY_TEST_DEBUG) || defined(BAREND_EXTRA_DEBUGGING)
#include <iostream>

#include <boost/geometry/io/dsv/write.hpp>
#include <boost/geometry/algorithms/detail/overlay/debug_turn_info.hpp>
#endif


namespace boost { namespace geometry { namespace detail { namespace debug
{

#if defined(BOOST_GEOMETRY_TEST_DEBUG) || defined(BAREND_EXTRA_DEBUGGING)
template <typename Turn>
inline void debug_print_turn(Turn const& turn)
{
    std::cout << " ["
              << method_char(turn.method)
              << ", " << operation_char(turn.operations[0].operation)
            << " " << turn.operations[0].seg_id
           // << " " << turn.operations[0].fraction
              << " // "
              << operation_char(turn.operations[1].operation)
            << " " << turn.operations[1].seg_id
            //<< " " << turn.operations[1].fraction
              << "} pnt "
              << geometry::dsv(turn.point)
              << (turn.discarded ? " discarded" : "")
              << (turn.cluster_id > 0 ? " cluster:" + std::to_string(turn.cluster_id) : "")
              << "]";
}

template <typename TurnIterator>
inline void debug_print_turns(TurnIterator first, TurnIterator beyond)
{
    std::cout << "turns:";
    for (TurnIterator tit = first; tit != beyond; ++tit)
    {
        debug_print_turn(*tit);
    }
    std::cout << std::endl << std::endl;
}
#else
template <typename Turn>
inline void debug_print_turn(Turn const&)
{}

template <typename TurnIterator>
inline void debug_print_turns(TurnIterator, TurnIterator)
{}
#endif // BOOST_GEOMETRY_TEST_DEBUG

}}}} // namespace boost::geometry::detail::debug

#endif // BOOST_GEOMETRY_ALGORITHMS_DETAIL_DEBUG_PRINT_TURNS_HPP
