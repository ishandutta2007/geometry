// Boost.Geometry (aka GGL, Generic Geometry Library)
// Unit Test

// Copyright (c) 2007-2012 Barend Gehrels, Amsterdam, the Netherlands.

// This file was modified by Oracle on 2013.
// Modifications copyright (c) 2013, Oracle and/or its affiliates.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_TEST_TOUCHES_HPP
#define BOOST_GEOMETRY_TEST_TOUCHES_HPP


#include <geometry_test_common.hpp>

#include <boost/geometry/core/ring_type.hpp>
#include <boost/geometry/algorithms/touches.hpp>
#include <boost/geometry/strategies/strategies.hpp>
#include <boost/geometry/geometries/ring.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/linestring.hpp>

#include <boost/geometry/io/wkt/read.hpp>

#include <boost/geometry/multi/core/point_order.hpp>
#include <boost/geometry/multi/core/geometry_id.hpp>
#include <boost/geometry/multi/algorithms/detail/for_each_range.hpp>
#include <boost/geometry/multi/algorithms/detail/within/point_in_geometry.hpp>

#include <boost/geometry/multi/geometries/multi_linestring.hpp>
#include <boost/geometry/multi/geometries/multi_polygon.hpp>

#include <boost/geometry/multi/io/wkt/read.hpp>

template <typename Geometry1, typename Geometry2>
void test_touches(std::string const& wkt1,
        std::string const& wkt2, bool expected)
{
    Geometry1 geometry1;
    Geometry2 geometry2;

    bg::read_wkt(wkt1, geometry1);
    bg::read_wkt(wkt2, geometry2);

    bool detected = bg::touches(geometry1, geometry2);

    BOOST_CHECK_MESSAGE(detected == expected,
        "touches: " << wkt1
        << " with " << wkt2
        << " -> Expected: " << expected
        << " detected: " << detected);

    detected = bg::touches(geometry2, geometry1);

    BOOST_CHECK_MESSAGE(detected == expected,
        "touches: " << wkt2
        << " with " << wkt1
        << " -> Expected: " << expected
        << " detected: " << detected);
}


template <typename Geometry>
void test_self_touches(std::string const& wkt, bool expected)
{
    Geometry geometry;

    bg::read_wkt(wkt, geometry);

    bool detected = bg::touches(geometry);

    BOOST_CHECK_MESSAGE(detected == expected,
        "touches: " << wkt
        << " -> Expected: " << expected
        << " detected: " << detected);
}



#endif
