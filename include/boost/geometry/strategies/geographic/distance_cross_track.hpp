// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2023-2025 Adam Wulkiewicz, Lodz, Poland.

// Copyright (c) 2016-2022, Oracle and/or its affiliates.
// Contributed and/or modified by Vissarion Fysikopoulos, on behalf of Oracle
// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_STRATEGIES_GEOGRAPHIC_DISTANCE_CROSS_TRACK_HPP
#define BOOST_GEOMETRY_STRATEGIES_GEOGRAPHIC_DISTANCE_CROSS_TRACK_HPP

#include <algorithm>

#include <boost/config.hpp>
#include <boost/concept_check.hpp>

#include <boost/geometry/core/cs.hpp>
#include <boost/geometry/core/access.hpp>
#include <boost/geometry/core/coordinate_promotion.hpp>
#include <boost/geometry/core/radian_access.hpp>
#include <boost/geometry/core/tags.hpp>

#include <boost/geometry/strategies/distance.hpp>
#include <boost/geometry/strategies/concepts/distance_concept.hpp>
#include <boost/geometry/strategies/spherical/distance_cross_track.hpp>
#include <boost/geometry/strategies/spherical/distance_haversine.hpp>
#include <boost/geometry/strategies/spherical/point_in_point.hpp>
#include <boost/geometry/strategies/geographic/azimuth.hpp>
#include <boost/geometry/strategies/geographic/distance.hpp>
#include <boost/geometry/strategies/geographic/parameters.hpp>
#include <boost/geometry/strategies/geographic/intersection.hpp>

#include <boost/geometry/formulas/vincenty_direct.hpp>

#include <boost/geometry/util/math.hpp>
#include <boost/geometry/util/select_calculation_type.hpp>
#include <boost/geometry/util/normalize_spheroidal_coordinates.hpp>

#include <boost/geometry/formulas/result_direct.hpp>
#include <boost/geometry/formulas/mean_radius.hpp>
#include <boost/geometry/formulas/spherical.hpp>

#ifdef BOOST_GEOMETRY_DEBUG_GEOGRAPHIC_CROSS_TRACK
#include <boost/geometry/io/dsv/write.hpp>
#endif

#ifndef BOOST_GEOMETRY_DETAIL_POINT_SEGMENT_DISTANCE_MAX_STEPS
#define BOOST_GEOMETRY_DETAIL_POINT_SEGMENT_DISTANCE_MAX_STEPS 100
#endif

#ifdef BOOST_GEOMETRY_DEBUG_GEOGRAPHIC_CROSS_TRACK
#include <iostream>
#endif

namespace boost { namespace geometry
{

namespace strategy { namespace distance
{
namespace detail
{

template <bool EnableClosestPoint>
struct set_result
{
    template <typename CT, typename ResultType>
    static void apply(CT const& distance,
                      CT const&,
                      CT const&,
                      ResultType& result)
    {
        result.distance = distance;
    }
};

template<>
struct set_result<true>
{
    template <typename CT, typename ResultType>
    static void apply(CT const&,
                      CT const& lon,
                      CT const& lat,
                      ResultType& result)
    {
        result.lon = lon;
        result.lat = lat;
    }
};


/*!
\brief Strategy functor for distance point to segment calculation on ellipsoid
       Algorithm uses direct and inverse geodesic problems as subroutines.
       The algorithm approximates the distance by an iterative Newton method.
\ingroup strategies
\details Class which calculates the distance of a point to a segment, for points
on the ellipsoid
\see C.F.F.Karney - Geodesics on an ellipsoid of revolution,
      https://arxiv.org/abs/1102.1215
\tparam FormulaPolicy underlying point-point distance strategy
\tparam Spheroid is the spheroidal model used
\tparam CalculationType \tparam_calculation
\tparam EnableClosestPoint computes the closest point on segment if true
*/
template
<
    typename FormulaPolicy = strategy::andoyer,
    typename Spheroid = srs::spheroid<double>,
    typename CalculationType = void,
    bool Bisection = false,
    bool EnableClosestPoint = false
>
class geographic_cross_track
{

public:

    geographic_cross_track() = default;

    explicit geographic_cross_track(Spheroid const& spheroid)
        : m_spheroid(spheroid)
    {}

    Spheroid const& model() const
    {
        return m_spheroid;
    }

    template <typename Point, typename PointOfSegment>
    struct return_type
        : promote_floating_point
          <
              typename select_calculation_type
                  <
                      Point,
                      PointOfSegment,
                      CalculationType
                  >::type
          >
    {};

    template <typename Point, typename PointOfSegment>
    auto apply(Point const& p,
               PointOfSegment const& sp1,
               PointOfSegment const& sp2) const
    {
        return apply(get_as_radian<0>(sp1), get_as_radian<1>(sp1),
                     get_as_radian<0>(sp2), get_as_radian<1>(sp2),
                     get_as_radian<0>(p), get_as_radian<1>(p),
                     m_spheroid).distance;
    }

    // points on a meridian not crossing poles
    template <typename CT>
    inline CT vertical_or_meridian(CT const& lat1, CT const& lat2) const
    {
        using meridian_inverse = typename formula::meridian_inverse
        <
            CT,
            strategy::default_order<FormulaPolicy>::value
        >;

        return meridian_inverse::meridian_not_crossing_pole_dist(lat1, lat2, m_spheroid);
    }

private:

    template <typename CT>
    struct result_type
    {
        result_type()
            : distance(0)
            , lon(0)
            , lat(0)
        {}

        CT distance;
        CT lon;
        CT lat;
    };

    template <typename CT>
    static inline CT normalize(CT const& g4, CT& der)
    {
        CT const pi = math::pi<CT>();
        CT const c_1_5 = 1.5;

        if (g4 < -1.25*pi)//close to -270
        {
#ifdef BOOST_GEOMETRY_DEBUG_GEOGRAPHIC_CROSS_TRACK
            std::cout << "g4=" << g4 * math::r2d<CT>() <<  ", close to -270" << std::endl;
#endif
            return g4 + c_1_5 * pi;
        }
        else if (g4 > 1.25*pi)//close to 270
        {
#ifdef BOOST_GEOMETRY_DEBUG_GEOGRAPHIC_CROSS_TRACK
            std::cout << "g4=" << g4 * math::r2d<CT>() <<  ", close to 270" << std::endl;
#endif
            der = -der;
            return - g4 + c_1_5 * pi;
        }
        else if (g4 < 0 && g4 > -0.75*pi)//close to -90
        {
#ifdef BOOST_GEOMETRY_DEBUG_GEOGRAPHIC_CROSS_TRACK
            std::cout << "g4=" << g4 * math::r2d<CT>() <<  ", close to -90" << std::endl;
#endif
            der = -der;
            return -g4 - pi/2;
        }
        return g4 - pi/2;
    }


    template <typename CT>
    static void bisection(CT const& lon1, CT const& lat1, //p1
                          CT const& lon2, CT const& lat2, //p2
                          CT const& lon3, CT const& lat3, //query point p3
                          Spheroid const& spheroid,
                          CT const& s14_start, CT const& a12,
                          result_type<CT>& result)
    {
        using direct_distance_type =
            typename FormulaPolicy::template direct<CT, true, false, false, false>;
        using inverse_distance_type =
            typename FormulaPolicy::template inverse<CT, true, false, false, false, false>;

        geometry::formula::result_direct<CT> res14;

        int counter = 0; // robustness
        bool dist_improve = true;

        CT pl_lon = lon1;
        CT pl_lat = lat1;
        CT pr_lon = lon2;
        CT pr_lat = lat2;

        CT s14 = s14_start;

        do {
            // Solve the direct problem to find p4 (GEO)
            res14 = direct_distance_type::apply(lon1, lat1, s14, a12, spheroid);

#ifdef BOOST_GEOMETRY_DEBUG_GEOGRAPHIC_CROSS_TRACK
            std::cout << "dist(pl,p3)="
                      << inverse_distance_type::apply(lon3, lat3, pr_lon, pr_lat, spheroid).distance
                      << std::endl;
            std::cout << "dist(pr,p3)="
                      << inverse_distance_type::apply(lon3, lat3, pr_lon, pr_lat, spheroid).distance
                      << std::endl;
#endif
            CT const dist_l =
                inverse_distance_type::apply(lon3, lat3, pl_lon, pl_lat, spheroid).distance;
            CT const dist_r =
                inverse_distance_type::apply(lon3, lat3, pr_lon, pr_lat, spheroid).distance;

            if (dist_l < dist_r)
            {
                s14 -= inverse_distance_type::apply(res14.lon2, res14.lat2, pl_lon,
                    pl_lat, spheroid).distance/2;
                pr_lon = res14.lon2;
                pr_lat = res14.lat2;
            }
            else
            {
                s14 += inverse_distance_type::apply(res14.lon2, res14.lat2, pr_lon,
                    pr_lat, spheroid).distance/2;
                pl_lon = res14.lon2;
                pl_lat = res14.lat2;
            }

            CT const new_distance = inverse_distance_type::apply(lon3, lat3, res14.lon2, res14.lat2,
                spheroid).distance;

            dist_improve = new_distance != result.distance;

#ifdef BOOST_GEOMETRY_DEBUG_GEOGRAPHIC_CROSS_TRACK
            std::cout << "p4=" << res14.lon2 * math::r2d<CT>() <<
                         "," << res14.lat2 * math::r2d<CT>() << std::endl;
            std::cout << "pl=" << pl_lon * math::r2d<CT>() << ","
                      << pl_lat * math::r2d<CT>()<< std::endl;
            std::cout << "pr=" << pr_lon * math::r2d<CT>() << ","
                      << pr_lat * math::r2d<CT>() << std::endl;
            std::cout << "new_s14=" << s14 << std::endl;
            std::cout << std::setprecision(16) << "result.distance ="
                      << result.distance << std::endl;
            std::cout << std::setprecision(16) << "new_distance ="
                      << new_distance << std::endl;
            std::cout << "---------end of step " << counter
                      << std::endl<< std::endl;
            if (!dist_improve)
            {
                std::cout << "Stop msg: res34.distance >= prev_distance" << std::endl;
            }
            if (counter == BOOST_GEOMETRY_DETAIL_POINT_SEGMENT_DISTANCE_MAX_STEPS)
            {
                std::cout << "Stop msg: counter" << std::endl;
            }
#endif
            set_result<EnableClosestPoint>::apply(new_distance, res14.lon2, res14.lat2, result);

        } while (dist_improve && counter++
                 < BOOST_GEOMETRY_DETAIL_POINT_SEGMENT_DISTANCE_MAX_STEPS);
    }

    template <typename CT>
    static void newton(CT const& lon1, CT const& lat1, //p1
                       CT const& lon2, CT const& lat2, //p2
                       CT const& lon3, CT const& lat3, //query point p3
                       Spheroid const& spheroid,
                       CT const& s14_start, CT const& a12,
                       result_type<CT>& result)
    {
        using inverse_distance_azimuth_quantities_type =
            typename FormulaPolicy::template inverse<CT, true, true, false, true, true>;

        using inverse_dist_azimuth_type =
            typename FormulaPolicy::template inverse<CT, false, true, false, false, false>;

        using direct_distance_type =
            typename FormulaPolicy::template direct<CT, true, false, false, false>;

        CT const half_pi = math::pi<CT>() / CT(2);
        geometry::formula::result_direct<CT> res14;
        geometry::formula::result_inverse<CT> res34;
        res34.distance = -1;

        int counter = 0; // robustness
        CT g4;
        CT delta_g4 = 0;
        bool dist_improve = true;
        CT s14 = s14_start;

        do {
            auto prev_distance = res34.distance;
            auto prev_res = res14;

            // Solve the direct problem to find p4 (GEO)
            res14 = direct_distance_type::apply(lon1, lat1, s14, a12, spheroid);

            // Solve an inverse problem to find g4
            // g4 is the angle between segment (p1,p2) and segment (p3,p4)
            // that meet on p4 (GEO)

            CT const a4 = inverse_dist_azimuth_type::apply(res14.lon2, res14.lat2,
                                                           lon2, lat2, spheroid).azimuth;
            res34 = inverse_distance_azimuth_quantities_type::apply(res14.lon2, res14.lat2,
                                                                    lon3, lat3, spheroid);
            g4 = res34.azimuth - a4;

            // cos(s14/earth_radius) is the spherical limit
            CT M43 = res34.geodesic_scale;
            CT m34 = res34.reduced_length;

            if (m34 != 0)
            {
                CT der = (M43 / m34) * sin(g4);
                delta_g4 = normalize(g4, der);
                s14 -= der != 0 ? delta_g4 / der : 0;
            }

            dist_improve = prev_distance > res34.distance || prev_distance == -1;

            if (dist_improve)
            {
                set_result<EnableClosestPoint>::apply(res34.distance, res14.lon2, res14.lat2,
                                                      result);
            }
            else
            {
                set_result<EnableClosestPoint>::apply(prev_distance, prev_res.lon2, prev_res.lat2,
                                                      result);
            }

#ifdef BOOST_GEOMETRY_DEBUG_GEOGRAPHIC_CROSS_TRACK
            std::cout << "p4=" << res14.lon2 * math::r2d<CT>() <<
                         "," << res14.lat2 * math::r2d<CT>() << std::endl;
            std::cout << "a34=" << res34.azimuth * math::r2d<CT>() << std::endl;
            std::cout << "a4=" << a4 * math::r2d<CT>() << std::endl;
            std::cout << "g4(normalized)=" << g4 * math::r2d<CT>() << std::endl;
            std::cout << "delta_g4=" << delta_g4 * math::r2d<CT>()  << std::endl;
            std::cout << "M43=" << M43 << std::endl;
            std::cout << "m34=" << m34 << std::endl;
            std::cout << "new_s14=" << s14 << std::endl;
            std::cout << std::setprecision(16) << "dist     ="
                      << res34.distance << std::endl;
            std::cout << "---------end of step " << counter
                      << std::endl<< std::endl;
            if (g4 == half_pi)
            {
                std::cout << "Stop msg: g4 == half_pi" << std::endl;
            }
            if (!dist_improve)
            {
                std::cout << "Stop msg: res34.distance >= prev_distance" << std::endl;
            }
            if (delta_g4 == 0)
            {
                std::cout << "Stop msg: delta_g4 == 0" << std::endl;
            }
            if (counter == BOOST_GEOMETRY_DETAIL_POINT_SEGMENT_DISTANCE_MAX_STEPS)
            {
                std::cout << "Stop msg: counter" << std::endl;
            }
#endif

        } while (g4 != half_pi
                 && dist_improve
                 && delta_g4 != 0
                 && counter++ < BOOST_GEOMETRY_DETAIL_POINT_SEGMENT_DISTANCE_MAX_STEPS);
#ifdef BOOST_GEOMETRY_DEBUG_GEOGRAPHIC_CROSS_TRACK
        std::cout << "distance=" << res34.distance << std::endl;

        std::cout << "s34(geo) ="
                  << inverse_distance_azimuth_quantities_type
                     ::apply(res14.lon2, res14.lat2,
                             lon3, lat3, spheroid).distance
                  << ", p4=(" << res14.lon2 * math::r2d<CT>() << ","
                              << res14.lat2 * math::r2d<CT>() << ")"
                  << std::endl;

        CT s31 = inverse_distance_azimuth_quantities_type::apply(lon3, lat3, lon1, lat1, spheroid)
            .distance;
        CT s32 = inverse_distance_azimuth_quantities_type::apply(lon3, lat3, lon2, lat2, spheroid)
            .distance;

        CT a4 = inverse_dist_azimuth_type::apply(res14.lon2, res14.lat2, lon2, lat2, spheroid)
            .azimuth;

        geometry::formula::result_direct<CT> res4 =
            direct_distance_type::apply(res14.lon2, res14.lat2, .04, a4, spheroid);

        CT p4_plus = inverse_distance_azimuth_quantities_type::apply(res4.lon2, res4.lat2, lon3,
            lat3, spheroid).distance;

        geometry::formula::result_direct<CT> res1 =
            direct_distance_type::apply(lon1, lat1, s14-.04, a12, spheroid);

        CT p4_minus = inverse_distance_azimuth_quantities_type::apply(res1.lon2, res1.lat2, lon3,
            lat3, spheroid).distance;

        std::cout << "s31=" << s31 << "\ns32=" << s32
                  << "\np4_plus=" << p4_plus << ", p4=("
                  << res4.lon2 * math::r2d<CT>() << ","
                  << res4.lat2 * math::r2d<CT>() << ")"
                  << "\np4_minus=" << p4_minus << ", p4=("
                  << res1.lon2 * math::r2d<CT>() << ","
                  << res1.lat2 * math::r2d<CT>() << ")"
                  << std::endl;

        if (res34.distance <= p4_plus && res34.distance <= p4_minus)
        {
            std::cout << "Closest point computed" << std::endl;
        }
        else
        {
            std::cout << "There is a closer point nearby" << std::endl;
        }
#endif
    }


    template <typename CT>
    static inline auto non_iterative_case(CT const&     , CT const&     , //p1
                                          CT const& lon2, CT const& lat2, //p2
                                          CT const& distance)
    {
        result_type<CT> result;

        set_result<EnableClosestPoint>::apply(distance, lon2, lat2, result);

        return result;
    }

    template <typename CT>
    static inline auto non_iterative_case(CT const& lon1, CT const& lat1, //p1
                                          CT const& lon2, CT const& lat2, //p2
                                          Spheroid const& spheroid)
    {
        CT distance = geometry::strategy::distance::geographic
            <
                FormulaPolicy,
                Spheroid,
                CT
            >::apply(lon1, lat1, lon2, lat2, spheroid);

        return non_iterative_case(lon1, lat1, lon2, lat2, distance);
    }

protected:

    template <typename CT>
    static inline auto apply(CT const& lo1, CT const& la1, //p1
                             CT const& lo2, CT const& la2, //p2
                             CT const& lo3, CT const& la3, //query point p3
                             Spheroid const& spheroid)
    {
        using inverse_dist_azimuth_type =
            typename FormulaPolicy::template inverse<CT, true, true, false, false, false>;

        using inverse_dist_azimuth_reverse_type =
            typename FormulaPolicy::template inverse<CT, true, true, true, false, false>;

        CT const earth_radius = geometry::formula::mean_radius<CT>(spheroid);

        result_type<CT> result;

        // if the query points coincide with one of segments' endpoints
        if ((lo1 == lo3 && la1 == la3) || (lo2 == lo3 && la2 == la3))
        {
            result.lon = lo3;
            result.lat = la3;
            return result;
        }

        // Constants
        //CT const f = geometry::formula::flattening<CT>(spheroid);
        CT const pi = math::pi<CT>();
        CT const half_pi = pi / CT(2);
        CT const c0 = CT(0);

        CT lon1 = lo1;
        CT lat1 = la1;
        CT lon2 = lo2;
        CT lat2 = la2;
        CT const lon3 = lo3;
        CT const lat3 = la3;

        if (lon1 > lon2)
        {
            std::swap(lon1, lon2);
            std::swap(lat1, lat2);
        }

#ifdef BOOST_GEOMETRY_DEBUG_GEOGRAPHIC_CROSS_TRACK
        std::cout << ">>\nSegment=(" << lon1 * math::r2d<CT>();
        std::cout << "," << lat1 * math::r2d<CT>();
        std::cout << "),(" << lon2 * math::r2d<CT>();
        std::cout << "," << lat2 * math::r2d<CT>();
        std::cout << ")\np=(" << lon3 * math::r2d<CT>();
        std::cout << "," << lat3 * math::r2d<CT>();
        std::cout << ")" << std::endl;
#endif

        //segment on equator
        //Note: antipodal points on equator does not define segment on equator
        //but pass by the pole
        CT const diff = geometry::math::longitude_distance_signed<geometry::radian>(lon1, lon2);

        using meridian_inverse = typename formula::meridian_inverse<CT>;

        bool const meridian_not_crossing_pole =
            meridian_inverse::meridian_not_crossing_pole(lat1, lat2, diff);

        bool const meridian_crossing_pole =
            meridian_inverse::meridian_crossing_pole(diff);

        if (math::equals(lat1, c0) && math::equals(lat2, c0)
            && !meridian_crossing_pole)
        {
#ifdef BOOST_GEOMETRY_DEBUG_GEOGRAPHIC_CROSS_TRACK
            std::cout << "Equatorial segment" << std::endl;
            std::cout << "segment=(" << lon1 * math::r2d<CT>();
            std::cout << "," << lat1 * math::r2d<CT>();
            std::cout << "),(" << lon2 * math::r2d<CT>();
            std::cout << "," << lat2 * math::r2d<CT>();
            std::cout << ")\np=(" << lon3 * math::r2d<CT>();
            std::cout << "," << lat3 * math::r2d<CT>() << ")\n";
#endif
            if (lon3 <= lon1)
            {
                return non_iterative_case(lon3, lat3, lon1, lat1, spheroid);
            }
            if (lon3 >= lon2)
            {
                return non_iterative_case(lon3, lat3, lon2, lat2, spheroid);
            }
            return non_iterative_case(lon3, lat3, lon3, lat1, spheroid);
        }

        if ( (meridian_not_crossing_pole || meridian_crossing_pole )
            && std::abs(lat1) > std::abs(lat2))
        {
#ifdef BOOST_GEOMETRY_DEBUG_GEOGRAPHIC_CROSS_TRACK
            std::cout << "Meridian segment not crossing pole" << std::endl;
#endif
            std::swap(lat1,lat2);
        }

        if (meridian_crossing_pole)
        {
#ifdef BOOST_GEOMETRY_DEBUG_GEOGRAPHIC_CROSS_TRACK
            std::cout << "Meridian segment crossing pole" << std::endl;
#endif
            CT const sign_non_zero = lat3 >= c0 ? 1 : -1;

            auto const res13 = apply(lon1, lat1, lon1, half_pi * sign_non_zero, lon3, lat3, spheroid);

            auto const res23 = apply(lon2, lat2, lon2, half_pi * sign_non_zero, lon3, lat3, spheroid);

            return (res13.distance) < (res23.distance) ? res13 : res23;
        }

        auto const res12 = inverse_dist_azimuth_reverse_type::apply(lon1, lat1, lon2, lat2, spheroid);

        auto const res13 = inverse_dist_azimuth_type::apply(lon1, lat1, lon3, lat3, spheroid);

        if (geometry::math::equals(res12.distance, c0))
        {
#ifdef BOOST_GEOMETRY_DEBUG_GEOGRAPHIC_CROSS_TRACK
            std::cout << "Degenerate segment" << std::endl;
            std::cout << "distance between points="
                      << res13.distance << std::endl;
#endif
            auto const res = meridian_inverse::apply(lon1, lat1, lon3, lat3, spheroid);

            return non_iterative_case(lon3, lat3, lon1, lat2,
                res.meridian ? res.distance : res13.distance);
        }

        // Compute a12 (GEO)
        CT const a312 = res13.azimuth - res12.azimuth;

        // TODO: meridian case optimization
        if (geometry::math::equals(a312, c0) && meridian_not_crossing_pole)
        {
            auto const minmax_elem = std::minmax(lat1, lat2);

            if (lat3 >= std::get<0>(minmax_elem) &&
                lat3 <= std::get<1>(minmax_elem))
            {
#ifdef BOOST_GEOMETRY_DEBUG_GEOGRAPHIC_CROSS_TRACK
                std::cout << "Point on meridian segment" << std::endl;
#endif
                return non_iterative_case(lon3, lat3, lon3, lat3, c0);
            }
        }

        CT const projection1 = cos( a312 ) * res13.distance / res12.distance;

#ifdef BOOST_GEOMETRY_DEBUG_GEOGRAPHIC_CROSS_TRACK
        std::cout << "a1=" << res12.azimuth * math::r2d<CT>() << std::endl;
        std::cout << "a13=" << res13.azimuth * math::r2d<CT>() << std::endl;
        std::cout << "a312=" << a312 * math::r2d<CT>() << std::endl;
        std::cout << "cos(a312)=" << cos(a312) << std::endl;
        std::cout << "projection 1=" << projection1 << std::endl;
#endif

        if (projection1 < c0)
        {
#ifdef BOOST_GEOMETRY_DEBUG_GEOGRAPHIC_CROSS_TRACK
            std::cout << "projection closer to p1" << std::endl;
#endif
            // projection of p3 on geodesic spanned by segment (p1,p2) fall
            // outside of segment on the side of p1

            return non_iterative_case(lon3, lat3, lon1, lat1, spheroid);
        }

        auto const res23 = inverse_dist_azimuth_type::apply(lon2, lat2, lon3, lat3, spheroid);

        CT const a321 = res23.azimuth - res12.reverse_azimuth + pi;
        CT const projection2 = cos( a321 ) * res23.distance / res12.distance;

#ifdef BOOST_GEOMETRY_DEBUG_GEOGRAPHIC_CROSS_TRACK
        std::cout << "a21=" << res12.reverse_azimuth * math::r2d<CT>()
                  << std::endl;
        std::cout << "a23=" << res23.azimuth * math::r2d<CT>() << std::endl;
        std::cout << "a321=" << a321 * math::r2d<CT>() << std::endl;
        std::cout << "cos(a321)=" << cos(a321) << std::endl;
        std::cout << "projection 2=" << projection2 << std::endl;
#endif

        if (projection2 < c0)
        {
#ifdef BOOST_GEOMETRY_DEBUG_GEOGRAPHIC_CROSS_TRACK
            std::cout << "projection closer to p2" << std::endl;
#endif
            // projection of p3 on geodesic spanned by segment (p1,p2) fall
            // outside of segment on the side of p2
            return non_iterative_case(lon3, lat3, lon2, lat2, spheroid);
        }

        // Guess s14 (SPHERICAL) aka along-track distance
        using point = geometry::model::point
            <
                CT,
                2,
                geometry::cs::spherical_equatorial<geometry::radian>
            >;

        point const p1(lon1, lat1);
        point const p2(lon2, lat2);
        point const p3(lon3, lat3);

        using haversine_t = geometry::strategy::distance::haversine<CT>;
        using cross_track_t = geometry::strategy::distance::cross_track<void, haversine_t>;

        cross_track_t const cross_track(earth_radius);
        CT const s34_sph = cross_track.apply(p3, p1, p2);

        haversine_t const str(earth_radius);
        CT const s13_sph = str.apply(p1, p3);

        //CT s14 = acos( cos(s13/earth_radius) / cos(s34/earth_radius) ) * earth_radius;
        CT const cos_frac = cos(s13_sph / earth_radius) / cos(s34_sph / earth_radius);
        CT const s14_sph = cos_frac >= 1 ? CT(0)
            : cos_frac <= -1 ? pi * earth_radius
            : acos(cos_frac) * earth_radius;

        CT const a12_sph = geometry::formula::spherical_azimuth<>(lon1, lat1, lon2, lat2);

        auto const res = geometry::formula::spherical_direct<true, false>(lon1, lat1,
            s14_sph, a12_sph, srs::sphere<CT>(earth_radius));

        // this is what postgis (version 2.5) returns
        // geometry::strategy::distance::geographic<FormulaPolicy, Spheroid, CT>
        //                     ::apply(lon3, lat3, res.lon2, res.lat2, spheroid);

#ifdef BOOST_GEOMETRY_DEBUG_GEOGRAPHIC_CROSS_TRACK
        std::cout << "s34=" << s34_sph << std::endl;
        std::cout << "s13=" << res13.distance << std::endl;
        std::cout << "s14=" << s14_sph << std::endl;
        std::cout << "===============" << std::endl;
#endif

        // Update s14 (using Newton method)

        if (Bisection)
        {
            bisection(lon1, lat1, lon2, lat2, lon3, lat3, spheroid,
                res12.distance/2, res12.azimuth, result);
        }
        else
        {
            CT const s14_start = geometry::strategy::distance::geographic
                <
                    FormulaPolicy,
                    Spheroid,
                    CT
                >::apply(lon1, lat1, res.lon2, res.lat2, spheroid);

            newton(lon1, lat1, lon2, lat2, lon3, lat3, spheroid, s14_start, res12.azimuth, result);
        }

        return result;
    }

    Spheroid m_spheroid;
};

} // namespace detail

template
<
    typename FormulaPolicy = strategy::andoyer,
    typename Spheroid = srs::spheroid<double>,
    typename CalculationType = void
>
class geographic_cross_track
    : public detail::geographic_cross_track
        <
            FormulaPolicy,
            Spheroid,
            CalculationType,
            false,
            false
        >
{
    using base_t = detail::geographic_cross_track
        <
            FormulaPolicy,
            Spheroid,
            CalculationType,
            false,
            false
        >;

public :
    explicit geographic_cross_track(Spheroid const& spheroid = Spheroid())
        : base_t(spheroid)
        {}
};

#ifndef DOXYGEN_NO_STRATEGY_SPECIALIZATIONS
namespace services
{

//tags
template <typename FormulaPolicy>
struct tag<geographic_cross_track<FormulaPolicy> >
{
    typedef strategy_tag_distance_point_segment type;
};

template
<
        typename FormulaPolicy,
        typename Spheroid
>
struct tag<geographic_cross_track<FormulaPolicy, Spheroid> >
{
    typedef strategy_tag_distance_point_segment type;
};

template
<
    typename FormulaPolicy,
    typename Spheroid,
    typename CalculationType
>
struct tag<geographic_cross_track<FormulaPolicy, Spheroid, CalculationType> >
{
    typedef strategy_tag_distance_point_segment type;
};

template
<
    typename FormulaPolicy,
    typename Spheroid,
    typename CalculationType,
    bool Bisection,
    bool EnableClosestPoint
>
struct tag<detail::geographic_cross_track<FormulaPolicy, Spheroid, CalculationType, Bisection, EnableClosestPoint> >
{
    typedef strategy_tag_distance_point_segment type;
};

//return types
template <typename FormulaPolicy, typename P, typename PS>
struct return_type<geographic_cross_track<FormulaPolicy>, P, PS>
    : geographic_cross_track<FormulaPolicy>::template return_type<P, PS>
{};

template
<
    typename FormulaPolicy,
    typename Spheroid,
    typename P,
    typename PS
>
struct return_type<geographic_cross_track<FormulaPolicy, Spheroid>, P, PS>
    : geographic_cross_track<FormulaPolicy, Spheroid>::template return_type<P, PS>
{};

template
<
    typename FormulaPolicy,
    typename Spheroid,
    typename CalculationType,
    typename P,
    typename PS
>
struct return_type<geographic_cross_track<FormulaPolicy, Spheroid, CalculationType>, P, PS>
    : geographic_cross_track<FormulaPolicy, Spheroid, CalculationType>::template return_type<P, PS>
{};

template
<
    typename FormulaPolicy,
    typename Spheroid,
    typename CalculationType,
    bool Bisection,
    bool EnableClosestPoint,
    typename P,
    typename PS
>
struct return_type<detail::geographic_cross_track<FormulaPolicy, Spheroid, CalculationType, Bisection, EnableClosestPoint>, P, PS>
    : detail::geographic_cross_track<FormulaPolicy, Spheroid, CalculationType, Bisection, EnableClosestPoint>::template return_type<P, PS>
{};

//comparable types
template
<
    typename FormulaPolicy,
    typename Spheroid,
    typename CalculationType
>
struct comparable_type<geographic_cross_track<FormulaPolicy, Spheroid, CalculationType> >
{
    typedef geographic_cross_track
        <
            FormulaPolicy, Spheroid, CalculationType
        >  type;
};


template
<
    typename FormulaPolicy,
    typename Spheroid,
    typename CalculationType,
    bool Bisection,
    bool EnableClosestPoint
>
struct comparable_type<detail::geographic_cross_track<FormulaPolicy, Spheroid, CalculationType, Bisection, EnableClosestPoint> >
{
    typedef detail::geographic_cross_track
        <
            FormulaPolicy, Spheroid, CalculationType, Bisection, EnableClosestPoint
        >  type;
};

template
<
    typename FormulaPolicy,
    typename Spheroid,
    typename CalculationType
>
struct get_comparable<geographic_cross_track<FormulaPolicy, Spheroid, CalculationType> >
{
public :
    static inline geographic_cross_track<FormulaPolicy, Spheroid, CalculationType>
    apply(geographic_cross_track<FormulaPolicy, Spheroid, CalculationType> const& strategy)
    {
        return strategy;
    }
};

template
<
    typename FormulaPolicy,
    typename Spheroid,
    typename CalculationType,
    bool Bisection,
    bool EnableClosestPoint
>
struct get_comparable<detail::geographic_cross_track<FormulaPolicy, Spheroid, CalculationType, Bisection, EnableClosestPoint> >
{
public :
    static inline detail::geographic_cross_track<FormulaPolicy, Spheroid, CalculationType, Bisection, EnableClosestPoint>
    apply(detail::geographic_cross_track<FormulaPolicy, Spheroid, CalculationType, Bisection, EnableClosestPoint> const& strategy)
    {
        return strategy;
    }
};


template
<
    typename FormulaPolicy,
    typename Spheroid,
    typename CalculationType,
    typename P,
    typename PS
>
struct result_from_distance<geographic_cross_track<FormulaPolicy, Spheroid, CalculationType>, P, PS>
{
private :
    typedef typename geographic_cross_track
        <
            FormulaPolicy, Spheroid, CalculationType
        >::template return_type<P, PS>::type return_type;
public :
    template <typename T>
    static inline return_type
    apply(geographic_cross_track<FormulaPolicy, Spheroid, CalculationType> const& , T const& distance)
    {
        return distance;
    }
};

template
<
    typename FormulaPolicy,
    typename Spheroid,
    typename CalculationType,
    bool Bisection,
    bool EnableClosestPoint,
    typename P,
    typename PS
>
struct result_from_distance<detail::geographic_cross_track<FormulaPolicy, Spheroid, CalculationType, Bisection, EnableClosestPoint>, P, PS>
{
private :
    typedef typename detail::geographic_cross_track
        <
            FormulaPolicy, Spheroid, CalculationType, Bisection, EnableClosestPoint
        >::template return_type<P, PS>::type return_type;
public :
    template <typename T>
    static inline return_type
    apply(detail::geographic_cross_track<FormulaPolicy, Spheroid, CalculationType, Bisection, EnableClosestPoint> const& , T const& distance)
    {
        return distance;
    }
};


template <typename Point, typename PointOfSegment>
struct default_strategy
    <
        point_tag, segment_tag, Point, PointOfSegment,
        geographic_tag, geographic_tag
    >
{
    typedef geographic_cross_track<> type;
};


template <typename PointOfSegment, typename Point>
struct default_strategy
    <
        segment_tag, point_tag, PointOfSegment, Point,
        geographic_tag, geographic_tag
    >
{
    typedef typename default_strategy
        <
            point_tag, segment_tag, Point, PointOfSegment,
            geographic_tag, geographic_tag
        >::type type;
};

} // namespace services
#endif // DOXYGEN_NO_STRATEGY_SPECIALIZATIONS

}} // namespace strategy::distance

}} // namespace boost::geometry
#endif // BOOST_GEOMETRY_STRATEGIES_GEOGRAPHIC_DISTANCE_CROSS_TRACK_HPP
