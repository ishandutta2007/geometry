// Boost.Geometry

// Copyright (c) 2018-2019 Barend Gehrels, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_ARITHMETIC_GENERAL_FORM_PREDICATES_HPP
#define BOOST_GEOMETRY_ARITHMETIC_GENERAL_FORM_PREDICATES_HPP

#include <boost/geometry/util/math.hpp>
#include <boost/geometry/core/access.hpp>
#include <boost/geometry/core/config.hpp>
#include <boost/geometry/geometries/infinite_line.hpp>
#include <boost/geometry/util/select_most_precise.hpp>
#include <boost/geometry/strategies/cartesian/detail/thresholds.hpp>

namespace boost { namespace geometry
{

namespace arithmetic
{

// Structure containing thresholds, per type, for denominator.
// Determined with corresponding unit test.
// It should not be replaced by machine epsilon or math::equals
template <typename Type>
struct general_threshold
{
    static Type get_cl() { return 0; } // used for collinearity
};

template <>
struct general_threshold<long double>
{
   static long double get_cl() { return 1.0e-4; } // used for collinearity
};

template <>
struct general_threshold<double>
{
   static double get_cl() { return 1.0e-4; } // used for collinearity
};

template <>
struct general_threshold<float>
{
   static float get_cl() { return 1.0e-4; }
};

// Calculates intersection point (in homogeneous coordinates)
// of two infinite lines. The so-called z-value (denominator) is also returned.
// It does not use the point type, because coordinate type might be promoted
template <typename Type>
inline
void get_intersection(model::infinite_line<Type> const& p,
                      model::infinite_line<Type> const& q,
                      Type& x, Type& y, Type& homogeneous_z)
{
    // Calculate x and y, and the z-value in homogeneous coordinates
    // Do not divide it yet: the undivided homogeneous coordinates can be
    // used to calculate the fraction, without division it might be more precise
    x = p.c * q.b - p.b * q.c;
    y = p.a * q.c - p.c * q.a;
    homogeneous_z =      p.b * q.a - p.a * q.b;
}

//! Normalize the line. For robustness reasons it is often better to
//! NOT use normalization. It uses sqrt and therefore the intersection
//! point, if calculated with normalization, might go,
//! for example, from 7 to 6.99...997
template <typename FloatingPointType, typename InputType>
inline
model::infinite_line<FloatingPointType> normalize_line(model::infinite_line<InputType> const& line)
{
    model::infinite_line<FloatingPointType> result;
    result.a = line.a;
    result.b = line.b;
    result.c = line.c;

    FloatingPointType const norm
            = std::sqrt(result.a * result.a + result.b * result.b);

    // Compare with 0 (even very small values like 1.0e-12 are supported)
    if (norm != 0)
    {
        result.a /= norm;
        result.b /= norm;
        result.c /= norm;
        result.normalized = true;
    }
    return result;
}

template <typename Type>
inline bool more_horizontal(model::infinite_line<Type> const& line)
{
    // If a=0, then line is 'by+c=0' so y=-c/b  which is horizontal
    // If a=0.1 and b=0.9, then the line is quite horizontal
    return std::abs(line.a) < std::abs(line.b);
}

// Returns a side
template <typename Type, typename CoordinateType>
inline int get_side(model::infinite_line<Type> const& line,
    CoordinateType const& x, CoordinateType const& y)
{
    // https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line#Line_defined_by_an_equation
    // Distance from point to line in general form is given as:
    // (a * x + b * y + c) / sqrt(a * a + b * b);
    // In most use cases comparisons are enough, saving the sqrt
    // (performance plus making it a bit more precise)
    typedef typename select_most_precise<Type, CoordinateType>::type ct;

    ct const num = line.a * x + line.b * y + line.c;
    if (num == 0)
    {
        // Point is located on the line
        return 0;
    }

    // TODO: can be precalculated
    Type const denom = line.a * line.a + line.b * line.b;
    if (denom == 0)
    {
        // Line is degenerate
        return 0;
    }

    if (! boost::is_integral<ct>::value)
    {
        // Floating point: compare with a threshold.
        if (std::fabs(num / denom) < strategy::intersection::general_distance_threshold<ct>::get())
        {
            return 0;
        }
#if 0
        // This is wrong for points lying close together!
        // But still fixes one of the diff.bugs
        ct const threshold = denom * strategy::intersection::general_distance_threshold<ct>::get();
        if (g_side_debug)
        {
            std::cout << threshold << " " << num << std::endl;
        }

        if (std::fabs(num * num) < threshold)
        {
            std::cout << threshold << " " << num
                      << " " << num * num << " " << denom
                      << " " << num / denom
                      << (collinear ? " col" : "")
                      << std::endl;
            return 0;
        }
#endif
    }

    // For integer coordinates or distances > threshold, just compare signs
    int const num_sign  = num > 0 ? 1 : -1;
    int const denom_sign = denom > 0 ? 1 : -1;
    return num_sign == denom_sign ? 1 : -1;
}

template <typename Type>
inline
bool lines_collinear(model::infinite_line<Type> const& p, model::infinite_line<Type> const& q)
{
    if (p.normalized && q.normalized)
    {
        bool const same_sign = more_horizontal(p)
                ? p.b * q.b > 0
                : p.a * q.a > 0;

        // c is the interception on x or y axis of normalized line
        // The normalized lign is still directed, if they have the same
        // direction (same_sign), check for intercept. If they are opposite,
        // then reverse one intercept
        Type const diff = same_sign ? p.c - q.c : p.c + q.c;
        Type const threshold = general_threshold<Type>::get_cl();
        return std::fabs(diff) < threshold;

    }

    // Non normalized lines are not (yet) implemented
    return false;
}

} // namespace arithmetic


}} // namespace boost::geometry


#endif // BOOST_GEOMETRY_ARITHMETIC_GENERAL_FORM_HPP
