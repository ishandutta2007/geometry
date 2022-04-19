// Boost.Geometry (aka GGL, Generic Geometry Library)
// Robustness Test

// Copyright (c) 2022 Barend Gehrels, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_GEOMETRY_NO_BOOST_TEST

#include <chrono>
#include <fstream>
#include <sstream>
#include <random>

#include <boost/program_options.hpp>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/geometries/point_xy.hpp>

#include <boost/geometry/strategy/cartesian/side_non_robust.hpp>
#include <boost/geometry/strategies/cartesian/side_rounded_input.hpp>

#include <geometry_test_common.hpp>
#include <robustness/common/common_settings.hpp>

struct test_settings : public common_settings
{
    long long seed{static_cast<long long>(std::time(0))};
    std::size_t count{1};
    long double size{10};
    long double epsilon_multiplier{1.0e15};
};

namespace bg = boost::geometry;

template <typename Point, typename Segment>
void create_svg(std::string const& filename,
    std::vector<Point> const& v1, std::vector<Point> const& v2,
    Segment const& sp, Segment const& sq, Point const& e)
{
    std::ofstream svg(filename);
    boost::geometry::svg_mapper<Point> mapper(svg, 800, 800, 0.95);

    for (const auto& p : v1) { mapper.add(p); }
    for (const auto& p : v2) { mapper.add(p); }
    mapper.add(e);

    mapper.map(sp, "stroke:rgb(255,128,0);stroke-width:1", 3);
    mapper.map(sq, "stroke:rgb(0,128,0);stroke-width:1", 3);
    for (const auto& p : v1) {
        mapper.map(p, "fill:rgb(255,128,0);stroke:rgb(0,0,0);stroke-width:0.5", 3);
    }
    for (const auto& p : v2) {
        mapper.map(p, "fill:rgb(0,128,0);stroke:rgb(0,0,0);stroke-width:0.5", 3);
    }
    mapper.map(e, "fill:rgb(0,0,255);stroke:rgb(0,0,0);stroke-width:0.5", 3);
}

template <typename P, typename Distribution >
P random_point(std::mt19937& gen, Distribution dis)
{
    P result;
    boost::geometry::set<0>(result, dis(gen));
    boost::geometry::set<1>(result, dis(gen));
    return result;
}

template <typename SideStrategy, typename Point>
bool verify_collinear(Point const& p1, Point const& p2, Point const& q1, Point const& q2, Point const& ip)
{
    int const side_p = SideStrategy::apply(p1, p2, ip);
    int const side_q = SideStrategy::apply(q1, q2, ip);
    return side_p == 0 && side_q == 0;
}

template <typename SideStrategy, typename Point>
bool verify_not_collinear(Point const& p1, Point const& p2, Point const& q1, Point const& q2, Point const& ip)
{
    int const side_p = SideStrategy::apply(p1, p2, ip);
    int const side_q = SideStrategy::apply(q1, q2, ip);
    return side_p != 0 && side_q != 0;
}

template <typename Point>
Point generate_bisecting_point(Point const& p1, Point const& p2, Point const& q1, Point const& q2, Point const& ip,
                               typename bg::coordinate_type<Point>::type epsilon_multiplier)
{
    using T = typename bg::coordinate_type<Point>::type;
    auto const dp1 = bg::distance(ip, p1);
    auto const dp2 = bg::distance(ip, p2);
    auto const dq1 = bg::distance(ip, q1);
    auto const dq2 = bg::distance(ip, q2);

    // Calculate the azimuth (0 = North) and modify it to Mathematical (0 = east)
    auto const azimuth_p = bg::azimuth(ip, dp1 < dp2 ? p1 : p2);// - bg::math::pi<T>() / 2.0;
    auto const azimuth_q = bg::azimuth(ip, dq1 < dq2 ? q1 : q2);// - bg::math::pi<T>() / 2.0;

    // Calculate the average circular azimuth
    auto const sum_sin = std::sin(azimuth_p) + std::sin(azimuth_q);
    auto const sum_cos = std::cos(azimuth_p) + std::cos(azimuth_q);
    auto const azimuth_avg = std::atan2(sum_sin, sum_cos);

    auto const d = std::numeric_limits<T>::epsilon() * epsilon_multiplier;

    // Calculate the point using sin,cos in this order (because North=0)
    Point result;
    bg::set<0>(result, bg::get<0>(ip) + d * std::sin(azimuth_avg));
    bg::set<1>(result, bg::get<1>(ip) + d * std::cos(azimuth_avg));

    // T const rd = 180 / bg::math::pi<T>();
    // std::cout << "Index " << azimuth_p * rd << " " << azimuth_q * rd << " -> " << azimuth_avg * rd << " dist=" << bg::distance(e, ip) << " d=" << d << std::endl;

    return result;
}

template <typename T, typename Settings>
void test_all(Settings const& settings)
{
    using point_t = bg::model::d2::point_xy<T>;
    using segment_type = bg::model::referring_segment<const point_t>;
    using result_type = bg::segment_intersection_points<point_t> ;

    using points_policy_type = bg::policies::relate::segments_intersection_points
        <
            result_type
        >;

    using side_by_triangle = bg::strategy::side::side_by_triangle<T>;
    using side_non_robust = bg::strategy::side::side_non_robust<T>;
    using side_robust_fp_3 = bg::strategy::side::side_robust<T, bg::strategy::side::fp_equals_policy, 3>;
    using side_robust_eq_3 = bg::strategy::side::side_robust<T, bg::strategy::side::epsilon_equals_policy, 3>;
    using side_rounded_input = bg::strategy::side::side_rounded_input<T>;

    auto const t0 = std::chrono::high_resolution_clock::now();

    std::size_t count = 0;
    std::size_t errors_triangle = 0;
    std::size_t errors_non_robust = 0;
    std::size_t errors_robust_fp_3 = 0;
    std::size_t errors_robust_eq_3 = 0;
    std::size_t errors_rounded = 0;
    std::size_t errors_bsp_triangle = 0;
    std::size_t errors_bsp_non_robust = 0;
    std::size_t errors_bsp_robust_fp_3 = 0;
    std::size_t errors_bsp_robust_eq_3 = 0;
    std::size_t errors_bsp_rounded = 0;

    std::mt19937 gen(settings.seed);
    std::uniform_real_distribution<T> dis(0, settings.size);

    T sum_distance_error = 0;
    T sum_distance_for_rounded_error = 0;

    for (std::size_t index = 0; index < settings.count; index++)
    {
        const point_t p1 = random_point<point_t>(gen, dis);
        const point_t p2 = random_point<point_t>(gen, dis);
        const point_t q1 = random_point<point_t>(gen, dis);
        const point_t q2 = random_point<point_t>(gen, dis);

        segment_type p(p1, p2);
        segment_type q(q1, q2);

        bg::detail::segment_as_subrange<segment_type> p_sub(p);
        bg::detail::segment_as_subrange<segment_type> q_sub(q);

        // Get the intersection point (or two points)
        const result_type is
            = bg::strategy::intersection::cartesian_segments<>
                ::apply(p_sub, q_sub, points_policy_type());

        if (is.count != 1) {
            continue;
        }

        count++;

        auto const& ip = is.intersections[0];

        bool const side_triangle_ok = verify_collinear<side_by_triangle>(p1, p2, q1, q2, ip);
        bool const side_non_robust_ok = verify_collinear<side_non_robust>(p1, p2, q1, q2, ip);
        bool const side_robust_fp_3_ok = verify_collinear<side_robust_fp_3>(p1, p2, q1, q2, ip);
        bool const side_robust_eq_3_ok = verify_collinear<side_robust_eq_3>(p1, p2, q1, q2, ip);
        bool const side_rounded_ok = verify_collinear<side_rounded_input>(p1, p2, q1, q2, ip);

        // Calculate the distance of the IP w.r.t. both lines
        auto const dm1 = bg::detail_dispatch::get_distance_measure<T, bg::cartesian_tag>::apply(p1, p2, ip);
        auto const dm2 = bg::detail_dispatch::get_distance_measure<T, bg::cartesian_tag>::apply(q1, q2, ip);
        auto const dm_distance = std::fabs(dm1.measure) + std::fabs(dm2.measure);

        sum_distance_error += dm_distance;
        if (! side_rounded_ok) {
            sum_distance_for_rounded_error += dm_distance;
        }

        if (side_triangle_ok
            && side_non_robust_ok
            && side_robust_fp_3_ok
            && side_robust_eq_3_ok
            && side_rounded_ok )
        {
            continue;
        }

        // Generate a point between the two lines and verify if it is NOT set as collinear
        auto const bsp = generate_bisecting_point(p1, p2, q1, q2, ip, settings.epsilon_multiplier);

        bool const side_triangle_bsp = verify_not_collinear<side_by_triangle>(p1, p2, q1, q2, bsp);
        bool const side_non_robust_bsp = verify_not_collinear<side_non_robust>(p1, p2, q1, q2, bsp);
        bool const side_robust_fp_3_bsp = verify_not_collinear<side_robust_fp_3>(p1, p2, q1, q2, bsp);
        bool const side_robust_eq_3_bsp = verify_not_collinear<side_robust_eq_3>(p1, p2, q1, q2, bsp);
        bool const side_rounded_bsp = verify_not_collinear<side_rounded_input>(p1, p2, q1, q2, bsp);

        if (! side_triangle_ok) { errors_triangle++; }
        if (! side_non_robust_ok) { errors_non_robust++; }
        if (! side_robust_fp_3_ok) { errors_robust_fp_3++; }
        if (! side_robust_eq_3_ok) { errors_robust_eq_3++; }
        if (! side_rounded_ok) { errors_rounded++; }

        if (! side_triangle_bsp) { errors_bsp_triangle++; }
        if (! side_non_robust_bsp) { errors_bsp_non_robust++; }
        if (! side_robust_fp_3_bsp) { errors_bsp_robust_fp_3++; }
        if (! side_robust_eq_3_bsp) { errors_bsp_robust_eq_3++; }
        if (! side_rounded_bsp) { errors_bsp_rounded++; }

        std::string filename;

        {
            std::ostringstream out;
            out << "/tmp/svg/ipas_" << index
                << "_" << string_from_type<T>::name()
                << ".";

            filename = out.str();
        }

        if (settings.svg)
        {
            create_svg(filename + "svg", std::vector<point_t>{p1, p2}, std::vector<point_t>{q1, q2}, p, q, bsp);
        }
        if (settings.wkt)
        {
            std::ofstream wkt(filename + "wkt");
        }
    }

    auto const t = std::chrono::high_resolution_clock::now();
    auto const elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t - t0).count();

    auto const avg = count > 0 ? sum_distance_error / (count * 2) : 0;
    auto const avg_rounded = count > 0 ? sum_distance_for_rounded_error / (count * 2) : 0;

    std::cout
        << "Results"
        << "  type: " << string_from_type<T>::name() << std::endl
        << "  intersections: " << count << std::endl
        << "  errors (triangle): " << errors_triangle << " " << errors_bsp_triangle << std::endl
        << "  errors (non robust): " << errors_non_robust << " " << errors_bsp_non_robust << std::endl
        << "  errors (side robust fp 3): " << errors_robust_fp_3 << " " << errors_bsp_robust_fp_3 << std::endl
        << "  errors (side robust eq 3): " << errors_robust_eq_3 << " " << errors_bsp_robust_eq_3 << std::endl
        << "  errors (rounded): " << errors_rounded << " " << errors_bsp_rounded << std::endl
        << "  distance (avg): " << std::setprecision(32) << avg
        << "  distance (avg, rnd): " << avg_rounded
        << "  time: " << std::setprecision(5) << elapsed_ms / 1000.0 << std::endl;
}

int main(int argc, char** argv)
{
    BoostGeometryWriteTestConfiguration();
    try
    {
        namespace po = boost::program_options;
        po::options_description description("=== gridded_point_buffer ===\nAllowed options");

        std::string type = "double";

        test_settings settings;

        description.add_options()
            ("help", "Help message")
            ("seed", po::value<decltype(settings.seed)>(&settings.seed), "Initialization seed for random generator")
            ("count", po::value<decltype(settings.count)>(&settings.count)->default_value(1), "Number of tests")
            ("size", po::value<decltype(settings.size)>(&settings.size)->default_value(10), "Size of the field")
            ("multiplier", po::value<decltype(settings.epsilon_multiplier)>(&settings.epsilon_multiplier), "Epsilon multiplier")
            ("type", po::value<std::string>(&type)->default_value("double"), "Type (int,float,double)")
            ("wkt", po::value<bool>(&settings.wkt)->default_value(false), "Create a WKT of the inputs, for all tests")
            ("svg", po::value<bool>(&settings.svg)->default_value(false), "Create a SVG for all tests")
        ;

        po::variables_map varmap;
        po::store(po::parse_command_line(argc, argv, description), varmap);
        po::notify(varmap);

        if (varmap.count("help") || settings.size <= 0)
        {
            std::cout << description << std::endl;
            return 1;
        }

        if (type == "float" || type == "f")
        {
            test_all<float>(settings);
        }
        else if (type == "double" || type == "d")
        {
            test_all<double>(settings);
        }
        else if (type == "long double" || type == "extended" || type == "e")
        {
            test_all<long double>(settings);
        }
        else
        {
            std::cout << description << std::endl;
            return 1;
        }
    }
    catch(std::exception const& e)
    {
        std::cout << "Exception " << e.what() << std::endl;
    }
    catch(...)
    {
        std::cout << "Other exception" << std::endl;
    }

    return 0;
}

