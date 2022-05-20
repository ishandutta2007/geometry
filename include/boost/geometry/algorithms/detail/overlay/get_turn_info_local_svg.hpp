
#ifndef BOOST_GEOMETRY_TEST_DEBUG_GET_TURN_INFO_SVG_HPP
#define BOOST_GEOMETRY_TEST_DEBUG_GET_TURN_INFO_SVG_HPP

#include <fstream>
#include <sstream>

#include <boost/geometry/core/config.hpp>
#include <boost/geometry/io/svg/svg_mapper.hpp>

namespace boost { namespace geometry { namespace debug
{

inline int get_counter()
{
    static int counter = 0;
    return counter++;
}

template <typename Range1, typename Range2, typename Point>
inline void get_turn_info_svg(Range1 const& range1, Range2 const& range2, Point const& point)
{
    std::ostringstream filename;
#ifdef BAREND_TEST_ONLY_ONE_WRONG_CASE
    const std::string suffix = "wrong_";
#else
    const std::string suffix = "";
#endif
    filename << "/tmp/svg/get_turn_info_" << suffix << get_counter() << ".svg";
    std::ofstream svg(filename.str().c_str());

    typedef geometry::svg_mapper<Point> mapper_type;
    typedef geometry::model::referring_segment<Point const> seg;

    mapper_type mapper(svg, 500, 500);

    mapper.add(point);
    for (std::size_t i = 0; i < 3; i++)
    {
        mapper.add(range1.at(i));
        mapper.add(range2.at(i));
    }



    for (std::size_t i = 0; i < 2; i++)
    {
        mapper.map(seg(range1.at(i), range1.at(i + 1)), "opacity:0.7;stroke:rgb(153,204,0);stroke-width:2;");
        mapper.map(seg(range2.at(i), range2.at(i + 1)), "opacity:0.7;stroke:rgb(51,51,153);stroke-width:2;");
    }

    for (std::size_t i = 0; i < 3; i++)
    {
        mapper.map(range1.at(i), "opacity:0.7;fill:rgb(153,204,0);", i == 2 ? 5 : 3);
        mapper.map(range2.at(i),  "opacity:0.7;fill:rgb(51,51,153);", i == 2 ? 5 : 3);
    }
    mapper.map(point, "opacity:0.7;fill:rgb(255,0,0);", 5);
}

}}} // namespace boost::geometry::debug

#endif // BOOST_GEOMETRY_TEST_DEBUG_GET_TURN_INFO_SVG_HPP
