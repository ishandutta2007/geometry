// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2022 Tinko Bartels, Berlin, Germany.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <array>
#include <cassert>
#include <limits>
#include <numeric>

#include <boost/geometry/util/precise_math.hpp>

namespace boost { namespace geometry
{
namespace detail { namespace precise_math
{

template <typename FPT>
FPT fast_two_sum_tail(FPT const& a, FPT const& b, FPT const& x)
{
    const FPT b_virtual = x - a;
    return b - b_virtual;
}

// Algorithm 7 in https://doi.org/10.1109/TC.2015.2441714 with m = 2
// available at https://hal.archives-ouvertes.fr/hal-01111551v2
// Normalizes a nonoverlapping (in the sense of https://doi.org/10.1007/PL00009321)
// expansion e = e[0] + ... + e[e_size - 1] to an expansion of equal value
// f = f[0] + ... + f[m], such that f[j + 1] <= ulp(f[j]), truncated at to m = 2.
template <typename FPT, std::size_t ESize>
std::array<FPT, 2> vec_sum_err_branch_2(const std::array<FPT, ESize>& e, int e_size = ESize)
{
    int j = 0;
    FPT ee = e[e_size - 1];
    std::array<FPT, 2> f {0, 0};
    for(int i = e_size - 1; i >= 1; i--)
    {
        f[j] = ee + e[i - 1];
        ee = fast_two_sum_tail(ee, e[i - 1], f[j]);
        if(ee != 0)
        {
            if(j >= 1)
            {
                assert(std::abs(f[1]) <= std::abs(f[0]) * std::numeric_limits<FPT>::epsilon());
                return f;
            }
            ++j;
        }
        else
        {
            ee = f[j];
        }
    }
    if (ee != 0 && j < 2)
        f[j] = ee;
    assert(std::abs(f[1]) <= std::abs(f[0]) * std::numeric_limits<FPT>::epsilon());
    return f;
}

// Algorithm 9 in https://doi.org/10.1145/3121432
// Computes an approximation z ~= (xh + xl) * y such that the relative error is less than
// or equal to 2 * u ^ 2 with // u = std::numeric_limits<FPT>::epsilon() / 2.
template <typename FPT>
std::array<FPT, 2> dwtimesfp3(FPT xh, FPT xl, FPT y)
{
    const auto c1 = two_product(xh, y);
    const FPT cl3 = std::fma(xl, y, c1[1]);
    const auto z = fast_two_sum(c1[0], cl3);
    return z;
}

// Algorithm 17 in https://doi.org/10.1145/3121432
// Computes an approximation z ~= (xh + xl) / (yh + yl) such that the relative error is
// in O(u^2) (precise error bound not stated in source). Returns only the more
// significant component of z, so final relative error is on the order of u.
template <typename FPT>
FPT dwdivdw2_h(FPT xh, FPT xl, FPT yh, FPT yl)
{
    const FPT th = xh / yh;
    const auto r = dwtimesfp3(yh, yl, th);
    const FPT pih = xh - r[0];
    const FPT deltal = xl - r[1];
    const FPT delta = pih + deltal;
    const FPT tl = delta / yh;
    const auto z = fast_two_sum(th, tl);
    return z[0];
}

// Computes remainder r = r[0] + ... + r[remainder_size - 1] such that n / d = q + r / d
// for an approximation q of the quotient n / d.
template <typename FPT, std::size_t NSize, std::size_t DSize>
int compute_remainder(
        const std::array<FPT, NSize>& n,
        const std::array<FPT, DSize>& d,
        const FPT q_approx,
        std::array<FPT, NSize + DSize * 2>& remainder,
        int n_size = NSize,
        int d_size = DSize)
{
    std::array<FPT, DSize * 2> tmp;
    int tmp_size = scale_expansion_zeroelim(d, -q_approx, tmp, d_size);
    int remainder_size = fast_expansion_sum_zeroelim(n, tmp, remainder, n_size, tmp_size);
    return remainder_size;
}

// Computes an approximation q_approx of the quotient n / d for expansions
// n = n[0] + ... + n[n_size - 1], d = d[0] + ... + d[d_size - 1], such that the absolute
// error of q_approx is smaller than half of the unit in the last place, i.e. the 
// correctly rounded quotient under the assumption of round-to-nearest, tie-break to even
// rounding mode. Expects d_approx to be an approximation of d with relative error in
// O(u^2). Assumes that d is not zero.
template <bool CorrectRounding = true, typename FPT, std::size_t NSize, std::size_t DSize>
FPT divide_approx(const std::array<FPT, NSize>& n,
                  const std::array<FPT, DSize>& d,
                  const std::array<FPT, 2>& d_approx,
                  int n_size = NSize,
                  int d_size = DSize)
{
    const std::array<FPT, 2> n_approx = vec_sum_err_branch_2(n, n_size); //error in O(u²)
    FPT q_approx = dwdivdw2_h(n_approx[0], n_approx[1], d_approx[0], d_approx[1]);
    // This is an approximation with an error in O(u^2), which is then rounded again.
    // I think the error could be up to the unit in the last place (rather than half of
    // the unit in the last place), so will compute the remainder and based on its sign
    // a second candidate (nextafter(q_approx, +-INFINITY)) and test which candidate has
    // the smaller remainder.
    if (CorrectRounding)
    {
        std::array<FPT, NSize + DSize * 2> remainder;
        const int remainder_size = compute_remainder(n, d, q_approx, remainder, n_size, d_size);
        // We can skip all of this, if the result is exact <=> the remainder is zero.
        if(remainder[remainder_size - 1] != 0)
        {
            FPT q_approx2;
            // neither d[d_size - 1] nor remainder[remainder_size - 1] should be zero, so
            // the product in the second argument below should never be NaN.
            q_approx2 = std::nextafter(q_approx,
                                       INFINITY * d[d_size - 1] * remainder[remainder_size - 1]);
            std::array<FPT, NSize + DSize * 2> remainder2;
            int remainder2_size = compute_remainder(n, d, q_approx2, remainder2, n_size, d_size);
            // The two approximations should be the closest floating-point numbers to the true
            // result. The corresponding remainders should therefore have different signs.
            assert( (remainder[remainder_size - 1] > 0) - (remainder[remainder_size - 1] < 0) !=
                    (remainder2[remainder2_size - 1] > 0) - (remainder2[remainder2_size - 1] < 0));
            if (remainder[remainder_size - 1] < 0)
            {
                for(int i = remainder_size - 1; i >= 0; --i)
                    remainder[i] = -remainder[i]; // get abs(r)
            }
            if(remainder2[remainder2_size - 1] == 0)
            {
                return q_approx2; // I suspect this case may be impossible.
            }
            else
            {
                if(remainder2[remainder2_size - 1] > 0)
                {
                    for(int i = remainder2_size - 1; i >= 0; --i)
                        remainder2[i] = -remainder2[i]; // get abs(r2)
                }
                std::array<FPT, 2 * NSize + 4 * DSize> tmp2;
                int tmp2_size = fast_expansion_sum_zeroelim(remainder, remainder2, tmp2, remainder_size, remainder2_size);
                if(tmp2[tmp2_size - 1] > 0)
                {
                    return q_approx2; // second remainder is smaller, return second approx.
                }
                else if(tmp2[tmp2_size - 1] == 0)
                {
                    // The remainders are equal so the true quotient is exactly in the middle.
                    // We return the middle which will be correctly rounded towards one of the
                    // candidate by the FPU.
                    const FPT eps = (q_approx2 - q_approx) / 2;
                    return q_approx + eps;
                }
            }
        }
    }
    return q_approx;
}

template <typename FPT>
void intersection_common(FPT x1,
                             FPT y1,
                             FPT x2,
                             FPT y2,
                             FPT x3,
                             FPT y3,
                             FPT x4,
                             FPT y4,
                             std::array<FPT, 2>& x1y3,
                             std::array<FPT, 2>& x1y4,
                             std::array<FPT, 2>& x2y3,
                             std::array<FPT, 2>& x2y4,
                             std::array<FPT, 2>& x3y1,
                             std::array<FPT, 2>& x3y2,
                             std::array<FPT, 2>& x4y1,
                             std::array<FPT, 2>& x4y2)
{
    x1y3 = two_product(x1, y3);
    x1y4 = two_product(x1, y4);
    x2y3 = two_product(x2, y3);
    x2y4 = two_product(x2, y4);
    x3y1 = two_product(x3, y1);
    x3y2 = two_product(x3, y2);
    x4y1 = two_product(x4, y1);
    x4y2 = two_product(x4, y2);
}

template <typename FPT>
int intersection_denominator(FPT x1,
                             FPT y1,
                             FPT x2,
                             FPT y2,
                             FPT x3,
                             FPT y3,
                             FPT x4,
                             FPT y4,
                             const std::array<FPT, 2>& x1y3,
                             const std::array<FPT, 2>& x1y4,
                             const std::array<FPT, 2>& x2y3,
                             const std::array<FPT, 2>& x2y4,
                             const std::array<FPT, 2>& x3y1,
                             const std::array<FPT, 2>& x3y2,
                             const std::array<FPT, 2>& x4y1,
                             const std::array<FPT, 2>& x4y2,
                             std::array<FPT, 16>& d)
{
    const auto d1 = two_two_expansion_diff(x1y3, x1y4);
    const auto d2 = two_two_expansion_diff(x2y4, x2y3);
    const auto d3 = two_two_expansion_diff(x3y2, x3y1);
    const auto d4 = two_two_expansion_diff(x4y1, x4y2);
    int d_size;
    std::array<FPT, 8> d12;
    const auto d12_size = fast_expansion_sum_zeroelim(d1, d2, d12);
    std::array<FPT, 8> d34;
    const auto d34_size = fast_expansion_sum_zeroelim(d3, d4, d34);
    d_size = fast_expansion_sum_zeroelim(d12, d34, d, d12_size, d34_size);
    return d_size;
}

template <typename FPT>
int intersection_x_nominator(FPT x1,
                             FPT y1,
                             FPT x2,
                             FPT y2,
                             FPT x3,
                             FPT y3,
                             FPT x4,
                             FPT y4,
                             const std::array<FPT, 2>& x1y3,
                             const std::array<FPT, 2>& x1y4,
                             const std::array<FPT, 2>& x2y3,
                             const std::array<FPT, 2>& x2y4,
                             const std::array<FPT, 2>& x3y1,
                             const std::array<FPT, 2>& x3y2,
                             const std::array<FPT, 2>& x4y1,
                             const std::array<FPT, 2>& x4y2,
                             std::array<FPT, 32>& pxn)
{
    const auto r = [](const std::array<FPT, 2>& a){ return std::array<FPT, 2> { a[1], a[0] }; };
    int pxn_size;
    {
        std::array<FPT, 4> x1x3y2;
        const auto x1x3y2_size = scale_expansion_zeroelim(r(x3y2), x1, x1x3y2);
        std::array<FPT, 4> x1x4y2;
        const auto x1x4y2_size = scale_expansion_zeroelim(r(x4y2), -x1, x1x4y2);
        std::array<FPT, 8> pxn1;
        const auto pxn1_size =
            fast_expansion_sum_zeroelim(x1x3y2, x1x4y2, pxn1, x1x3y2_size, x1x4y2_size);
        std::array<FPT, 4> x2x3y1;
        const auto x2x3y1_size = scale_expansion_zeroelim(r(x3y1), -x2, x2x3y1);
        std::array<FPT, 4> x2x4y1;
        const auto x2x4y1_size = scale_expansion_zeroelim(r(x4y1), x2, x2x4y1);
        std::array<FPT, 8> pxn2;
        const auto pxn2_size =
            fast_expansion_sum_zeroelim(x2x3y1, x2x4y1, pxn2, x2x3y1_size, x2x4y1_size);
        std::array<FPT, 16> pxn12;
        const auto pxn12_size =
            fast_expansion_sum_zeroelim(pxn1, pxn2, pxn12, pxn1_size, pxn2_size);

        std::array<FPT, 4> x1x3y4;
        const auto x1x3y4_size = scale_expansion_zeroelim(r(x1y4), -x3, x1x3y4);
        std::array<FPT, 4> x1x4y3;
        const auto x1x4y3_size = scale_expansion_zeroelim(r(x1y3), x4, x1x4y3);
        std::array<FPT, 8> pxn3;
        const auto pxn3_size =
            fast_expansion_sum_zeroelim(x1x3y4, x1x4y3, pxn3, x1x3y4_size, x1x4y3_size);
        std::array<FPT, 4> x2x3y4;
        const auto x2x3y4_size = scale_expansion_zeroelim(r(x2y4), x3, x2x3y4);
        std::array<FPT, 4> x2x4y3;
        const auto x2x4y3_size = scale_expansion_zeroelim(r(x2y3), -x4, x2x4y3);
        std::array<FPT, 8> pxn4;
        const auto pxn4_size =
            fast_expansion_sum_zeroelim(x2x3y4, x2x4y3, pxn4, x2x3y4_size, x2x4y3_size);
        std::array<FPT, 16> pxn34;
        const auto pxn34_size =
            fast_expansion_sum_zeroelim(pxn3, pxn4, pxn34, pxn3_size, pxn4_size);
        pxn_size =
            fast_expansion_sum_zeroelim(pxn12, pxn34, pxn, pxn12_size, pxn34_size);
    }
    return pxn_size;
}

template <typename FPT>
int intersection_y_nominator(FPT x1,
                             FPT y1,
                             FPT x2,
                             FPT y2,
                             FPT x3,
                             FPT y3,
                             FPT x4,
                             FPT y4,
                             const std::array<FPT, 2>& x1y3,
                             const std::array<FPT, 2>& x1y4,
                             const std::array<FPT, 2>& x2y3,
                             const std::array<FPT, 2>& x2y4,
                             const std::array<FPT, 2>& x3y1,
                             const std::array<FPT, 2>& x3y2,
                             const std::array<FPT, 2>& x4y1,
                             const std::array<FPT, 2>& x4y2,
                             std::array<FPT, 32>& pyn)
{
    const auto r = [](const std::array<FPT, 2>& a){ return std::array<FPT, 2> { a[1], a[0] }; };
    int pyn_size;
    {
        std::array<FPT, 4> x1y2y3;
        const auto x1y2y3_size = scale_expansion_zeroelim(r(x1y3), y2, x1y2y3);
        std::array<FPT, 4> x1y2y4;
        const auto x1y2y4_size = scale_expansion_zeroelim(r(x1y4), -y2, x1y2y4);
        std::array<FPT, 8> pyn1;
        const auto pyn1_size =
            fast_expansion_sum_zeroelim(x1y2y3, x1y2y4, pyn1, x1y2y3_size, x1y2y4_size);
        std::array<FPT, 4> x2y1y3;
        const auto x2y1y3_size = scale_expansion_zeroelim(r(x2y3), -y1, x2y1y3);
        std::array<FPT, 4> x2y1y4;
        const auto x2y1y4_size = scale_expansion_zeroelim(r(x2y4), y1, x2y1y4);
        std::array<FPT, 8> pyn2;
        const auto pyn2_size =
            fast_expansion_sum_zeroelim(x2y1y3, x2y1y4, pyn2, x2y1y3_size, x2y1y4_size);
        std::array<FPT, 16> pyn12;
        const auto pyn12_size =
            fast_expansion_sum_zeroelim(pyn1, pyn2, pyn12, pyn1_size, pyn2_size);

        std::array<FPT, 4> x3y1y4;
        const auto x3y1y4_size = scale_expansion_zeroelim(r(x3y1), -y4, x3y1y4);
        std::array<FPT, 4> x4y1y3;
        const auto x4y1y3_size = scale_expansion_zeroelim(r(x4y1), y3, x4y1y3);
        std::array<FPT, 8> pyn3;
        const auto pyn3_size =
            fast_expansion_sum_zeroelim(x3y1y4, x4y1y3, pyn3, x3y1y4_size, x4y1y3_size);
        std::array<FPT, 4> x3y2y4;
        const auto x3y2y4_size = scale_expansion_zeroelim(r(x3y2), y4, x3y2y4);
        std::array<FPT, 4> x4y2y3;
        const auto x4y2y3_size = scale_expansion_zeroelim(r(x4y2), -y3, x4y2y3);
        std::array<FPT, 8> pyn4;
        const auto pyn4_size =
            fast_expansion_sum_zeroelim(x3y2y4, x4y2y3, pyn4, x3y2y4_size, x4y2y3_size);
        std::array<FPT, 16> pyn34;
        const auto pyn34_size =
            fast_expansion_sum_zeroelim(pyn3, pyn4, pyn34, pyn3_size, pyn4_size);
        pyn_size = fast_expansion_sum_zeroelim(pyn12, pyn34, pyn, pyn12_size, pyn34_size);
    }
    return pyn_size;
}

// Returns true iff lines (x1, y1), (x2, y2) and (x3, y3), (x4, y4) are guaranteed to be
// neither parallel nor correctly rounded versions of parallel lines, returns false otherwise.
// No guarantees in case of overflow or underflow.
template <typename FPT>
bool non_parallel_rounded(FPT x1,
                          FPT y1,
                          FPT x2,
                          FPT y2,
                          FPT x3,
                          FPT y3,
                          FPT x4,
                          FPT y4)
{
    const FPT d = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
    constexpr FPT eps = std::numeric_limits<FPT>::epsilon() / 2;
    const FPT error =   (5 * eps + 32 * eps * eps)
                      * (   (std::abs(x1) + std::abs(x2)) * (std::abs(y3) + std::abs(y4))
                          + (std::abs(y1) + std::abs(y2)) * (std::abs(x3) + std::abs(x4)));
    return std::abs(d) > std::abs(error);
}

// Returns true iff lines (x1, y1), (x2, y2) and (x3, y3), (x4, y4) are parallel, false iff the lines are not parallel.
// No guarantees in case of overflow or underflow.
template <typename FPT>
bool non_parallel(FPT x1, FPT y1, FPT x2, FPT y2, FPT x3, FPT y3, FPT x4, FPT y4)
{
    const FPT d = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
    constexpr FPT eps = std::numeric_limits<FPT>::epsilon() / 2;
    const FPT error =   (3 * eps + 16 * eps * eps)
                      * (   std::abs((x1 - x2) * (y3 - y4))
                          + std::abs((y1 - y2) * (x3 - x4)) );
    if ( std::abs(d) > error )
        return true;
    if ( d == 0 && error == 0 )
        return false;
    std::array<FPT, 2> x1y3, x1y4, x2y3, x2y4, x3y1, x3y2, x4y1, x4y2;
    intersection_common(x1, y1, x2, y2, x3, y3, x4, y4, x1y3, x1y4, x2y3, x2y4, x3y1, x3y2, x4y1,
                        x4y2);
    std::array<FPT, 16> exact_d;
    auto ds = intersection_denominator(x1, y1, x2, y2, x3, y3, x4, y4, x1y3, x1y4, x2y3, x2y4,
                                       x3y1, x3y2, x4y1, x4y2, exact_d);
    return exact_d[ds - 1] != 0;
}

// Compute the intersection (px, py) of lines (x1, y1), (x2, y2) and (x3, y3), (x4, y4) with
// correctly rounded coordinates, i.e. up to the available precision.
// Assumes that a unique intersection exists (i.e. d != 0).
template <typename FPT, bool CorrectRounding = true>
std::array<FPT, 2> intersection_robust(FPT x1,
                                       FPT y1,
                                       FPT x2,
                                       FPT y2,
                                       FPT x3,
                                       FPT y3,
                                       FPT x4,
                                       FPT y4)
{
    assert(non_parallel(x1, y1, x2, y2, x3, y3, x4, y4));
    std::array<FPT, 2> x1y3, x1y4, x2y3, x2y4, x3y1, x3y2, x4y1, x4y2;
    intersection_common(x1, y1, x2, y2, x3, y3, x4, y4, x1y3, x1y4, x2y3, x2y4, x3y1, x3y2, x4y1,
                        x4y2);
    std::array<FPT, 16> d;
    int d_size = intersection_denominator(x1, y1, x2, y2, x3, y3, x4, y4, x1y3, x1y4, x2y3, x2y4,
                                          x3y1, x3y2, x4y1, x4y2, d);
    std::array<FPT, 2> d_approx = vec_sum_err_branch_2(d, d_size); //error in O(u²)

    FPT px_approx;
    {
        std::array<FPT, 32> pxn;
        int pxn_size = intersection_x_nominator(x1, y1, x2, y2, x3, y3, x4, y4, x1y3, x1y4, x2y3,
                                                x2y4, x3y1, x3y2, x4y1, x4y2, pxn);
        px_approx = divide_approx<CorrectRounding>(pxn, d, d_approx, pxn_size, d_size);
    }

    FPT py_approx;
    {
        std::array<FPT, 32> pyn;
        int pyn_size = intersection_y_nominator(x1, y1, x2, y2, x3, y3, x4, y4, x1y3, x1y4, x2y3,
                                                x2y4, x3y1, x3y2, x4y1, x4y2, pyn);
        py_approx = divide_approx<CorrectRounding>(pyn, d, d_approx, pyn_size, d_size);
    }
    return std::array<FPT, 2> {px_approx, py_approx};
}

// Computes line intersections with no guarantees wrt to accuracy.
template <typename FPT>
std::array<FPT, 2> intersection_nonrobust(FPT x1,
                                          FPT y1,
                                          FPT x2,
                                          FPT y2,
                                          FPT x3,
                                          FPT y3,
                                          FPT x4,
                                          FPT y4)
{
    const FPT d = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
    const FPT xn = (x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4);
    const FPT yn = (x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4);
    return std::array<FPT, 2> { xn / d, yn / d };
}

// Computes line intersections with an relative error of less than 200 * u. This is fast in most
// cases (assuming mostly non-degenerate input) but slow in some cases.
// Should be consistent with side_rounded_input<FPT, 404, 0>. Assumes that a unique intersection
// exists (i.e. d != 0).
template <typename FPT>
std::array<FPT, 2> intersection_filtered(FPT x1,
                                         FPT y1,
                                         FPT x2,
                                         FPT y2,
                                         FPT x3,
                                         FPT y3,
                                         FPT x4,
                                         FPT y4)
{
    assert(non_parallel(x1, y1, x2, y2, x3, y3, x4, y4));
    FPT d = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
    FPT xn = (x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4);
    FPT yn = (x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4);
    constexpr FPT u = std::numeric_limits<FPT>::epsilon() / 2;
    constexpr FPT d_coeff = 4 * u + 12 * u * u;
    constexpr FPT n_coeff = 5 * u + 24 * u * u;
    const FPT d_mag = std::abs((x1 - x2) * (y3 - y4)) + std::abs((y1 - y2) * (x3 - x4));
    const FPT xn_mag = (std::abs(x1 * y2) + std::abs(y1 * x2)) * std::abs(x3 - x4) + std::abs(x1 - x2) * (std::abs(x3 * y4) + std::abs(y3 * x4));
    const FPT yn_mag = (std::abs(x1 * y2) + std::abs(y1 * x2)) * std::abs(y3 - y4) + std::abs(y1 - y2) * (std::abs(x3 * y4) + std::abs(y3 * x4));
    const FPT d_err = d_coeff * d_mag;
    const FPT xn_err = n_coeff * xn_mag;
    const FPT yn_err = n_coeff * yn_mag;
    constexpr FPT bound = 99 * u;
    if (    std::abs(d_err / d) < bound
         && ( xn_err == 0. || std::abs(xn_err / xn) < bound)
         && ( yn_err == 0. || std::abs(yn_err / yn) < bound) )
        return std::array<FPT, 2> { xn / d, yn / d };
    else
    {
        std::array<FPT, 2> x1y3, x1y4, x2y3, x2y4, x3y1, x3y2, x4y1, x4y2;
        intersection_common(x1, y1, x2, y2, x3, y3, x4, y4, x1y3, x1y4, x2y3, x2y4, x3y1, x3y2,
                            x4y1, x4y2);
        if ( std::abs(d_err / d) >= bound )
        {
            std::array<FPT, 16> d_exact;
            int d_size = intersection_denominator(x1, y1, x2, y2, x3, y3, x4, y4, x1y3, x1y4, x2y3,
                                                  x2y4, x3y1, x3y2, x4y1, x4y2, d_exact);
            d = std::accumulate( d_exact.begin(), d_exact.begin() + d_size, FPT(0) );
        }
        if ( std::abs(xn_err / xn) >= bound && xn_err != 0. )
        {
            std::array<FPT, 32> xn_exact;
            int xn_size = intersection_x_nominator(x1, y1, x2, y2, x3, y3, x4, y4, x1y3, x1y4, x2y3,
                                                   x2y4, x3y1, x3y2, x4y1, x4y2, xn_exact);
            xn = std::accumulate( xn_exact.begin(), xn_exact.begin() + xn_size, FPT(0) );
        }
        if ( std::abs(yn_err / yn) >= bound && yn_err != 0. )
        {
            std::array<FPT, 32> yn_exact;
            int yn_size = intersection_y_nominator(x1, y1, x2, y2, x3, y3, x4, y4, x1y3, x1y4, x2y3,
                                                   x2y4, x3y1, x3y2, x4y1, x4y2, yn_exact);
            yn = std::accumulate( yn_exact.begin(), yn_exact.begin() + yn_size, FPT(0) );
        }
        return {{ xn / d, yn / d }};
    }
}

}}

}}

