#ifndef NUMM_LEGENDRE_HPP
#define NUMM_LEGENDRE_HPP

#include <array>
#include <utility>
#include <numbers>
#include "roots.hpp"

namespace numm{

    template<std::size_t n, std::size_t first_order = 0, number_type num_t>
	requires (first_order <= n)
    constexpr
    auto
    legendre(num_t x)
    {
	return ortho_poly<n, first_order, num_t>(x, {num_t{1.0}, x},
		[](num_t x, auto& res)
		{
		    for(std::size_t i = 1; i < n; ++i)
			res[i+1] = ( (2 * i + 1) * x * res[i] - i * res[i-1] ) / (i + 1);
		});
    }

    template<std::size_t n, std::size_t first_order = 0, std::size_t m, number_type num_t>
	requires (first_order <= n)
    constexpr
    auto
    legendre(const std::array<num_t, m>& xs)
    {
	return ortho_poly<n, first_order, m, num_t>(xs, legendre<n, first_order, num_t>);
    }

    template<std::size_t n, std::size_t first_order = 0, number_type num_t>
	requires (first_order <= n)
    constexpr
    auto
    dlegendre(num_t x)
    {
	return ortho_poly<n, first_order, num_t>(x, {num_t{0.0}, num_t{1.0}},
		[](num_t x, auto& res)
		{
		    auto l = legendre<n - 1>(x);
		    for(std::size_t i = 1; i < n; ++i)
			res[i+1] = (i + 1) * l[i] + x * res[i];
		});
    }

    template<std::size_t n, std::size_t first_order = 0, std::size_t m, number_type num_t>
	requires (first_order <= n)
    constexpr
    std::array<std::array<num_t, n - first_order + 1>, m>
    dlegendre(const std::array<num_t, m>& xs)
    {
	return ortho_poly<n, first_order, m, num_t>(xs, dlegendre<n, first_order, num_t>);
    }

    template<std::size_t n, std::size_t first_order = 0, number_type num_t>
	requires (first_order <= n)
    constexpr
    auto
    d2legendre(num_t x)
    {
	return ortho_poly<n, first_order, num_t>(x, {num_t{0.0}, num_t{0.0}},
		[](num_t x, auto& res)
		{
		    auto dl = dlegendre<n - 1>(x);
		    for(std::size_t i = 1; i < n; ++i)
			res[i+1] = (i + 2) * dl[i] + x * res[i];
		});
    }

    template<std::size_t n, std::size_t first_order = 0, std::size_t m, number_type num_t>
	requires (first_order <= n)
    constexpr
    std::array<std::array<num_t, n - first_order + 1>, m>
    d2legendre(const std::array<num_t, m>& xs)
    {
	return ortho_poly<n, first_order, m, num_t>(xs, d2legendre<n, first_order, num_t>);
    }

    template<std::size_t n, std::size_t first_order = 0, number_type num_t>
	requires (first_order <= n)
    constexpr
    auto
    legendre_sh(num_t y) // x = 2 * y - 1
    {
	return legendre<n, first_order, num_t>(num_t{2.0} * y - num_t{1.0});
    }

    template<std::size_t n, std::size_t first_order = 0, std::size_t m, number_type num_t>
	requires (first_order <= n)
    constexpr
    std::array<std::array<num_t, n - first_order + 1>, m>
    legendre_sh(const std::array<num_t, m>& ys)
    {
	return ortho_poly<n, first_order, m, num_t>(ys, legendre_sh<n, first_order, num_t>);
    }

    template<std::size_t n, std::size_t first_order = 0, number_type num_t>
	requires (first_order <= n)
    constexpr
    auto
    dlegendre_sh(num_t y) // x = 2 * y - 1
    {
	return ortho_poly<n, first_order, num_t>(y, {num_t{0.0}, num_t{2.0}},
		[](num_t y, auto& res)
		{
		    auto l = legendre_sh<n - 1>(y);
		    for(std::size_t i = 1; i < n; ++i)
			res[i+1] = 2 * (i + 1) * l[i] + (2 * y - 1) * res[i];
		});
    }

    template<std::size_t n, std::size_t first_order = 0, std::size_t m, number_type num_t>
	requires (first_order <= n)
    constexpr
    std::array<std::array<num_t, n - first_order + 1>, m>
    dlegendre_sh(const std::array<num_t, m>& ys)
    {
	return ortho_poly<n, first_order, m, num_t>(ys, dlegendre_sh<n, first_order, num_t>);
    }

    template<std::size_t n, std::size_t first_order = 0, number_type num_t>
	requires (first_order <= n)
    constexpr
    auto
    d2legendre_sh(num_t y) // x = 2 * y - 1
    {
	return ortho_poly<n, first_order, num_t>(y, {num_t{0.0}, num_t{0.0}},
		[](num_t y, auto& res)
		{
		    auto dl = dlegendre_sh<n - 1>(y);
		    for(std::size_t i = 1; i < n; ++i)
			res[i+1] = 2 * (i + 2) * dl[i] + (2 * y - 1) * res[i];
		});
    }

    template<std::size_t n, std::size_t first_order = 0, std::size_t m, number_type num_t>
	requires (first_order <= n)
    constexpr
    std::array<std::array<num_t, n - first_order + 1>, m>
    d2legendre_sh(const std::array<num_t, m>& ys)
    {
	return ortho_poly<n, first_order, m, num_t>(ys, d2legendre_sh<n, first_order, num_t>);
    }

    template<std::size_t n, number_type num_t = double>
    consteval
    std::array<num_t, n>
    roots_legendre()
    {
	if constexpr (n == 1) return { num_t{0.0} };
	else
	{
	    num_t c{1.0};
	    c -= num_t{1.0} / 8 / n / n;
	    c += num_t{1.0} / 8 / n / n / n;
	    auto pi = std::numbers::pi_v<num_t>;
	    std::array<num_t, n> roots;
	    for(std::size_t k = 0; k < n; ++k)
		roots[k] = newton( c * cos( pi * (4 * (n - k) - 1) / (4 * n + 2) ),
					[](num_t x){return legendre<n>(x)[n];},
					[](num_t x){return dlegendre<n>(x)[n];} );
	    return roots;
	}
    }

    template<std::size_t n, number_type num_t = double>
    consteval
    std::array<num_t, n-1>
    roots_dlegendre()
    {
	if constexpr (n == 2) return { num_t{0.0} };
	else
	{
	    auto pi = std::numbers::pi_v<num_t>;
	    std::array<num_t, n - 1> roots;
	    for(std::size_t k = 0; k < n - 1; ++k)
		roots[k] = newton( cos( pi * (n - k - 1) / n ),
					[](num_t x){return dlegendre<n>(x)[n];},
					[](num_t x){return d2legendre<n>(x)[n];} );
	    return roots;
	}
    }

    template<std::size_t n, number_type num_t = double>
    consteval
    std::array<num_t, n+1>
    roots_ilegendre()
    {
	if constexpr (n == 1) return { num_t{-1.0}, num_t{1.0} };
	else
	{
	    std::array<num_t, n + 1> roots;
	    auto droots = roots_dlegendre<n, num_t>();
	    roots[0] = num_t{-1.0};
	    for(std::size_t k = 0; k < n - 1; ++k) roots[k + 1] = droots[k];
	    roots[n] = num_t{1.0};
	    return roots;
	}
    }

    template<std::size_t n, number_type num_t = double>
    consteval
    std::array<num_t, n>
    roots_legendre_sh()
    {
	if constexpr (n == 1) return { num_t{0.5} };
	else
	{
	    std::array<num_t, n> roots = roots_legendre<n, num_t>();
	    for(std::size_t k = 0; k < n; ++k)
		roots[k] = newton( ( roots[k] + 1 ) / 2,
					[](num_t x){return legendre_sh<n>(x)[n];},
					[](num_t x){return dlegendre_sh<n>(x)[n];} );
	    return roots;
	}
    }

    template<std::size_t n, number_type num_t = double>
    consteval
    std::array<num_t, n - 1>
    roots_dlegendre_sh()
    {
	if constexpr (n == 2) return { num_t{0.5} };
	else
	{
	    std::array<num_t, n-1> roots = roots_dlegendre<n, num_t>();
	    for(std::size_t k=0; k<n-1; ++k)
		roots[k] = newton( ( roots[k] + 1 ) / 2,
					[](num_t x){return dlegendre_sh<n>(x)[n];},
					[](num_t x){return d2legendre_sh<n>(x)[n];} );
	    return roots;
	}
    }

    template<std::size_t n, number_type num_t = double>
    consteval
    std::array<num_t, n + 1>
    roots_ilegendre_sh()
    {
	if constexpr (n == 1) return { num_t{0.0}, num_t{1.0} };
	else
	{
	    std::array<num_t, n + 1> roots;
	    auto droots = roots_dlegendre_sh<n, num_t>();
	    roots[0] = num_t{0.0};
	    for(std::size_t k = 0; k < n - 1; ++k) roots[k + 1] = droots[k];
	    roots[n] = num_t{1.0};
	    return roots;
	}
    }

}

#endif //NUMM_LEGENDRE_HPP
