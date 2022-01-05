#ifndef NUMM_BASE_HPP
#define NUMM_BASE_HPP

#include <concepts>

namespace numm{

    template < typename T >
    concept number_type = std::floating_point<T>;

    /*template<number_type num_t = double>
    consteval num_t pi()
    {
	num_t res{3.0}, corr{1.0/8.0};
	for(std::size_t n=2; res+corr != res; ++n )
	{
	    res += corr;
	    corr /= 8 * n * (2*n+1);
	    corr *= (2*n-1) * (2*n-1);
	}
	return res;
    }*/

    template<number_type num_t>
    constexpr num_t cos(const num_t& x)
    {
	num_t res{1.0}, x2{x * x}, corr{-x2 / 2.0};
	std::size_t n = 2;
	for(std::size_t n = 2; res + corr != res; ++n)
	{
	    res += corr;
	    corr /= (2*n) * (2*n-1);
	    corr *= -x2;
	}
	return res;
    }

    template<std::size_t k, number_type num_t, std::size_t n>
	requires (k < n)
    constexpr
    std::array<num_t, n - k>
    cutoff_first(const std::array<num_t, n>& arr)
    {
	std::array<num_t, n - k> res;
	for(std::size_t i = 0; i < n - k; ++i)
	    res[i] = arr[i + k];
	return res;
    }

    template<std::size_t n, std::size_t first_order = 0, number_type num_t>
	requires (first_order <= n)
    constexpr
    std::array<num_t, n - first_order + 1>
    ortho_poly(num_t x, std::array<num_t, 2>&& base, auto&& fill)
    {
	if constexpr (n == 0) return { base[0] };
	else if constexpr (n == 1 and first_order == 0) return { base[0], base[1] };
	else if constexpr (n == 1 and first_order == 1) return { base[1] };
	else
	{
	    std::array<num_t, n + 1> res{base[0], base[1]};
	    fill(x, res);
	    if constexpr (first_order > 0)
		return cutoff_first<first_order>(res);
	    else
		return res;
	}
    }

    template<std::size_t n, std::size_t first_order = 0, std::size_t m, number_type num_t>
	requires (first_order <= n)
    constexpr
    std::array<std::array<num_t, n - first_order + 1>, m>
    ortho_poly(const std::array<num_t, m>& xs, auto&& op)
    {
	std::array<std::array<num_t, n - first_order + 1>, m> res;
	for(std::size_t i = 0; i < m; ++i)
	    res[i] = op(xs[i]);
	return res;
    }

}

#endif //NUMM_BASE_HPP
