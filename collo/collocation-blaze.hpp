#ifndef COLLOCATION_HPP
#define COLLOCATION_HPP

#include <blaze/math/StaticVector.h>
#include <blaze/math/StaticMatrix.h>
#include <tuple>

namespace collocation
{

    template<typename num_t,
	     std::size_t system_order,
	     std::size_t method_stage>
    class Collocation
    {
	protected:
	    
	    using state_vector_t = blaze::StaticVector<num_t, system_order>;
	    using vector_t = blaze::StaticVector<num_t, method_stage>;
	    using matrix_t = blaze::StaticMatrix<num_t, method_stage, method_stage>;

	private:

	    using sva_t = blaze::StaticMatrix<num_t, method_stage, system_order>; // type of state_vector's array

	    static num_t distance(const sva_t& first, const sva_t& second)
	    {
		return blaze::norm(first - second);
	    }

	    static state_vector_t shift(const sva_t& alphas)
	    {
		return blaze::trans( blaze::sum<blaze::columnwise>( alphas ) );
	    }

	    static sva_t iteration_step(const sva_t& alphas)
	    {
		sva_t alphas_prev = alphas;
		return alphas_prev;
	    }

	protected:

	    template<typename rhs_t>
	    static auto do_step_impl(const state_vector_t& y0, const rhs_t& rhs, const sva_t& alphas_init)
	    {
		auto alphas = make_initial_alphas(y0, alphas_init);

		std::size_t iternum = 0;
		num_t dist;

		do
		{
		    auto alphas_prev = iteration_step(alphas);
		    dist = distance(alphas, alphas_prev);
		    ++iternum;
		}
		while(dist + num_t{1.0} != num_t{1.0});

		auto y_shift = shift(alphas);

		return std::make_tuple(
		           y_shift,
			   y0 + y_shift,
			   iternum,
			   alphas
		       );
	    }

    };

}

#endif //COLLOCATION_HPP
