#pragma once

#include <vector>
#include <complex>

namespace jeanbaptiste::tools
{
    template<typename Complex>
	class Real2Complex
	{
		using ValueType = typename Complex::value_type;

	public:
		std::vector<Complex> operator()(const std::vector<ValueType>& realIn) const
		{
			std::vector<Complex> result(realIn.size());

			for (std::size_t i = 0, iEnd = realIn.size(); i < iEnd; ++i)
				result[i].real(realIn[i]);

			return result;
		}
	};

	template<typename Complex>
	class Complex2Real
	{
		using ValueType = typename Complex::value_type;

	public:
		std::vector<ValueType> operator()(const std::vector<Complex>& complexIn) const
		{
			std::vector<ValueType> result(complexIn.size());

			for (std::size_t i = 0, iEnd = complexIn.size(); i < iEnd; ++i)
				result[i] = complexIn[i].real();

			return result;
		}
	};
}