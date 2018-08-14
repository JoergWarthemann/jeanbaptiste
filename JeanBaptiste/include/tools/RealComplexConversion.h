#pragma once

//#include "../../../gsl/include/gsl/span"
#include "../../../gsl/span"

namespace jeanbaptiste::tools
{
    template<typename SampleCnt,
             typename Complex>
	class Real2Complex
	{
		using ValueType = typename Complex::value_type;

	public:
		void operator()(gsl::span<const ValueType> realIn, gsl::span<Complex> complexOut) const
		{
			int i = 0;
			//if ((realIn != nullptr) && (complexOut != nullptr))
			//	for (unsigned i = 0; i < SampleCnt; ++i, ++realIn)
			//		complexOut[i].real(*realIn);
		}
	};
}