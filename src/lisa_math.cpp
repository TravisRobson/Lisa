

#include "lisa_errors.hpp"
#include "lisa_math.hpp"

#include <assert.h>
#include <stdexcept>
#include <sstream>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>


static const gsl_rng_type* gslRngType { gsl_rng_default };
static gsl_rng* pGslRng { gsl_rng_alloc(gslRngType) };


namespace lisa {


///////////////////////////////////////////////////////////////////////////////
/// \brief Wrapper function for thirdparty libraries.
///
/// Will allow us to easily swap out for other libraries in the future.
///////////////////////////////////////////////////////////////////////////////
void forward_fft(std::vector<double>& v) {

  assert(v.size() != 0); // Else the following check won't work.
  assert(isPowerOfTwo(v.size()));

  constexpr int stride { 1 };

  int status = gsl_fft_complex_radix2_forward(v.data(), stride, v.size() / 2);
  if (status) {
    std::stringstream ss;
    ss << "gsl_fft_complex_radix2_forward failure, exit code " << status << ".";
    THROW_ERROR(ss.str().c_str());
  }

}


/// https://stackoverflow.com/questions/108318/whats-the-simplest-way-to-test-whether-a-number-is-a-power-of-2-in-c
bool isPowerOfTwo(int i) {
  return (i & (i - 1)) == 0;
}


double gaussian_variate(double sigma)
{
  return gsl_ran_gaussian(pGslRng, sigma);
}


}


/// \todo I need to have some object attached to this probably...
//gsl_rng_free(pGslRng);





