

#include "lisa_math.hpp"

#include <stdexcept>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>


namespace lisa {


void forward_fft(std::vector<double>& v) {

  assert(v.size() != 0); // Else the following check won't work.
  assert(isPowerOfTwo(v.size()));

  int status = gsl_fft_complex_radix2_forward(v.data(), 1, v.size() / 2);
  if (status) { throw std::runtime_error("GSL error"); }

}


/// https://stackoverflow.com/questions/108318/whats-the-simplest-way-to-test-whether-a-number-is-a-power-of-2-in-c
bool isPowerOfTwo(int i) {
  return (i & (i - 1)) == 0;
}


}

