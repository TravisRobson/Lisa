

#include <cmath>
#include <vector>

#include "lisa_constants.hpp"
#include "lisa_detector.hpp"
#include "lisa_noise.hpp"
#include "lisa_types.hpp"

#include "gtest/gtest.h"


using namespace lisa;


constexpr int samples_1_year {2'097'152};
constexpr double obs_period { samples_1_year * lisa::constants::delta_t };


TEST(TDIPowerFitTests, Dummy)
{

  // double min_freq {  0.1e-3 }; // 0.1 mHz
  // double max_freq { 20.0e-3 }; // 20 mHz

  double min_freq { 1.0e-3 }; // 1.0 mHz
  double max_freq { 2.0e-3 }; // 2.0 mHz

  // Find the freq domain index corresponding to the above frequency bounds.
  size_t min_index { static_cast<size_t>(std::floor(min_freq * obs_period)) };
  size_t max_index { static_cast<size_t>(std::ceil(max_freq  * obs_period)) };

  assert(min_index < man_index);
  size_t num_samples { max_index - min_index + 1 };

  size_t samples_per_segment { 100 + 1 }; // We demand it's odd because we're performing a median
                                          // using samples_per_segment values.

  int num_segments { static_cast<int>(static_cast<double>(num_samples) / static_cast<double>(samples_per_segment)) };


  FitState state {
    std::vector<double>(num_samples),
    std::vector<double>(num_samples),
    std::vector<double>(samples_per_segment),
    std::vector<double>(num_segments),
    std::vector<double>(num_segments),
    std::vector<double>(num_segments),
    std::vector<double>(num_segments),
    std::vector<double>(num_segments),
    min_index,
    max_index
  };

  /// \todo Put in instrument noise itself.

  size_t half_segment { (samples_per_segment - 1) / 2 };

  RealFreqVector tdi_power(min_index - half_segment, max_index + half_segment);

  for (size_t i = tdi_power.min(); i <= tdi_power.max(); ++i) {
    double freq { static_cast<double>(i) / obs_period };
    tdi_power[ i ] = analytic_instrument_noise_legacy(freq).XYZ;
  }

  fit_tdi_power(tdi_power, state, obs_period);

  for (size_t i = 0; i < state.means.size(); ++i) {
    //std::cout << std::scientific << state.means[i] << "\n";
  }

}
