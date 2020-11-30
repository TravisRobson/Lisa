

#include <vector>

#include "lisa_constants.hpp"
#include "lisa_detector.hpp"
#include "lisa_types.hpp"
#include "lisa_utils.hpp"

#include "gtest/gtest.h"


using namespace lisa;


#ifdef BUILD_LEGACY_LISA
extern "C" void FAST_LISA(double *params, long N, long M, double *XLS, double *ALS, double *ELS);
#endif


///////////////////////////////////////////////////////////////////////////////
/// \brief
///////////////////////////////////////////////////////////////////////////////
#ifdef BUILD_LEGACY_LISA
TEST(GBSignalTests, LegacyRegression)
{

  GB_params params {
    9.252046e-05, // freq
    6.704986e-22, // freq_dot
    0.5 * constants::pi - (-7.769250e-01), // theta
    4.129078e+00, // phi
    7.813862e-24, // amp
    1.720600e-01, // iota
    5.658070e-01, // psi
    3.390864e+00, // phase
    0.0
  };

  params.freq_double_dot = 11.0 / 3.0 * params.freq_dot * params.freq_dot / params.freq;

  int num_samples { 512 };

  std::vector<double> X(2 * num_samples);
  std::vector<double> A(2 * num_samples);
  std::vector<double> E(2 * num_samples);

  double obs_time { 2'097'152 * lisa::constants::delta_t };

  calc_gb_signal(params, X, A, E, obs_time);

  std::vector<double> X_legacy(2 * num_samples + 1);
  std::vector<double> A_legacy(2 * num_samples + 1);
  std::vector<double> E_legacy(2 * num_samples + 1);

  double params_legacy[9] = { 
    params.freq,
    params.theta,
    params.phi,
    params.amp,
    params.iota,
    params.psi,
    params.phase,
    params.freq_dot,
    params.freq_double_dot
  };

  FAST_LISA(
    params_legacy, 
    num_samples, 
    num_samples, 
    X_legacy.data(), 
    A_legacy.data(), 
    E_legacy.data());


  // First two expect, third assert to kill test. Don't want
  // too many printouts.
  for (int i = 0; i < X.size(); ++i) {
    EXPECT_FLOAT_EQ(X[i], X_legacy[i + 1]) << "[i] = [" << i << "]";
    EXPECT_FLOAT_EQ(A[i], A_legacy[i + 1]) << "[i] = [" << i << "]";
    ASSERT_FLOAT_EQ(E[i], E_legacy[i + 1]) << "[i] = [" << i << "]";
  }

}
#endif // #ifdef BUILD_LEGACY_LISA


#ifdef ENABLE_TIMING_TESTS
///////////////////////////////////////////////////////////////////////////////
/// \brief 
///
/// For this specific source the legacy median runtime (compiled -O3 NDEBUG) 
/// on my 2 GHz Macbook Pro was 18.229 ms.
///////////////////////////////////////////////////////////////////////////////
TEST(GBSignalTests, Time)
{

  GB_params params {
    9.252046e-05, // freq
    6.704986e-22, // freq dot
    0.5 * constants::pi - (-7.769250e-01), // theta
    4.129078e+00, // phi
    7.813862e-24, // amp
    1.720600e-01, // iota
    5.658070e-01, // psi
    3.390864e+00, // phase
    0.0           // freq double dot
  };

  params.freq_double_dot = 11.0 / 3.0 * params.freq_dot * params.freq_dot / params.freq;

  int num_samples { 2 * 8192 };

  std::vector<double> X(num_samples);
  std::vector<double> A(num_samples);
  std::vector<double> E(num_samples);

  double obs_time { 2'097'152 * lisa::constants::delta_t };

  int num_iterations { 100 };

  Timer timer;

  for (int i = 0; i < num_iterations; ++i) {
    timer.start();
    calc_gb_signal(params, X, A, E, obs_time);
    timer.stop();
  }

  double to_milliseconds { 1.0e-6 };

  std::cout << "Median GB time (" << num_samples / 2 << " time samples) " 
            << timer.median_time() * to_milliseconds << " ms \n";

  ASSERT_LE(timer.median_time() * to_milliseconds, 10.0);


#ifdef BUILD_LEGACY_LISA
  double params_legacy[9] = { 
    params.freq,
    params.theta,
    params.phi,
    params.amp,
    params.iota,
    params.psi,
    params.phase,
    params.freq_dot,
    params.freq_double_dot
  };


  Timer timer_legacy;

  for (int i = 0; i < num_iterations; ++i) {
    timer_legacy.start();

    FAST_LISA(
      params_legacy, 
      num_samples, 
      num_samples, 
      X.data(), 
      A.data(), 
      E.data());

    timer_legacy.stop();
  }

  std::cout << "Median GB LEGACY time (" << num_samples / 2 << " time samples) " 
            << timer_legacy.median_time() * to_milliseconds << " ms \n";
#endif // BUILD_LEGACY_LISA


}
#endif // ENABLE_TIMING_TESTS
