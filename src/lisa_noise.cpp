

#include "lisa_noise.hpp"
#include "lisa_utils.hpp"

#include <cmath>
#include <stdexcept>

#include "lisa_constants.hpp"
#include "lisa_detector.hpp"
#include "lisa_utils.hpp"


#include <iostream> /// \todo remove.
namespace lisa {


void fit_tdi_power(const RealFreqVector& power, FitState& state, double obs_period)
{

  //
  // First we must calculate the mean of log-periodogram, it's median, and
  // the standard deviation of the mean itself.
  //
  size_t num_segments { state.log_freqs.size() };
  size_t samples_per_segment { state.segment.size() };
  size_t half_samples { (samples_per_segment - 1) / 2 };
  size_t max_freq_index { state.min_freq_index + state.segment.size() * num_segments };

  /// \todo Do something more sophisitcated if this error is hit.
  if (max_freq_index > state.max_freq_index)
    throw std::runtime_error("We can't deal with this!!!");

  for (size_t seg_idx = 0; seg_idx < num_segments; ++seg_idx) {

    size_t freq_index { state.min_freq_index + seg_idx * samples_per_segment };

    for (size_t i = 0; i < state.segment.size(); ++i) {
      // Segments are centered in frequencies of interest. For example, the
      // first frequenciy of interest is at state.min_freq_index so we must 
      // subtract off (samples_per_segment - 1) / 2 here. For segments 
      // higher in frequency we must add on the total amount of segment samples
      // previously used, i.e. # segment * samples_per_segment.
      size_t lower_segment_index { freq_index - half_samples};

      /// \todo give magic number a name.
      state.segment[i] = power[lower_segment_index + i] * 1.0e40;
    }

    /// \todo create a function to calculate mean and standard deviation together.
    double mean { };
    double standard_dev { };

    for (size_t i = 0; i < samples_per_segment; ++i) {
      double x { std::log(state.segment[i]) };
      mean += x;
      standard_dev += x * x;
    }

    double norm { 1.0 / static_cast<double>(state.segment.size()) };
    mean *= norm;
    mean += constants::euler; /// The log-periodogram. See Wahba 1980.
    standard_dev *= norm;
    standard_dev = std::sqrt(standard_dev - mean * mean); 

    // We wish to calculate the standard deviation of the mean itself.
    double stand_dev_of_mean { standard_dev * std::sqrt(norm) };


    state.means[seg_idx]     = mean;
    state.stand_dev[seg_idx] = stand_dev_of_mean;

    double freq { freq_index / obs_period };
    state.log_freqs[seg_idx] = std::log(freq);

    /// 2 degrees of freedom chi squared PDF median.
    /// The median of log(chi-squared) is log(chisquared median)
    constexpr double chi_square_2dof_median { 2.0 * std::log(2.0) };
    double med { median(state.segment) * chi_square_2dof_median };

    state.log_medians[seg_idx] = std::log(med);


 
  }

}


}
