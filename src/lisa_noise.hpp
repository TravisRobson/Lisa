

#ifndef lisa_noise_hpp
#define lisa_noise_hpp


#include <assert.h>
#include <vector>


namespace lisa {



///////////////////////////////////////////////////////////////////////////////
/// This is a vector of strictly real numbers (e.g. power of X TDI channel)
/// which has a sense of it's location in frequency space.
///
/// \todo Maybe I need a better name.
///////////////////////////////////////////////////////////////////////////////
class RealFreqVector {

  std::vector<double> data_;
  size_t min_; ///< Minumum frequency index.
  size_t max_; ///< Maxium frequency index. (inclusive)

public:

  RealFreqVector(size_t min, size_t max) 
    : data_(max - min + 1), min_( min ), max_( max )
  { assert(max > min); }

  RealFreqVector(size_t min, size_t max, double value) 
    : data_(max - min + 1, value), min_( min ), max_( max )
  { assert(max > min); }

  double& operator[](int i) { 
    assert(i > 0 && i <= max_); /// \todo put back to normal indexing after sufficiently tested.
    return data_.at(i - min_); //[i - min_]; 
  }

  const double& operator[](int i) const { 
    assert(i > 0 && i <= max_); /// \todo put back to normal indexing after sufficiently tested.
    return data_.at(i - min_); //[i - min_]; 
  }

  size_t min() const { return min_; }
  size_t max() const { return max_; }

};


///////////////////////////////////////////////////////////////////////////////
/// Currently, the collection of parameters and memory needed for the fit.
///
/// \todo Find a better name for this. Perhaps split it up logically.
///////////////////////////////////////////////////////////////////////////////
struct FitState {
  std::vector<double> instrument_noise; ///< Estimate of the LISA instrument noise.
  std::vector<double> confusion_noise;  ///< Estimate of the GB confusion noise.

  std::vector<double> segment;   ///< Segment of TDI power to take median over.

  std::vector<double> log_freqs;      ///< log frequency for each segment.
  std::vector<double> means;          ///< Mean power in each segment.
  std::vector<double> stand_dev;      ///< Standard deviation of mean power in segment.
  std::vector<double> log_medians;    ///< Median of power in each segment.
  std::vector<double> analytic_noise; ///< Analytic noise estimate in each segment.

  size_t min_freq_index;  ///< Minimum frequency index to perform estimate over. 
  size_t max_freq_index;  ///< Maximum frequency index to perform estimate over.
};


///////////////////////////////////////////////////////////////////////////////
///
///////////////////////////////////////////////////////////////////////////////
void fit_tdi_power(const RealFreqVector& power, FitState& state, double obs_period);


}


#endif

