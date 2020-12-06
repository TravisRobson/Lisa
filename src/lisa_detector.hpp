
#ifndef LISA_DETECTOR_HPP
#define LISA_DETECTOR_HPP


#include <assert.h>
#include <iostream>
#include <string>
#include <vector>
#include <utility>


namespace lisa 
{


// Forward declarations
struct GB_params;
struct XAE;
struct NoiseXAE;


///////////////////////////////////////////////////////////////////////////////
/// \brief  Fill output stream with GB_param values, " ", i.e. space, 
/// delimited.
///
/// \param os output stream. Typically std::cout. Stream to send GB_params.
/// \param params Galactic binary parameters to send to the output stream.
///////////////////////////////////////////////////////////////////////////////
std::ostream& operator<<(std::ostream& os, const GB_params& params);


///////////////////////////////////////////////////////////////////////////////
/// \brief Fill output stream with XAE data, " ", i.e. space, delimited.
///
/// \param os output stream. Typically std::cout, where to print XAE values.
/// \param xae XAE data channels.
///////////////////////////////////////////////////////////////////////////////
std::ostream& operator<<(std::ostream& os, const XAE& xae); 


///////////////////////////////////////////////////////////////////////////////
/// \brief Given a line of space delimited values store them in params.
///
/// The expected order of floating point values space delimited in the 
/// std::string "line" is: freq, freq_dot, theta, phi, amp, iota, psi, and
/// phase. The parameter freq_double_dot is calculated assume GW radiation
/// is fully responsibly for GB evolution. Additionally, theta is assumed to 
/// be in
/// \todo What coordinate system is theta assumed to be in?
/// And therefore theta gets transformed to heliocentric spherical polar 
/// coordinates.
///
/// \params params GB params struct to be filled out.
/// \params line string contain space delimited parameter values.
///////////////////////////////////////////////////////////////////////////////
void fill_gb_params(lisa::GB_params& params, const std::string& line);

void fill_gb_params_binary(GB_params& params, std::istream& is);


///////////////////////////////////////////////////////////////////////////////
/// \brief Verify GB params are within their mathematical and physical range.
///////////////////////////////////////////////////////////////////////////////
void validate(const GB_params& params);


///////////////////////////////////////////////////////////////////////////////
/// \brief Calculate the (analytic) LISA noise PSD at a given frequency.
///
/// \param freq Frequency at which to evaluate this function.
///////////////////////////////////////////////////////////////////////////////
NoiseXAE analytic_instrument_noise(double freq);


NoiseXAE analytic_instrument_noise_legacy(double freq);


///////////////////////////////////////////////////////////////////////////////
/// \brief  Calculate the number of samples in signal.
///
/// i.e. number of real samples for an FFT let's say such that the complex FFT 
/// would be 2x this number of samples.
///
/// \param params GB parameters for which the signal would be generated.
/// \param obs_period Observation period of LISA.
///////////////////////////////////////////////////////////////////////////////
std::pair<size_t, double> gb_bandwidth(const GB_params& params, double obs_period);


// void spacecraft_position(ThreeVec& pos_0, ThreeVec& pos_1, ThreeVec& pos_2, double time);
void spacecraft_position(double* pos_0, double* pos_1, double* pos_2, double time);


struct Strains {

  // Vector<double> y_01;
  // Vector<double> y_02;
  // Vector<double> y_12;
  // Vector<double> y_10;
  // Vector<double> y_20;
  // Vector<double> y_21;

  std::vector<double> y_01;
  std::vector<double> y_02;
  std::vector<double> y_12;
  std::vector<double> y_10;
  std::vector<double> y_20;
  std::vector<double> y_21;

  Strains(int size)
    :
    y_01(size),
    y_02(size),
    y_12(size),
    y_10(size),
    y_20(size),
    y_21(size)
  {}
  void resize(int new_size) {
    y_01.resize(new_size);
    y_02.resize(new_size);
    y_12.resize(new_size);
    y_10.resize(new_size);
    y_20.resize(new_size);
    y_21.resize(new_size);
  }

  int size() const { return y_01.size(); }

};


///////////////////////////////////////////////////////////////////////////////
/// \brief
///
/// \param
///////////////////////////////////////////////////////////////////////////////
void calc_gb_signal(
  const GB_params& params,
  std::vector<double>& X, 
  std::vector<double>& A,
  std::vector<double>& E,
  double obs_period);


}


#endif
