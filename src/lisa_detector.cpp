

#include "lisa_detector.hpp"

#include <cassert>
#include <cmath>
#include <complex>
#include <sstream>
#include <vector>

#include "lisa_constants.hpp"
#include "lisa_types.hpp"
#include "lisa_utils.hpp" /// Just for Timers (experimenting with global Timers object)


namespace lisa 
{


void fill_gb_params(GB_params& params, const std::string& line)
{

  std::stringstream ss{line};

  ss >> params.freq;
  ss >> params.freq_dot;
  ss >> params.theta;
  ss >> params.phi;
  ss >> params.amp;
  ss >> params.iota;
  ss >> params.psi;
  ss >> params.phase;
  
  // The below equation for freq_double_dot assume gravitation wave emission causes 
  // frequency evolution.
  params.freq_double_dot = 11.0 / 3.0 * params.freq_dot * params.freq_dot / params.freq;
  params.theta = 0.5 * constants::pi - params.theta; // convert coordinates from file.

}


void fill_gb_params_binary(GB_params& params, std::istream& is)
{

  is.read( (char*)&params.freq,     sizeof(params.freq));
  is.read( (char*)&params.freq_dot, sizeof(params.freq_dot));
  is.read( (char*)&params.theta,    sizeof(params.theta));
  is.read( (char*)&params.phi,      sizeof(params.phi));
  is.read( (char*)&params.amp,      sizeof(params.amp));
  is.read( (char*)&params.iota,     sizeof(params.iota));
  is.read( (char*)&params.psi,      sizeof(params.psi));
  is.read( (char*)&params.phase,    sizeof(params.phase));

  // The below equation for freq_double_dot assume gravitation wave emission causes 
  // frequency evolution.
  params.freq_double_dot  = 11.0 / 3.0 * params.freq_dot * params.freq_dot / params.freq;
  params.theta = 0.5 * lisa::constants::pi - params.theta; // convert coordinates from file.

}


///////////////////////////////////////////////////////////////////////////////
/// \todo Need to double check some of the angle parameter domains.
///////////////////////////////////////////////////////////////////////////////
void validate(const GB_params& params) {

  assert(params.freq > 0.0);
  assert(0.0 <= params.theta && params.theta <= constants::pi / 2.0);
  assert(0.0 <= params.phi && params.phi < 2.0 * constants::pi);
  assert(params.amp > 0.0);
  assert(0.0 <= params.iota && params.iota <= constants::pi / 2.0);
  assert(0.0 <= params.psi && params.psi <= constants::pi / 2.0); /// \todo should it pi / 4.0?
  assert(0.0 <= params.phase && params.phase < 2.0 * constants::pi);

}


std::ostream& operator<<(std::ostream& os, const XAE& xae) {
  os << xae.X << " " << xae.A << " " << xae.E;
  return os;
}


std::ostream& operator<<(std::ostream& os, const GB_params& params) {
  os << params.freq << " " << params.freq_dot << " " << params.theta << " " 
     << params.phi << " " << params.amp << " " << params.iota << " " << params.psi
     << " " << params.phase << " " << params.freq_double_dot;
  return os;
}


///////////////////////////////////////////////////////////////////////////////
/// \todo find better names for these variables. I need to reread some 
/// documents
/// \todo I need to look at my up to date code for this. These values must be 
/// out of date.
///
/// I have tried profiling this code with g_timers but it bloats the runtime of
/// this function by a factor of 10! I think I have pinpointed most of that 
/// time to be indexing into the container of Timer's within Timers class.
///
/// https://iopscience.iop.org/article/10.1088/1742-6596/840/1/012024/pdf
///////////////////////////////////////////////////////////////////////////////
NoiseXAE analytic_instrument_noise(double freq)
{

  NoiseXAE result{};

  // Red noise component to the acceleration noise.
  double red_noise { 16.0 * std::pow(2.0e-5 / freq, 10.0) + 1.0e-8 / freq / freq };
  double accel_term { constants::accel_noise / std::pow(2.0 * constants::pi * freq, 4.0) 
    * (1.0 + red_noise) };

  double freq_ratio{ freq / constants::freq_star };
  double cos_freq{ std::cos(freq_ratio) };
  double sin_freq{ std::sin(freq_ratio) };

  double x_accel_factor { 8.0 * (1.0 + cos_freq * cos_freq) };
  result.XYZ = 4.0 * constants::position_noise + x_accel_factor * accel_term;

  double ae_accel_factor { 2.0 * (3.0 + 2.0 * cos_freq + std::cos(2.0 * freq_ratio)) };
  result.AE = 4.0 * (constants::position_noise * (2.0 + cos_freq) + ae_accel_factor * accel_term);

  double transfer_func { 4.0 * sin_freq * sin_freq };
  double factor{ transfer_func / (4.0 * constants::arm_length * constants::arm_length) };

  result.XYZ *= factor;
  result.AE  *= factor;

  return result;

}


NoiseXAE analytic_instrument_noise_legacy(double freq)
{

  double inv_freq { 1.0 / freq };
  double red_noise { 16.0 * (std::pow(2.0e-5 * inv_freq, 10.0) + std::pow(1.e-4 * inv_freq, 2.0)) };

  constexpr double loc_noise { 2.890e-24 };
  constexpr double pos_noise { 8.321e-23 };

  double arm_factor { 1.0 / std::pow(2.0 * constants::arm_length, 2.0) };
  double freq_ratio { freq / constants::freq_star };
  double cos_freq{ std::cos(freq_ratio) };
  double sin_freq{ std::sin(freq_ratio) };
  double transfer_func { 4.0  * std::pow(sin_freq, 2.0) * arm_factor };

  NoiseXAE result;

  double accel_term { constants::accel_noise * (1.0 + red_noise) / 
    std::pow(2.0 * constants::pi * freq, 4.0)  };

  constexpr double xyz_1 { 4.0 * (loc_noise + pos_noise) };
  double xyz_2 { 8.0 * (1.0 + std::pow(cos_freq, 2.0)) };
  result.XYZ = transfer_func * (xyz_1 + xyz_2 * (0.5 * loc_noise + accel_term));

  // terms broken up for clarity.
  double ae_1 { (2.0 + cos_freq) * (loc_noise + pos_noise) };
  double ae_2 { 3.0 + 2.0 * cos_freq + std::cos(2.0 * freq_ratio) };

  result.AE = 4.0 / 3.0 * transfer_func * (ae_1 + 2.0 * ae_2 * (0.5 * loc_noise + accel_term));

  return result;

}



///////////////////////////////////////////////////////////////////////////////
/// \todo  This is a pretty out of date function. I am currently use it while
/// I reproduce the functionality of the original code. 
///////////////////////////////////////////////////////////////////////////////
std::pair<int, double> gb_bandwidth(const GB_params& params, double obs_period)
{

  //
  // Make an adjustment based on how many years LISA observes.
  //
  int factor { 8 };
  double years { obs_period / constants::year };
  if (years <= 4.0) factor = 4;
  if (years <= 2.0) factor = 2;
  if (years <= 1.0) factor = 1;

  int samples { 32 * factor };

  //
  // Make an adjustment for carrier frequency of the GB. Higher frequency
  // sources experiece more frequency evolution and more Doppler shift
  // as LISA orbits the Sun.
  //
  if (params.freq > 0.001) samples = 64 * factor;
  if (params.freq > 0.01)  samples = 256 * factor;
  if (params.freq > 0.03)  samples = 512 * factor;
  if (params.freq > 0.1)   samples = 1024 * factor;

  double freq_ratio { params.freq / constants::freq_star };

  // NoiseXAE noise = analytic_instrument_noise(params.freq);
  NoiseXAE noise = analytic_instrument_noise_legacy(params.freq);

  double transfer_func { 4.0 * std::pow(std::sin(freq_ratio), 2.0) };
  double michelson_noise { noise.XYZ / transfer_func };

  double cutoff_amp { params.amp * std::sqrt(obs_period / michelson_noise) };

  int new_samples { static_cast<int>(std::pow(2.0, std::rint(std::log2(cutoff_amp)) + 1)) };

  // I can't say I understand this logic. YES, they're if's
  if (new_samples < samples) new_samples = samples;
  if (samples < new_samples) samples = new_samples;
  if (new_samples > 8192) new_samples = 8192; 

  samples = new_samples;

  return std::pair<int, double>(samples, cutoff_amp);

}

void spacecraft_position(ThreeVec& pos_0, ThreeVec& pos_1, ThreeVec& pos_2, double time)
{

  constexpr double beta_0 { 0.0 + constants::lambda };
  constexpr double beta_1 { 2.0 * constants::pi / 3.0 + constants::lambda };
  constexpr double beta_2 { 4.0 * constants::pi / 3.0 + constants::lambda };

  constexpr double sin_beta_0 { std::sin(beta_0) };
  constexpr double cos_beta_0 { std::cos(beta_0) };

  constexpr double sin_beta_1 { std::sin(beta_1) };
  constexpr double cos_beta_1 { std::cos(beta_1) };

  constexpr double sin_beta_2 { std::sin(beta_2) };
  constexpr double cos_beta_2 { std::cos(beta_2) };

  double phase { 2.0 * constants::pi * constants::mod_freq * time + constants::kappa };

  double sin_phase { std::sin(phase) };
  double cos_phase { std::cos(phase) };

  // X coordinate
  pos_0[0] = cos_phase + constants::eccentricity * 
    (sin_phase * cos_phase * sin_beta_0 - (1.0 + sin_phase * sin_phase) * cos_beta_0);

  pos_1[0] = cos_phase + constants::eccentricity * 
    (sin_phase * cos_phase * sin_beta_1 - (1.0 + sin_phase * sin_phase) * cos_beta_1);

  pos_2[0] = cos_phase + constants::eccentricity * 
    (sin_phase * cos_phase * sin_beta_2 - (1.0 + sin_phase * sin_phase) * cos_beta_2);

  // Y coordinate
  pos_0[1] = sin_phase + constants::eccentricity *
    (sin_phase * cos_phase * cos_beta_0 - (1.0 + cos_phase * cos_phase) * sin_beta_0);

  pos_1[1] = sin_phase + constants::eccentricity *
    (sin_phase * cos_phase * cos_beta_1 - (1.0 + cos_phase * cos_phase) * sin_beta_1);

  pos_2[1] = sin_phase + constants::eccentricity *
    (sin_phase * cos_phase * cos_beta_2 - (1.0 + cos_phase * cos_phase) * sin_beta_2);

  // Z coordinate
  pos_0[2] = -std::sqrt(3.0) * constants::eccentricity *
    (cos_phase * cos_beta_0 + sin_phase * sin_beta_0);

  pos_1[2] = -std::sqrt(3.0) * constants::eccentricity *
    (cos_phase * cos_beta_1 + sin_phase * sin_beta_1);

  pos_2[2] = -std::sqrt(3.0) * constants::eccentricity *
    (cos_phase * cos_beta_2 + sin_phase * sin_beta_2);

  // Multiply each term by 1.0 AU
  pos_0 *= constants::au;
  pos_1 *= constants::au;
  pos_2 *= constants::au;
   
}


void spacecraft_position(double* pos_0, double* pos_1, double* pos_2, double time) 
{
  constexpr double beta_0 { 0.0 + constants::lambda };
  constexpr double beta_1 { 2.0 * constants::pi / 3.0 + constants::lambda };
  constexpr double beta_2 { 4.0 * constants::pi / 3.0 + constants::lambda };

  constexpr double sin_beta_0 { std::sin(beta_0) };
  constexpr double cos_beta_0 { std::cos(beta_0) };

  constexpr double sin_beta_1 { std::sin(beta_1) };
  constexpr double cos_beta_1 { std::cos(beta_1) };

  constexpr double sin_beta_2 { std::sin(beta_2) };
  constexpr double cos_beta_2 { std::cos(beta_2) };

  double phase { 2.0 * constants::pi * constants::mod_freq * time + constants::kappa };

  double cos_phase { std::cos(phase) };
  double sin_phase { std::sin(phase) };

  // X coordinate
  pos_0[0] = cos_phase + constants::eccentricity * 
    (sin_phase * cos_phase * sin_beta_0 - (1.0 + sin_phase * sin_phase) * cos_beta_0);

  pos_1[0] = cos_phase + constants::eccentricity * 
    (sin_phase * cos_phase * sin_beta_1 - (1.0 + sin_phase * sin_phase) * cos_beta_1);

  pos_2[0] = cos_phase + constants::eccentricity * 
    (sin_phase * cos_phase * sin_beta_2 - (1.0 + sin_phase * sin_phase) * cos_beta_2);

  // Y coordinate
  pos_0[1] = sin_phase + constants::eccentricity *
    (sin_phase * cos_phase * cos_beta_0 - (1.0 + cos_phase * cos_phase) * sin_beta_0);

  pos_1[1] = sin_phase + constants::eccentricity *
    (sin_phase * cos_phase * cos_beta_1 - (1.0 + cos_phase * cos_phase) * sin_beta_1);

  pos_2[1] = sin_phase + constants::eccentricity *
    (sin_phase * cos_phase * cos_beta_2 - (1.0 + cos_phase * cos_phase) * sin_beta_2);

  // Z coordinate
  pos_0[2] = -std::sqrt(3.0) * constants::eccentricity *
    (cos_phase * cos_beta_0 + sin_phase * sin_beta_0);

  pos_1[2] = -std::sqrt(3.0) * constants::eccentricity *
    (cos_phase * cos_beta_1 + sin_phase * sin_beta_1);

  pos_2[2] = -std::sqrt(3.0) * constants::eccentricity *
    (cos_phase * cos_beta_2 + sin_phase * sin_beta_2);

  // Multiply each term by 1.0 AU
  pos_0[0] *= constants::au;
  pos_0[1] *= constants::au;
  pos_0[2] *= constants::au;

  pos_1[0] *= constants::au;
  pos_1[1] *= constants::au;
  pos_1[2] *= constants::au;

  pos_2[0] *= constants::au;  
  pos_2[1] *= constants::au;  
  pos_2[2] *= constants::au;  

}


// std::complex<double> transfer(
//   int first, 
//   int second, 
//   const GB_params& params,
//   double* freq_ratio,
//   double k_dot_r[][3],
//   double angular_carrier_freq,
//   double time,
//   double plus_ij[][3],
//   double cross_ij[][3],
//   double amp_plus_real,
//   double amp_plus_imag,
//   double amp_cross_real,
//   double amp_cross_imag,
//   double* xi)
// {


//   double transfer_arg = 0.5 * freq_ratio[first] * (1.0 - k_dot_r[first][second]);

//   double exp_arg = 2.0 * params.freq * xi[first] - angular_carrier_freq * time +
//                    params.phase + params.freq_dot * xi[first] * xi[first] +
//                    params.freq_double_dot * xi[first] * xi[first] * xi[first] / 3.0;
//   exp_arg *= constants::pi;

//   double sinc = 0.25 * std::sin(transfer_arg) / transfer_arg;

//   // Evolution of amplitude as GW frequency increases.
//   double amp_evolve = 1.0 + 2.0 / 3.0 * params.freq_dot / params.freq * xi[first];

//   double term_1_real = plus_ij[first][second] * amp_plus_real + cross_ij[first][second] * amp_cross_real;
//   double term_1_imag = plus_ij[first][second] * amp_plus_imag + cross_ij[first][second] * amp_cross_imag;

//   double term_2_real = std::cos(transfer_arg + exp_arg);
//   double term_2_imag = std::sin(transfer_arg + exp_arg);

  // std::complex<double> result { };

//   std::complex<double> result {
//     term_1_real * term_2_real - term_1_imag * term_2_imag,
//     term_1_real * term_2_imag + term_1_imag * term_2_real
//   };


//   // // signal_real[i][j] = term_1_real * term_2_real - term_1_imag * term_2_imag;
//   // // signal_imag[i][j] = term_1_real * term_2_imag + term_1_imag * term_2_real;

//   // result *= amp_evolve * sinc;
// /*
//   signal_real[i][j] *= amp_evolve * sinc;
//   signal_imag[i][j] *= amp_evolve * sinc;*/
//   return result;
// }


void calc_gb_signal(
  const GB_params& params,
  std::vector<double>& X, 
  std::vector<double>& A,
  std::vector<double>& E,
  double obs_period)
{

  // Data channels must be of the same size
  assert(X.size() == A.size());
  assert(X.size() == E.size());
  // assert(X.size() == strains_buffer.size());
  assert(obs_period > 0.0);

  // Calc GW polarization amplitudes
  double sin_psi { std::sin(2.0 * params.psi) };
  double cos_psi { std::cos(2.0 * params.psi) };

  double cos_iota { std::cos(params.iota) };

  double amp_plus { params.amp * (1.0 + cos_iota * cos_iota) };
  double amp_cross { -2.0 * params.amp * cos_iota };
  // Calculate the constant pieces of transfer functions
  double amp_plus_real { amp_plus * cos_psi }; // Real part of complex amplitude
  double amp_plus_imag { -amp_cross * sin_psi }; // Imaginary part (imag)
  double amp_cross_real { -amp_plus * sin_psi };
  double amp_cross_imag { -amp_cross * cos_psi };

  // Calculate thee tensors for the slowly evolving response.
  //  First, calculate the basis vectors the the GW.
  double sin_theta { std::sin(params.theta) };
  double cos_theta { std::cos(params.theta) };

  double sin_phi { std::sin(params.phi) };
  double cos_phi { std::cos(params.phi) };

  // ThreeVec u { cos_theta * cos_phi, cos_theta * sin_phi, -sin_theta }; // u vector
  // ThreeVec v { sin_phi, -cos_phi, 0.0 }; // v vector
  // ThreeVec k { -sin_theta * cos_phi, -sin_theta * sin_phi, -cos_theta }; // GW propagation vector

  double u[3] { cos_theta * cos_phi, cos_theta * sin_phi, -sin_theta }; // u vector
  double v[3] { sin_phi, -cos_phi, 0.0 }; // v vector
  double k[3] { -sin_theta * cos_phi, -sin_theta * sin_phi, -cos_theta }; // GW propagation vector


  // // The polarization tensors.
  // ThreeTensor e_plus { };
  // e_plus(0, 0) = u[0] * u[0] - v[0] * v[0];
  // e_plus(0, 1) = u[0] * u[1] - v[0] * v[1];
  // e_plus(0, 2) = u[0] * u[2] - v[0] * v[2];

  // e_plus(1, 0) = u[1] * u[0] - v[1] * v[0];
  // e_plus(1, 1) = u[1] * u[1] - v[1] * v[1];
  // e_plus(1, 2) = u[1] * u[2] - v[1] * v[2];

  // e_plus(2, 0) = u[2] * u[0] - v[2] * v[0];
  // e_plus(2, 1) = u[2] * u[1] - v[2] * v[1];
  // e_plus(2, 2) = u[2] * u[2] - v[2] * v[2];

  // ThreeTensor e_cross { };
  // e_cross(0, 0) = u[0] * v[0] + v[0] * u[0];
  // e_cross(0, 1) = u[0] * v[1] + v[0] * u[1];
  // e_cross(0, 2) = u[0] * v[2] + v[0] * u[2];

  // e_cross(1, 0) = u[1] * v[0] + v[1] * u[0];
  // e_cross(1, 1) = u[1] * v[1] + v[1] * u[1];
  // e_cross(1, 2) = u[1] * v[2] + v[1] * u[2];

  // e_cross(2, 0) = u[2] * v[0] + v[2] * u[0];
  // e_cross(2, 1) = u[2] * v[1] + v[2] * u[1];
  // e_cross(2, 2) = u[2] * v[2] + v[2] * u[2];

  // The polarization tensors.
  double e_plus[3][3] { };
  // e_plus[0][0] = u[0] * u[0] - v[0] * v[0];
  // e_plus[0][1] = u[0] * u[1] - v[0] * v[1];
  // e_plus[0][2] = u[0] * u[2] - v[0] * v[2];

  // e_plus[1][0] = u[1] * u[0] - v[1] * v[0];
  // e_plus[1][1] = u[1] * u[1] - v[1] * v[1];
  // e_plus[1][2] = u[1] * u[2] - v[1] * v[2];

  // e_plus[2][0] = u[2] * u[0] - v[2] * v[0];
  // e_plus[2][1] = u[2] * u[1] - v[2] * v[1];
  // e_plus[2][2] = u[2] * u[2] - v[2] * v[2];

  double e_cross[3][3] { };
  // e_cross[0][0] = u[0] * v[0] + v[0] * u[0];
  // e_cross[0][1] = u[0] * v[1] + v[0] * u[1];
  // e_cross[0][2] = u[0] * v[2] + v[0] * u[2];

  // e_cross[1][0] = u[1] * v[0] + v[1] * u[0];
  // e_cross[1][1] = u[1] * v[1] + v[1] * u[1];
  // e_cross[1][2] = u[1] * v[2] + v[1] * u[2];

  // e_cross[2][0] = u[2] * v[0] + v[2] * u[0];
  // e_cross[2][1] = u[2] * v[1] + v[2] * u[1];
  // e_cross[2][2] = u[2] * v[2] + v[2] * u[2];


  for(int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      e_plus [i][j] = u[i] * u[j] - v[i] * v[j];
      e_cross[i][j] = u[i] * v[j] + v[i] * u[j];
    } 
  }


  //
  // Calucate the response for all times.
  //  

  // Calculate frequency bin of carrier frequency
  int carrier_freq_bin { static_cast<int>(params.freq * obs_period) };

  double angular_carrier_freq { 2.0 * constants::pi * 
    static_cast<double>(carrier_freq_bin) / obs_period };

  int num_time_samples { X.size() / 2 }; // because X holds real and imag values.

  std::vector<double> signal_01(X.size());
  std::vector<double> signal_02(X.size());
  std::vector<double> signal_12(X.size());
  std::vector<double> signal_10(X.size());
  std::vector<double> signal_20(X.size());
  std::vector<double> signal_21(X.size());

  // Vector<double> signal_01(X.size());
  // Vector<double> signal_02(X.size());
  // Vector<double> signal_12(X.size());
  // Vector<double> signal_10(X.size());
  // Vector<double> signal_20(X.size());
  // Vector<double> signal_21(X.size());

  // double* signal_01 = new double[X.size()];
  // double* signal_02 = new double[X.size()];
  // double* signal_12 = new double[X.size()];
  // double* signal_10 = new double[X.size()];
  // double* signal_20 = new double[X.size()];
  // double* signal_21 = new double[X.size()];

  // Strains strains(X.size());

  for (int time_idx = 0; time_idx < num_time_samples; ++time_idx) {

    // First time sample must be at 0.0 to have correct phasing.
    double time { obs_period * static_cast<double>(time_idx) 
      / static_cast<double>(num_time_samples) };

    // // Calucate the barycentric position of LISA spacecraft.
    // ThreeVec pos_0 {};
    // ThreeVec pos_1 {};
    // ThreeVec pos_2 {};
    // spacecraft_position(pos_0, pos_1, pos_2, time);

    // Calucate the barycentric position of LISA spacecraft.
    double pos_0[3] {};
    double pos_1[3] {};
    double pos_2[3] {};
    spacecraft_position(pos_0, pos_1, pos_2, time);


    // // Calculate unit vectors separating spacecraft
    // ThreeVec r_01 { (pos_1 - pos_0) / constants::arm_length };
    // ThreeVec r_02 { (pos_2 - pos_0) / constants::arm_length };
    // ThreeVec r_12 { (pos_2 - pos_1) / constants::arm_length };

    // Calculate unit vectors separating spacecraft
    double r_01[3] = {
      (pos_1[0] - pos_0[0]) / constants::arm_length,
      (pos_1[1] - pos_0[1]) / constants::arm_length,
      (pos_1[2] - pos_0[2]) / constants::arm_length
    };
    double r_02[3] = {
      (pos_2[0] - pos_0[0]) / constants::arm_length,
      (pos_2[1] - pos_0[1]) / constants::arm_length,
      (pos_2[2] - pos_0[2]) / constants::arm_length
    };
    double r_12[3] = {
      (pos_2[0] - pos_1[0]) / constants::arm_length,
      (pos_2[1] - pos_1[1]) / constants::arm_length,
      (pos_2[2] - pos_1[2]) / constants::arm_length
    };

    // Indices into these represent spacecraft indices not spatial indices;
    // hence why we didn't use ThreeTensor.
    double plus_ij[3][3]  { };
    double cross_ij[3][3] { };

    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        plus_ij [0][1] += e_plus [i][j] * r_01[i] * r_01[j];
        cross_ij[0][1] += e_cross[i][j] * r_01[i] * r_01[j];

        plus_ij [0][2] += e_plus [i][j] * r_02[i] * r_02[j];
        cross_ij[0][2] += e_cross[i][j] * r_02[i] * r_02[j];

        plus_ij [1][2] += e_plus [i][j] * r_12[i] * r_12[j];
        cross_ij[1][2] += e_cross[i][j] * r_12[i] * r_12[j];

      }
    }

    // plus_ij [0][1] = dot(e_plus,  r_01, r_01);
    // cross_ij[0][1] = dot(e_cross, r_01, r_01);

    // plus_ij [0][2] = dot(e_plus,  r_02, r_02);
    // cross_ij[0][2] = dot(e_cross, r_02, r_02);

    // plus_ij [1][2] = dot(e_plus,  r_12, r_12);
    // cross_ij[1][2] = dot(e_cross, r_12, r_12);


    // Use symmetry to cheaply calculate other components of this matrix
    plus_ij [1][0] = plus_ij [0][1];
    cross_ij[1][0] = cross_ij[0][1];

    plus_ij [2][0] = plus_ij [0][2];
    cross_ij[2][0] = cross_ij[0][2];

    plus_ij [2][1] = plus_ij [1][2];
    cross_ij[2][1] = cross_ij[1][2];


    // Matrix of dot products between spacecraft separation vectors and 
    // GW propagation vector.
    double k_dot_r[3][3] { };

    for (int i = 0; i < 3; ++i) {
      k_dot_r[0][1] += k[i] * r_01[i];
      k_dot_r[0][2] += k[i] * r_02[i];
      k_dot_r[1][2] += k[i] * r_12[i];
    }


    // k_dot_r[0][1] = dot(k, r_01);
    // k_dot_r[0][2] = dot(k, r_02);
    // k_dot_r[1][2] = dot(k, r_12);

    // Use antisymmetry to calculate the other elements.
    k_dot_r[1][0] = -k_dot_r[0][1];
    k_dot_r[2][0] = -k_dot_r[0][2];
    k_dot_r[2][1] = -k_dot_r[1][2];


    // GW wave vector dotted with spacecraft position (x).
    // double k_dot_x[3] { 
    //   dot(k, pos_0) / constants::c_light,
    //   dot(k, pos_1) / constants::c_light,
    //   dot(k, pos_2) / constants::c_light
    // }; 

    double k_dot_x[3] { 
      (k[0] * pos_0[0] + k[1] * pos_0[1] + k[2] * pos_0[2]) / constants::c_light,
      (k[0] * pos_1[0] + k[1] * pos_1[1] + k[2] * pos_1[2]) / constants::c_light,
      (k[0] * pos_2[0] + k[1] * pos_2[1] + k[2] * pos_2[2]) / constants::c_light
    }; 


    // Redshifted wave variable.
    double xi[3] { 
      time - k_dot_x[0],
      time - k_dot_x[1],
      time - k_dot_x[2]
    };

    // Redshifted frequency
    double freq_at_spacecraft[3] { 
      params.freq + params.freq_dot * xi[0] + 0.5 * params.freq_double_dot * xi[0] * xi[0],
      params.freq + params.freq_dot * xi[1] + 0.5 * params.freq_double_dot * xi[1] * xi[1],
      params.freq + params.freq_dot * xi[2] + 0.5 * params.freq_double_dot * xi[2] * xi[2]
    };


    double freq_ratio[3] {
      freq_at_spacecraft[0] / constants::freq_star,
      freq_at_spacecraft[1] / constants::freq_star,
      freq_at_spacecraft[2] / constants::freq_star
    };

    double signal_real[3][3] { };
    double signal_imag[3][3] { };

    // double signal_real[6] { };
    // double signal_imag[6] { };

    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {

        if (i != j) {

          double term_1_real = plus_ij[i][j] * amp_plus_real + cross_ij[i][j] * amp_cross_real;
          double term_1_imag = plus_ij[i][j] * amp_plus_imag + cross_ij[i][j] * amp_cross_imag;

          double transfer_arg = 0.5 * freq_ratio[i] * (1.0 - k_dot_r[i][j]);

          double exp_arg = 2.0 * params.freq * xi[i]  +
                           params.freq_dot * xi[i] * xi[i] +
                           params.freq_double_dot * xi[i] * xi[i] * xi[i] / 3.0;
          exp_arg *= constants::pi;
          exp_arg += params.phase - angular_carrier_freq * time;

          double term_2_real = std::cos(transfer_arg + exp_arg);
          double term_2_imag = std::sin(transfer_arg + exp_arg);    

          signal_real[i][j] = term_1_real * term_2_real - term_1_imag * term_2_imag;
          signal_imag[i][j] = term_1_real * term_2_imag + term_1_imag * term_2_real;

            // Evolution of amplitude as GW frequency increases.
          double amp_evolve = 1.0 + 2.0 / 3.0 * params.freq_dot / params.freq * xi[i];
          double sinc       = 0.25 * std::sin(transfer_arg) / transfer_arg;

          signal_real[i][j] *= amp_evolve * sinc;
          signal_imag[i][j] *= amp_evolve * sinc;

        }
      }
    }



//   // std::complex<double> result { };

//   std::complex<double> result {
//     term_1_real * term_2_real - term_1_imag * term_2_imag,
//     term_1_real * term_2_imag + term_1_imag * term_2_real
//   };


//   // // signal_real[i][j] = term_1_real * term_2_real - term_1_imag * term_2_imag;
//   // // signal_imag[i][j] = term_1_real * term_2_imag + term_1_imag * term_2_real;

//   // result *= amp_evolve * sinc;
// /*
//   signal_real[i][j] *= amp_evolve * sinc;
//   signal_imag[i][j] *= amp_evolve * sinc;*/


    // std::complex<double> out = transfer(
    //                                     0, 
    //                                     1, 
    //                                     params,
    //                                     freq_ratio,
    //                                     k_dot_r,
    //                                     angular_carrier_freq,
    //                                     time,
    //                                     plus_ij,
    //                                     cross_ij,
    //                                     amp_plus_real,
    //                                     amp_plus_imag,
    //                                     amp_cross_real,
    //                                     amp_cross_imag,
    //                                     xi);

    // signal_real[0] = out.real();
    // signal_imag[0] = out.imag();


    // out = transfer(
    //                                     0, 
    //                                     2, 
    //                                     params,
    //                                     freq_ratio,
    //                                     k_dot_r,
    //                                     angular_carrier_freq,
    //                                     time,
    //                                     plus_ij,
    //                                     cross_ij,
    //                                     amp_plus_real,
    //                                     amp_plus_imag,
    //                                     amp_cross_real,
    //                                     amp_cross_imag,
    //                                     xi);

    // signal_real[1] = out.real();
    // signal_imag[1] = out.imag();

    // out = transfer(
    //                                     1, 
    //                                     2, 
    //                                     params,
    //                                     freq_ratio,
    //                                     k_dot_r,
    //                                     angular_carrier_freq,
    //                                     time,
    //                                     plus_ij,
    //                                     cross_ij,
    //                                     amp_plus_real,
    //                                     amp_plus_imag,
    //                                     amp_cross_real,
    //                                     amp_cross_imag,
    //                                     xi);

    // signal_real[2] = out.real();
    // signal_imag[2] = out.imag();

    // out = transfer(
    //                                     1, 
    //                                     0, 
    //                                     params,
    //                                     freq_ratio,
    //                                     k_dot_r,
    //                                     angular_carrier_freq,
    //                                     time,
    //                                     plus_ij,
    //                                     cross_ij,
    //                                     amp_plus_real,
    //                                     amp_plus_imag,
    //                                     amp_cross_real,
    //                                     amp_cross_imag,
    //                                     xi);

    // signal_real[3] = out.real();
    // signal_imag[3] = out.imag();

    // out = transfer(
    //                                     2, 
    //                                     0, 
    //                                     params,
    //                                     freq_ratio,
    //                                     k_dot_r,
    //                                     angular_carrier_freq,
    //                                     time,
    //                                     plus_ij,
    //                                     cross_ij,
    //                                     amp_plus_real,
    //                                     amp_plus_imag,
    //                                     amp_cross_real,
    //                                     amp_cross_imag,
    //                                     xi);

    // signal_real[4] = out.real();
    // signal_imag[4] = out.imag();


    // out = transfer(
    //                                     2, 
    //                                     1, 
    //                                     params,
    //                                     freq_ratio,
    //                                     k_dot_r,
    //                                     angular_carrier_freq,
    //                                     time,
    //                                     plus_ij,
    //                                     cross_ij,
    //                                     amp_plus_real,
    //                                     amp_plus_imag,
    //                                     amp_cross_real,
    //                                     amp_cross_imag,
    //                                     xi);

    // signal_real[5] = out.real();
    // signal_imag[5] = out.imag();

    // Fill in time series with slowly evolving signal
    int real_idx { 2 * time_idx };
    int imag_idx { 2 * time_idx + 1 };

    // signal_01[real_idx] = signal_real[0];
    // signal_01[imag_idx] = signal_imag[0];

    // signal_02[real_idx] = signal_real[1];
    // signal_02[imag_idx] = signal_imag[1];

    // signal_12[real_idx] = signal_real[2];
    // signal_12[imag_idx] = signal_imag[2];

    // signal_10[real_idx] = signal_real[3];
    // signal_10[imag_idx] = signal_imag[3];

    // signal_20[real_idx] = signal_real[4];
    // signal_20[imag_idx] = signal_imag[4];

    // signal_21[real_idx] = signal_real[5];
    // signal_21[imag_idx] = signal_imag[5];

    signal_01[real_idx] = signal_real[0][1];
    signal_01[imag_idx] = signal_imag[0][1];

    signal_02[real_idx] = signal_real[0][2];
    signal_02[imag_idx] = signal_imag[0][2];

    signal_12[real_idx] = signal_real[1][2];
    signal_12[imag_idx] = signal_imag[1][2];

    signal_10[real_idx] = signal_real[1][0];
    signal_10[imag_idx] = signal_imag[1][0];

    signal_20[real_idx] = signal_real[2][0];
    signal_20[imag_idx] = signal_imag[2][0];

    signal_21[real_idx] = signal_real[2][1]; 
    signal_21[imag_idx] = signal_imag[2][1];


    // printf("signal_01[real_idx] %g\n", signal_01[real_idx]);

    // strains_buffer.y_01[real_idx] = signal_real[0];
    // strains_buffer.y_01[imag_idx] = signal_imag[0];

    // strains_buffer.y_02[real_idx] = signal_real[1];
    // strains_buffer.y_02[imag_idx] = signal_imag[1];

    // strains_buffer.y_12[real_idx] = signal_real[2];
    // strains_buffer.y_12[imag_idx] = signal_imag[2];

    // strains_buffer.y_10[real_idx] = signal_real[3];
    // strains_buffer.y_10[imag_idx] = signal_imag[3];

    // strains_buffer.y_20[real_idx] = signal_real[4];
    // strains_buffer.y_20[imag_idx] = signal_imag[4];

    // strains_buffer.y_21[real_idx] = signal_real[5];
    // strains_buffer.y_21[imag_idx] = signal_imag[5];

    // strains_buffer.y_01[real_idx] = signal_real[0][1];
    // strains_buffer.y_01[imag_idx] = signal_imag[0][1];

    // strains_buffer.y_02[real_idx] = signal_real[0][2];
    // strains_buffer.y_02[imag_idx] = signal_imag[0][2];

    // strains_buffer.y_12[real_idx] = signal_real[1][2];
    // strains_buffer.y_12[imag_idx] = signal_imag[1][2];

    // strains_buffer.y_10[real_idx] = signal_real[1][0];
    // strains_buffer.y_10[imag_idx] = signal_imag[1][0];

    // strains_buffer.y_20[real_idx] = signal_real[2][0];
    // strains_buffer.y_20[imag_idx] = signal_imag[2][0];

    // strains_buffer.y_21[real_idx] = signal_real[2][1];
    // strains_buffer.y_21[imag_idx] = signal_imag[2][1];

    // strains.y_01[real_idx] = signal_real[0];
    // strains.y_01[imag_idx] = signal_imag[0];

    // strains.y_02[real_idx] = signal_real[1];
    // strains.y_02[imag_idx] = signal_imag[1];

    // strains.y_12[real_idx] = signal_real[2];
    // strains.y_12[imag_idx] = signal_imag[2];

    // strains.y_10[real_idx] = signal_real[3];
    // strains.y_10[imag_idx] = signal_imag[3];

    // strains.y_20[real_idx] = signal_real[4];
    // strains.y_20[imag_idx] = signal_imag[4];

    // strains.y_21[real_idx] = signal_real[5];
    // strains.y_21[imag_idx] = signal_imag[5];

  } // end time_idx for loop

  // Numerical Fourier transform of slowly evolving signal 
  int status = gsl_fft_complex_radix2_forward(signal_01.data(), 1, num_time_samples);
  if (status) { throw std::runtime_error("GSL error"); }

  status = gsl_fft_complex_radix2_forward(signal_02.data(), 1, num_time_samples);
  if (status) { throw std::runtime_error("GSL error"); }

  status = gsl_fft_complex_radix2_forward(signal_12.data(), 1, num_time_samples);
  if (status) { throw std::runtime_error("GSL error"); }

  status = gsl_fft_complex_radix2_forward(signal_10.data(), 1, num_time_samples);
  if (status) { throw std::runtime_error("GSL error"); }

  status = gsl_fft_complex_radix2_forward(signal_20.data(), 1, num_time_samples);
  if (status) { throw std::runtime_error("GSL error"); }

  status = gsl_fft_complex_radix2_forward(signal_21.data(), 1, num_time_samples);
  if (status) { throw std::runtime_error("GSL error"); }


  std::vector<double> y_01(X.size());
  std::vector<double> y_02(X.size());
  std::vector<double> y_12(X.size());
  std::vector<double> y_10(X.size());
  std::vector<double> y_20(X.size());
  std::vector<double> y_21(X.size());

  // The 0.5 is to make the time series real.
  double norm = 1.0 / static_cast<double>(num_time_samples) * 0.5;

  // "Unwrap" the FFT.
  // Recall FFT frequencies organized like:
  // neg[N/2-1], neg[N/2-2], ... , pos[0], ..., pos[N/2] 
  // Move negative frequncy part to "right" of positive.
  for (int i = 0; i < num_time_samples; ++i) {

    y_01[i] = signal_01[num_time_samples + i] * norm;
    y_02[i] = signal_02[num_time_samples + i] * norm;
    y_12[i] = signal_12[num_time_samples + i] * norm;
    y_10[i] = signal_10[num_time_samples + i] * norm;
    y_20[i] = signal_20[num_time_samples + i] * norm;
    y_21[i] = signal_21[num_time_samples + i] * norm;

    y_01[i + num_time_samples] = signal_01[i] * norm;
    y_02[i + num_time_samples] = signal_02[i] * norm;
    y_12[i + num_time_samples] = signal_12[i] * norm;
    y_10[i + num_time_samples] = signal_10[i] * norm;
    y_20[i + num_time_samples] = signal_20[i] * norm;
    y_21[i + num_time_samples] = signal_21[i] * norm;

  }



  double phi_lisa_simulator { 2.0 * constants::pi * params.freq * 
    (0.5 * constants::delta_t - constants::arm_length / constants::c_light) };

  double sin_offset { std::sin(phi_lisa_simulator) };
  double cos_offset { std::cos(phi_lisa_simulator) };

  std::vector<double> X_temp(X.size());
  std::vector<double> Y_temp(X.size());
  std::vector<double> Z_temp(X.size());
  std::vector<double> Y_lisa_sim(X.size());
  std::vector<double> Z_lisa_sim(X.size());


  for (int i = 0; i < num_time_samples; ++i) {

    double freq { static_cast<double>(carrier_freq_bin + i - num_time_samples / 2) / obs_period };
    double freq_on_fstar { freq / constants::freq_star };

    double sin_1 { std::sin(1.0 * freq_on_fstar) };
    double cos_1 { std::cos(1.0 * freq_on_fstar) };

    double sin_2 { std::sin(2.0 * freq_on_fstar) };
    double cos_2 { std::cos(2.0 * freq_on_fstar) };

    double sin_3 { std::sin(3.0 * freq_on_fstar) };
    double cos_3 { std::cos(3.0 * freq_on_fstar) };

    int real { 2 * i };
    int imag { 2 * i + 1 };

    // X channel
    X_temp[real] = (y_01[real] - y_02[real]) * cos_3 + (y_01[imag] - y_02[imag]) * sin_3 +
                   (y_10[real] - y_20[real]) * cos_2 + (y_10[imag] - y_20[imag]) * sin_2 +
                   (y_02[real] - y_01[real]) * cos_1 + (y_02[imag] - y_01[imag]) * sin_1 +
                   (y_20[real] - y_10[real]);

    X_temp[imag] = (y_01[imag] - y_02[imag]) * cos_3 - (y_01[real] - y_02[real]) * sin_3 +
                   (y_10[imag] - y_20[imag]) * cos_2 - (y_10[real] - y_20[real]) * sin_2 +
                   (y_02[imag] - y_01[imag]) * cos_1 - (y_02[real] - y_01[real]) * sin_1 +
                   (y_20[imag] - y_10[imag]);

    // Y channel
    Y_temp[real] = (y_12[real] - y_10[real]) * cos_3 + (y_12[imag] - y_10[imag]) * sin_3 +
                   (y_21[real] - y_01[real]) * cos_2 + (y_21[imag] - y_01[imag]) * sin_2 +
                   (y_10[real] - y_12[real]) * cos_1 + (y_10[imag] - y_12[imag]) * sin_1 +
                   (y_01[real] - y_21[real]);

    Y_temp[imag] = (y_12[imag] - y_10[imag]) * cos_3 - (y_12[real] - y_10[real]) * sin_3 +
                   (y_21[imag] - y_01[imag]) * cos_2 - (y_21[real] - y_01[real]) * sin_2 +
                   (y_10[imag] - y_12[imag]) * cos_1 - (y_10[real] - y_12[real]) * sin_1 +
                   (y_01[imag] - y_21[imag]);

    // Z channel
    Z_temp[real] = (y_20[real] - y_21[real]) * cos_3 + (y_20[imag] - y_21[imag]) * sin_3 +
                   (y_02[real] - y_12[real]) * cos_2 + (y_02[imag] - y_12[imag]) * sin_2 +
                   (y_21[real] - y_20[real]) * cos_1 + (y_21[imag] - y_20[imag]) * sin_1 +
                   (y_12[real] - y_02[real]);

    Z_temp[imag] = (y_20[imag] - y_21[imag]) * cos_3 - (y_20[real] - y_21[real]) * sin_3 +
                   (y_02[imag] - y_12[imag]) * cos_2 - (y_02[real] - y_12[real]) * sin_2 +
                   (y_21[imag] - y_20[imag]) * cos_1 - (y_21[real] - y_20[real]) * sin_1 +
                   (y_12[imag] - y_02[imag]);


    X[real] =  X_temp[real] * cos_offset - X_temp[imag] * sin_offset;
    X[imag] = -X_temp[real] * sin_offset - X_temp[imag] * cos_offset;

    Y_lisa_sim[real] =  Y_temp[real] * cos_offset - Y_temp[imag] * sin_offset;
    Y_lisa_sim[imag] = -Y_temp[real] * sin_offset - Y_temp[imag] * cos_offset;

    Z_lisa_sim[real] =  Z_temp[real] * cos_offset - Z_temp[imag] * sin_offset;
    Z_lisa_sim[imag] = -Z_temp[real] * sin_offset - Z_temp[imag] * cos_offset;

  }

  for (int i = 0; i < 2 * num_time_samples; ++i) {
    A[i] = (2.0 * X[i] - Y_lisa_sim[i] - Z_lisa_sim[i]) / 3.0;
    E[i] = (Z_lisa_sim[i] - Y_lisa_sim[i]) / std::sqrt(3.0);
  }

}


}
