

#ifndef LISA_BINFILES_HPP
#define LISA_BINFILES_HPP


#include <iostream>
#include <stdint.h>
#include <string>
#include <vector>

#include "lisa_config.hpp"
#include "lisa_types.hpp"


namespace lisa {



///////////////////////////////////////////////////////////////////////////////
/// \brief Parameter-Binfile footer metadata.
///
/// By parameter we mean MCMC parameters which would influence the model of 
/// LISA's TDI data.
///////////////////////////////////////////////////////////////////////////////
struct ParamMetadataFooter {
  char version[40]; ///< Version of this footer
  char app[40];     ///< Name of application "Lisa"
  int32_t gb_count; ///< Number of GBs in source file.
};


///////////////////////////////////////////////////////////////////////////////
/// \brief Footer at the end of Parameter-Binfiles
///////////////////////////////////////////////////////////////////////////////
struct ParamFileFooter {
  int64_t metdata_footer_size; ///< i.e. sizeof(ParamMetadataFooter)
  char    magic_string[15];    ///< "PARAMFILEFMT" Used to verify correctness of file.

  /// Get size of metadata
  static constexpr int64_t meta_size { sizeof(ParamMetadataFooter) };
  static constexpr char expected_string[15] { "PARAMSFILEFMT" };
};


ParamMetadataFooter create_params_metadata_footer();


///////////////////////////////////////////////////////////////////////////////
/// \brief Write the Parameter binfile metadata footer to output stream as 
///        bytes.
///
/// \param os Output stream, likely a filestream. Binary file you're 
///           printing parameters to.
///////////////////////////////////////////////////////////////////////////////
void write_params_metadata_footer(std::ostream& os, ParamMetadataFooter& footer);


///////////////////////////////////////////////////////////////////////////////
/// \brief Write the Parameter binfile footer to output stream as bytes.
///
/// \param os Output stream, likely a filestream. Binary file you're 
///           printing parameters to.
///////////////////////////////////////////////////////////////////////////////
void write_params_file_footer(std::ostream& os);


///////////////////////////////////////////////////////////////////////////////
/// \breif Write GB parameters to output stream as bytes.
/// 
/// \param params GB parameters to be written to binary file.
/// \param os Output stream, likely a filestream, where parameter are written.
///////////////////////////////////////////////////////////////////////////////
void write_gb_params(const GB_params& params, std::ostream& os, ParamMetadataFooter& footer);


// void write_gb_params(
//   const std::vector<GB_params>& params, 
//   std::ostream& os,
//   ParamMetadataFooter& meta_footer);


///////////////////////////////////////////////////////////////////////////////
///
///////////////////////////////////////////////////////////////////////////////
void read_params_file_footer(ParamFileFooter& footer, std::istream& is);


///////////////////////////////////////////////////////////////////////////////
///
///////////////////////////////////////////////////////////////////////////////
void read_params_metadata_footer(ParamMetadataFooter& footer, std::istream& is);


///////////////////////////////////////////////////////////////////////////////
///
///////////////////////////////////////////////////////////////////////////////
void read_gb_params(GB_params& params, std::istream& is);


std::vector<GB_params> read_gb_params(const std::string& filename);


} // end lisa namespace


///
/// C interface for Python consumption.
///
// #ifdef __cplusplus
// extern "C" {
// #endif


// void lisa_test();


// double lisa_add(double a, double b);


// struct C_GB_params {
//   double freq;            ///< Frequency (HZ)
//   double freq_dot;        ///< Frequency time derivative, df/dt, ("f-dot")
//   double theta;           ///< Polar angle on sky (solar system barycenter coord)
//   double phi;             ///< Azimuthal angle on sky
//   double amp;             ///< Amplitude
//   double iota;            ///< Binary system inclination angle
//   double psi;             ///< GW polarization angle
//   double phase;           ///< Initial GW phase angle
//   double freq_double_dot; ///< Second time derivative
// };

// C_GB_params* lisa_read_gb_params(char* const filename);

// double freq(C_GB_params* p); 


// #ifdef __cplusplus
// } // end extern "C"
// #endif


#endif
