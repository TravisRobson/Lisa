

#include "lisa_binfiles.hpp"

#include <cstring>
#include <fstream>
#include <stdexcept>

#include <iostream> /// \todo Remove.
#include <stdio.h>
namespace lisa {


// Definitions of static members.
constexpr int64_t ParamFileFooter::meta_size;
constexpr char ParamFileFooter::expected_string[15];


ParamMetadataFooter create_params_metadata_footer()
{

  ParamMetadataFooter footer {
    Lisa_VERSION,
    Lisa_PROJECT_NAME,
    0 
  };

  return footer;
}


void write_params_metadata_footer(std::ostream& os, ParamMetadataFooter& footer)
{

  os.write((char*)&footer, sizeof(footer));

}


void write_params_file_footer(std::ostream& os)
{

  ParamFileFooter footer;
  footer.metdata_footer_size = ParamFileFooter::meta_size;
  strcpy(footer.magic_string, ParamFileFooter::expected_string);

  os.write((char*)&footer, sizeof(footer));

}


void write_gb_params(const GB_params& params, std::ostream& os, ParamMetadataFooter& meta_footer)
{

  os.write((char*)&GB_params::id_,         sizeof(GB_params::id_));

  os.write((char*)&params.freq,            sizeof(params.freq));
  os.write((char*)&params.freq_dot,        sizeof(params.freq_dot));
  os.write((char*)&params.theta,           sizeof(params.theta));
  os.write((char*)&params.phi,             sizeof(params.phi));
  os.write((char*)&params.amp,             sizeof(params.amp));
  os.write((char*)&params.iota,            sizeof(params.iota));
  os.write((char*)&params.psi,             sizeof(params.psi));
  os.write((char*)&params.phase,           sizeof(params.phase));
  os.write((char*)&params.freq_double_dot, sizeof(params.freq_double_dot));

  ++meta_footer.gb_count;

}


// void write_gb_params(
//   const std::vector<GB_params>& params, 
//   std::ostream& os, 
//   ParamMetadataFooter& meta_footer)
// {
//   for ( auto& p : params ) {
//     write_gb_params(p, os);
//   }
// }


void read_params_file_footer(ParamFileFooter& footer, std::istream& is)
{

  is.seekg(-sizeof(ParamFileFooter), is.end);

  is.read((char*)&footer.metdata_footer_size, sizeof(footer.metdata_footer_size));
  is.read((char*)&footer.magic_string,        sizeof(footer.magic_string));

}


void read_params_metadata_footer(ParamMetadataFooter& footer, std::istream& is)
{

  is.seekg(-sizeof(ParamMetadataFooter) - sizeof(ParamFileFooter), is.end);

  is.read((char*)&footer.version,  sizeof(footer.version));
  is.read((char*)&footer.app,      sizeof(footer.app));
  is.read((char*)&footer.gb_count, sizeof(footer.gb_count));

}


void read_gb_params(GB_params& params, std::istream& is)
{

  is.read((char*)&params.freq,            sizeof(params.freq));
  is.read((char*)&params.freq_dot,        sizeof(params.freq_dot));
  is.read((char*)&params.theta,           sizeof(params.theta));
  is.read((char*)&params.phi,             sizeof(params.phi));
  is.read((char*)&params.amp,             sizeof(params.amp));
  is.read((char*)&params.iota,            sizeof(params.iota));
  is.read((char*)&params.psi,             sizeof(params.psi));
  is.read((char*)&params.phase,           sizeof(params.phase));
  is.read((char*)&params.freq_double_dot, sizeof(params.freq_double_dot));

}


std::vector<GB_params> read_gb_params(const std::string& filename)
{

  std::ifstream infile(filename, std::ios::in | std::ios::binary);

  // Read the file footer.
  ParamFileFooter file_footer;
  read_params_file_footer(file_footer, infile);

  // Ensure the footer is correct.
  if (file_footer.magic_string != ParamFileFooter::expected_string) {
    throw std::runtime_error(filename + " footer had wrong magic_string. "
      "Had '" + std::string(file_footer.magic_string) + "', expected '" +
      std::string(ParamFileFooter::expected_string) + "'");
  }

  if (file_footer.metdata_footer_size != ParamFileFooter::meta_size) {
    throw std::runtime_error(filename + " had wrong metdata_footer_size");
  }


  // Read the metadata footer.
  ParamMetadataFooter metadata;
  read_params_metadata_footer(metadata, infile);

  /// \todo Do error checking on metadata.


  // Go to the top of the file
  infile.seekg(0);
  int param_type;
  infile.read((char*)&param_type, sizeof(GB_params::id_));

  std::vector<GB_params> params(metadata.gb_count);

  for (auto& p : params) {
    read_gb_params(p, infile);
  }

  infile.close();

  return params;

}


} // end lisa namespace


// ///
// /// C interface for Python consumption.
// ///
// #ifdef __cplusplus
// extern "C" {
// #endif


// void lisa_test() {
//   // std::cout << "I'm in test!\n";
//   // puts("Hello");
//   printf("Hello\n");
// }


// double lisa_add(double a, double b) {
//   return a + b;
// }


// // struct C_GB_params_array {

// // }


// // C_GB_params* lisa_read_gb_params(char* const filename)
// // {

// //   // std::vector<C_GB_params> params(1);

// //   // return params.data();

// //   C_GB_params* pParams = new C_GB_params[1]();

// //   pParams->freq = 1.23456;

// //   return pParams;

// // }


// double freq(C_GB_params* p) { return p->freq; }


// #ifdef __cplusplus
// } // end extern "C"
// #endif

