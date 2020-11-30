

#include <cstdio> // For remove (i.e. to remove files created in tests).
#include <fstream>
#include <stdint.h>

#include "lisa_binfiles.hpp"
#include "lisa_types.hpp"

#include "gtest/gtest.h"


using namespace lisa;


///////////////////////////////////////////////////////////////////////////////
///
/// Design of this code was helped by this SO post:
/// https://stackoverflow.com/questions/416436/what-to-put-in-a-binary-data-files-header
///////////////////////////////////////////////////////////////////////////////
TEST( ParamBinfileTests, WriteReadOneGBparam )
{

  // Initialize GB parameters
  GB_params out_params { };

  out_params.freq            = 1.0e-3;
  out_params.freq_dot        = 1.0e-10;
  out_params.theta           = 1.1;
  out_params.phi             = 2.3;
  out_params.amp             = 1.0e-21;
  out_params.iota            = 0.234;
  out_params.psi             = 0.123456;
  out_params.phase           = 0.5;
  out_params.freq_double_dot = 11.0 / 3.0 * out_params.freq_dot 
                               * out_params.freq_dot / out_params.freq;

  // Write the associated parameter binary file
  std::string filename { "WriteReadFile.par" };                       
  std::ofstream outfile(filename, std::ios::out | std::ios::binary);

  ParamMetadataFooter out_metadata_footer = create_params_metadata_footer();
  write_gb_params(out_params, outfile, out_metadata_footer);
  write_params_metadata_footer(outfile, out_metadata_footer);
  write_params_file_footer(outfile);

  outfile.close();


  // Read the created binary file.
  std::ifstream infile(filename, std::ios::in | std::ios::binary);

  ParamFileFooter file_footer;
  read_params_file_footer(file_footer, infile);

  ParamMetadataFooter metadata;
  read_params_metadata_footer(metadata, infile);

  infile.seekg(0);
  int param_type;
  infile.read((char*)&param_type, sizeof(GB_params::id_));

  GB_params read_params { };
  read_gb_params(read_params, infile);


  infile.close();
  remove(filename.c_str());


  // validate the data read.

  // File footer
  EXPECT_EQ(file_footer.metdata_footer_size, ParamFileFooter::meta_size);
  EXPECT_STREQ("PARAMSFILEFMT", file_footer.magic_string);

  // Metadata footer
  EXPECT_STREQ("0.1.0", metadata.version);
  EXPECT_STREQ("Lisa", metadata.app);
  ASSERT_EQ(metadata.gb_count, 1);

  // GB params
  ASSERT_EQ(param_type, GB_params::id_);

  // These should be EQ, not FLOAT_EQ, since we demand binary equivalence 
  // for reading our custom binary file format.
  EXPECT_EQ(read_params.freq,            out_params.freq);
  EXPECT_EQ(read_params.freq_dot,        out_params.freq_dot);
  EXPECT_EQ(read_params.theta,           out_params.theta);
  EXPECT_EQ(read_params.phi,             out_params.phi);
  EXPECT_EQ(read_params.amp,             out_params.amp);
  EXPECT_EQ(read_params.iota,            out_params.iota);
  EXPECT_EQ(read_params.psi,             out_params.psi);
  EXPECT_EQ(read_params.phase,           out_params.phase);
  EXPECT_EQ(read_params.freq_double_dot, out_params.freq_double_dot);

}


TEST( ParamBinfileTests, NeilToParFormatConverter )
{


}

