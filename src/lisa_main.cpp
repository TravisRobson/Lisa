

#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "lisa_config.hpp"
#include "lisa_constants.hpp"
#include "lisa_detector.hpp"
#include "lisa_errors.hpp"
#include "lisa_math.hpp"
#include "lisa_types.hpp"
#include "lisa_utils.hpp"

#ifdef LISA_USE_OPENMP
# include <omp.h>
#endif

#ifdef ENABLE_LISA_LOGGING
#include <plog/Log.h> 
#endif

//#define USE_BINARY


void process_args(int argc, char* argv[]);
void run();


/// Number of time samples for a year's worth of data at 15 second cadence.
///
/// Note, int (4 Bytes most likely on your machine) caps out at
/// 2,147,483,647. Largest number of samples we expect for LISA at full
/// cadence of 1s is 
/// 1s * 365.25 days / yr * 24 hours / day * 3600 sec/hour * 10yrs
/// = 315,576,000. So there isn't a need for long int here.
///
/// Seconds in a year = 2,103,840 dividing by 15 seconds = 2,103,840
/// and 2,097,152 = 2^21, i.e. closest power of 2.
constexpr int samples_1_year {2'097'152};

///
constexpr double observation_period { samples_1_year * lisa::constants::delta_t };

///
static std::vector<lisa::XAE> XAE_channels;

///
static std::vector<lisa::XAE> XAE_source_data;

///
static std::vector<lisa::GB_params> gb_params; 

///
static std::vector<lisa::GB_params> bright_gb_params;

///
constexpr const char* root_path(LISA_ROOT_PATH);

///
#ifndef USE_BINARY
static std::string source_filename { "data/sources/Galaxy_196810p.dat" };
#else
static std::string source_filename { "data/sources/Galaxy_196810p.bin" };
#endif

///
static std::string log_filename { std::string(FILENAME_BASE) + ".log" };

///
static std::string bright_filename { "bright_gbs.dat" };

///
static std::string tdi_filename { "galaxy_tdi.dat" };

///
constexpr double amp_upper_limit { 2.0 };


///
constexpr double freq_upper_limit { 20.0e-3 }; // 20 mHz


/// Memory pool for single GB signals (complex)
std::vector<double> X_pool;
std::vector<double> A_pool;
std::vector<double> E_pool;


std::vector<double> X_channel;
std::vector<double> A_channel;
std::vector<double> E_channel;


int main(int argc, char* argv[])
{

#ifdef ENABLE_LISA_LOGGING
  // https://github.com/SergiusTheBest/plog/tree/1.1.5
  plog::init(plog::debug, log_filename.c_str()); 
#endif

  try {

    process_args(argc, argv);

    auto start_time { std::chrono::system_clock::now() };

    run();

    auto end_time { std::chrono::system_clock::now() };

    std::chrono::duration<double> elapsed_time{end_time - start_time};
    std::cout << "Elapse time: " << elapsed_time.count() << " seconds." << std::endl;

  }
  catch(std::exception& e) {
    std::cerr << "Exception caught at " << __FILENAME__ << ":" << __LINE__   
              << "\n\t message: " << e.what() << std::endl;

    return 1;
  }
  catch(...) {
    std::cerr << "Unknown exception caught at " << __FILENAME__ << ":" 
              << __LINE__ << std::endl;    

    return 2;
  }

  return 0;

}

///////////////////////////////////////////////////////////////////////////////
/// \brief  Process command line arguments.
///
/// Options:
///   --source-file  This is the file which contains GB parameters.
///   --help  Print a help message describing these options to the user.
///
/// \todo Error handling for bad inputs.
///////////////////////////////////////////////////////////////////////////////
void process_args(int argc, char* argv[])
{

#ifdef ENABLE_LISA_LOGGING
  for (int i = 0; i < argc; ++i) {
    PLOG(plog::info) << "argv[" << i << "]: " << argv[i]; 
  }
#endif

  for (int i = 1; i < argc; ++i) {

    std::string arg { argv[i] };

    if (arg == "--source-file" || arg == "-s" ) {
      source_filename = std::string(argv[i + 1]);
      ++i;
    }
    else if (arg == "--help" || arg == "-h") {

      std::cout << std::endl;
      std::cout << "Usage:\n";
      std::cout << "\t$ ./bin/lisa_driver <args> [options]\n";

      std::cout << "Example Usage:\n";
      std::cout << "\t$ ./bin/lisa_driver --source-file example.txt\n";

      std::cout << "Options:\n";
      std::cout << "\t--source-file, -s    GB parameter source file.\n";
      std::cout << "\t--help,        -h    Print this help message.\n"; 
      std::cout << std::endl;


      throw std::invalid_argument(std::string(__FILENAME__) + " won't run when --help is invoked.");

    }
    else{
      throw std::invalid_argument("User passed invalid command line argument: " + arg);
    }

  }

  if (!file_accessible(source_filename)) {
    std::string at{std::string(__FILENAME__) + ":" + std::string(TO_STRING(__LINE__))};
    throw std::runtime_error("At " + at + " file (" + source_filename + ") not accessible.\n");
  }

  std::cout << "LISA version: " << lisa::get_version() << "\n";
  std::cout << "Command line configuration:\n";
  std::cout << "\troot path .......................... " << root_path << "\n";
  std::cout << "\tsource file ........................ " << source_filename << "\n";
  std::cout << std::endl;

#ifdef ENABLE_LISA_LOGGING
  PLOG(plog::info) << "root path: " << root_path;
  PLOG(plog::info) << "source file: " << source_filename; 
#endif

}


///////////////////////////////////////////////////////////////////////////////
/// \brief
///////////////////////////////////////////////////////////////////////////////
void run() 
{

  //
  // Add timer tags.
  //
  lisa::g_timers.add_timer("allocate memory");
  lisa::g_timers.add_timer("read source file");
  lisa::g_timers.add_timer("simulate GBs");
  lisa::g_timers.add_timer("calc GB signal");
  lisa::g_timers.add_timer("calc summed signals");
  lisa::g_timers.add_timer("write signal file");


  //
  // Allocate memory needed by the application.
  //
  lisa::g_timers.start_timer("allocate memory");
  XAE_channels.resize(samples_1_year);
  XAE_source_data.reserve(1<<11); // 2048 (1<<11, 1<<10 is 1024) is a decent start size.
  lisa::g_timers.stop_timer("allocate memory");

  //
  // Read in source parameters from file into memory.
  // \todo Error handling for bad data in file.
  //
#ifndef USE_BINARY
  std::ifstream source_file { root_path + source_filename };
#else // Below code is used if using binary format
  std::ifstream source_file;
  source_file.open(root_path + source_filename, std::ios::binary | std::ios::in);
#endif
  std::string line { };

  lisa::g_timers.start_timer("read source file");

#ifndef USE_BINARY
  while (std::getline(source_file, line)) {

    lisa::GB_params params { };
    lisa::fill_gb_params(params, line);
    gb_params.push_back(params);

  }
 #else 
  while (  source_file.peek() != EOF ) {

    lisa::GB_params params { };
    lisa::fill_gb_params_binary(params, source_file);
    gb_params.push_back(params);

  } 
#endif


  lisa::g_timers.stop_timer("read source file");

  source_file.close();

  std::cout << "Read " << lisa::format_with_commans(gb_params.size()) 
            << " GB parameter sets.\n"; 

#ifdef ENABLE_LISA_LOGGING
  PLOG(plog::info) << "source_file (" << source_filename << ") read, "
                   << lisa::format_with_commans(gb_params.size()) << " GB params.";
#endif

  //
  // 
  //
  int AMCVn_count { 0 }; 
  int bright_gbs { 0 };
  std::cout << "Simulating galactic binaries..." << std::endl;
  lisa::g_timers.start_timer("simulate GBs");

  // Open up bright GB param file 
  std::ofstream bright_file;
  bright_file.open(root_path + bright_filename);

  // Reserve memory for GB signals
  X_pool.reserve(2 * 8192);
  A_pool.reserve(2 * 8192);
  E_pool.reserve(2 * 8192);


  // To hold the sum total of all galactic binary signals
  X_channel.resize(samples_1_year, 0.0);
  A_channel.resize(samples_1_year, 0.0);
  E_channel.resize(samples_1_year, 0.0);


  // // Strains
  // lisa::Strains strains_buffer(2 * 8192);

  for ( size_t param_idx = 0; param_idx < gb_params.size(); ++param_idx) {

    const lisa::GB_params& params = gb_params[ param_idx ];
    if (params.freq_dot < 0.0) {
      ++AMCVn_count;
    }

    std::pair<size_t, double> out { lisa::gb_bandwidth(params, observation_period) };

    size_t num_samples { out.first };
    double cutoff_amp { out.second };

    if (cutoff_amp > amp_upper_limit) {
      ++bright_gbs;
    }

    // Only simulate binaries below an upper limit (20 mHz).
    if (params.freq < freq_upper_limit) {
    
      X_pool.resize(2 * num_samples);
      A_pool.resize(2 * num_samples);
      E_pool.resize(2 * num_samples);

      //strains_buffer.resize(2 * num_samples);

      lisa::g_timers.start_timer("calc GB signal");
      lisa::calc_gb_signal(params, X_pool, A_pool, E_pool, observation_period);
      lisa::g_timers.stop_timer("calc GB signal");

      /////////////////////////////////////////////////////////////////////////
      /// \todo move this function into Lisa library.
      lisa::g_timers.start_timer("calc summed signals");
      for (size_t i = 0; i < num_samples; ++i) {

        size_t carrier_freq_bin { static_cast<size_t>(params.freq * observation_period) };

        size_t k { carrier_freq_bin + i - num_samples / 2 };

        size_t real { 2 * k };
        size_t imag { 2 * k + 1 };

        X_channel[real] += X_pool[2 * i];
        X_channel[imag] += X_pool[2 * i + 1];

        A_channel[real] += A_pool[2 * i];
        A_channel[imag] += A_pool[2 * i + 1];

        E_channel[real] += E_pool[2 * i];
        E_channel[imag] += E_pool[2 * i + 1];

      }
      lisa::g_timers.stop_timer("calc summed signals");
      /////////////////////////////////////////////////////////////////////////

    }

    //
    // Print prgress bar every 10,000 binaries simulated.
    //
    if ( param_idx % 10000 == 0 ) {
      double percent_done { (double)param_idx / gb_params.size() };
      lisa::print_progress(percent_done);
    }

  }
  lisa::print_progress(1.0);
  bright_file.close();
  lisa::g_timers.stop_timer("simulate GBs");


  /////////////////////////////////////////////////////////////////////////
  /// \todo move this function into Lisa library.
  std::ofstream tdi_file;
  tdi_file.open(root_path + tdi_filename, std::ios::out);
  tdi_file << std::scientific;

  int max_freq_index { static_cast<int>(std::ceil(20e-3 * observation_period)) };
  double sqrt_obs_period { std::sqrt(observation_period) };

  lisa::g_timers.start_timer("write signal file");
  for (int i = 1; i < max_freq_index; ++i) {

    double freq { i / observation_period };

    tdi_file << freq << " ";

    // Add in a realization of instrument noise
    lisa::NoiseXAE noise { lisa::analytic_instrument_noise_legacy(freq) };

    tdi_file << sqrt_obs_period * X_channel[2 * i + 0] + 0.5 * std::sqrt(noise.XYZ) * lisa::gaussian_variate() << " ";
    tdi_file << sqrt_obs_period * X_channel[2 * i + 1] + 0.5 * std::sqrt(noise.XYZ) * lisa::gaussian_variate() << " ";

    tdi_file << sqrt_obs_period * A_channel[2 * i + 0] + 0.5 * std::sqrt(noise.AE) * lisa::gaussian_variate() << " ";
    tdi_file << sqrt_obs_period * A_channel[2 * i + 1] + 0.5 * std::sqrt(noise.AE) * lisa::gaussian_variate() << " ";

    tdi_file << sqrt_obs_period * E_channel[2 * i + 0] + 0.5 * std::sqrt(noise.AE) * lisa::gaussian_variate() << " ";
    tdi_file << sqrt_obs_period * E_channel[2 * i + 1] + 0.5 * std::sqrt(noise.AE) * lisa::gaussian_variate() << " ";

    tdi_file << "\n";


  }
  tdi_file.close();
  lisa::g_timers.stop_timer("write signal file");
  /////////////////////////////////////////////////////////////////////////


  std::cout << std::endl;
  std::cout << std::fixed << "Read " << lisa::format_with_commans(AMCVn_count) << " AMCVn binaries.\n";
  std::cout << std::fixed << "Read " << lisa::format_with_commans(bright_gbs) << " bright GBs.\n";

#ifdef ENABLE_LISA_LOGGING
  PLOG(plog::info) << "GBs simulated, "
                   << lisa::format_with_commans(AMCVn_count) << " AMCVn binaries.";
  PLOG(plog::info) << lisa::format_with_commans(bright_gbs) << " bright GBs.\n";
#endif

  std::cout << std::endl;
  lisa::g_timers.report_median_times();

}


