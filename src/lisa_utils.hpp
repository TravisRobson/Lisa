

#ifndef LISA_UTILS_HPP
#define LISA_UTILS_HPP


#include <assert.h>
#include <chrono>
#include <iomanip>
#include <locale>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>

#include "lisa_config.hpp"
#include "lisa_errors.hpp"


static std::string remove_extension( std::string filename ) {

  size_t lastindex = filename.find_last_of("."); 
  return filename.substr(0, lastindex); 

}

#define FILENAME_BASE remove_extension(__FILENAME__)



namespace lisa
{


///////////////////////////////////////////////////////////////////////////////
/// \brief Count the number of lines in a file.
///
/// \param filename name of the file whose lines are to be counted.
///////////////////////////////////////////////////////////////////////////////
int find_num_lines(const std::string& filename); 


///////////////////////////////////////////////////////////////////////////////
/// \brief Format a large number with commans.
///
/// e.g. 10023421 ---> 10,023,421
///
/// https://stackoverflow.com/a/7276879/7228874
///////////////////////////////////////////////////////////////////////////////
class comma_numpunct : public std::numpunct<char>
{
protected:
  virtual char do_thousands_sep() const { return ','; }
  virtual std::string do_grouping() const { return "\03"; }
};


template<class T>
std::string format_with_commans(T value)
{

  std::locale comma_locale(std::locale(), new comma_numpunct());
  // std::cout.imbue(comma_locale);

  std::stringstream ss;
  ss.imbue(comma_locale); //std::locale(""));
  ss << std::fixed << value;

  return ss.str();

}


///////////////////////////////////////////////////////////////////////////////
/// \brief Calculate the median value of a vector of doubles
/// 
/// \param v vector holding onto doubles.
///////////////////////////////////////////////////////////////////////////////
double median(std::vector<double>& v);



///////////////////////////////////////////////////////////////////////////////
/// \brief Simple stopwatch like timer.
///
/// Start and stop the timer. It collects the time differences between start
/// and stop times.
///////////////////////////////////////////////////////////////////////////////
class Timer
{

  std::vector<double> times; ///< vector of times (differences between 
                             ///< start_time and end_time ) in milliseconds.

  /// The start time (like starting a stopwatch)
  std::chrono::time_point<std::chrono::high_resolution_clock> start_time;

  bool started; ///< Has the timer been started?

public:

  Timer() : started(false) {} ///< default constructor

  /// Start the timer.
  void start() { 
    assert(!started); // clock cannot have been started.
    start_time = std::chrono::high_resolution_clock::now(); 
    started = true;
  }

  /// Stop the timer. Append the time difference vector of times.
  void stop();

  /// Return the median value of the times.
  double median_time() {
    return median(times);
  }

  double total_time();

  size_t count() {
    return times.size();
  }


};


///////////////////////////////////////////////////////////////////////////////
///
///////////////////////////////////////////////////////////////////////////////
class Timers
{

#ifdef ENABLE_LISA_PROFILE ///< No need to waste memory if not profiling
  /// Map of Timer specified by a name (i.e. char const*)
  /// Unordered map chosen for faster search times (searching is
  /// is the operation we wish to be fastest so that we don't 
  /// affect the appliction profiling more than we have to).
  std::unordered_map<char const*, std::unique_ptr<Timer>> timers;

  /// Width alloted for std::string of Timer when printed out.
  static constexpr int name_width  { 27 };
  static constexpr int time_width  { 18 };
  static constexpr int count_width { 12 };
  static constexpr int num_fields { 4 }; ///< Number of fields printed when
                                         ///< timing statistics are reported.

  /// Number of decimals of precision to print times.
  static constexpr int precision { 7 }; 
#endif


public:

  Timers() { ///< Default constructor.

#ifdef ENABLE_LISA_PROFILE
    // add_timer("analytic_instrument_noise");
#endif

  }

  /// Add a Timer specified by a std::string.
  /// \param name Name to be given to the Timer
  void add_timer(char const* name) {
#ifdef ENABLE_LISA_PROFILE
    timers.insert(std::make_pair(name, std::make_unique<Timer>()));
#endif
  }

  /// Start a specific timer
  /// \param name Name of the Timer to start.
  void start_timer(char const* name) {
#ifdef ENABLE_LISA_PROFILE
    timers[name]->start();
#endif
  }

  /// Stop a specific timer.
  /// \param name Name of the Timer to stop.
  void stop_timer(char const* name) {
#ifdef ENABLE_LISA_PROFILE
    timers[name]->stop();
#endif
  }

  /// Print out profiling statistics. For example,
  ///
  /// |----------------------------------------------------------------------------------|
  /// |                       name |  median time (ms) |   total time (ms) |  call count |
  /// |----------------------------------------------------------------------------------|
  /// |            allocate memory |         41.190114 |         41.190114 |           1 |
  /// |  analytic_instrument_noise |          0.000220 |        851.840053 |     2985953 |
  /// |           read source file |      17640.363465 |      17640.363465 |           1 |
  /// |----------------------------------------------------------------------------------|
  void report_median_times() const;

};

extern Timers g_timers; ///< Global object declaration used for timing this application.
                        ///< Being explicit about extern.

/// Convenience macros to be used inside functions for profiling.
#define START_FUNC_TIMER() g_timers.start_timer(__func__)
#define STOP_FUNC_TIMER() g_timers.stop_timer(__func__)


///////////////////////////////////////////////////////////////////////////////
/// \brief Print a progress bar to the terminal.
///
/// This function uses a return carriage. This means that it doesn't print out 
/// new lines. However, if you print to screen something before the bar has
/// completed the progress bar will get broken up.
///
/// The following code has been borrowed from the stack overflow post:
/// https://stackoverflow.com/a/36315819/7228874
///
/// \param percentage Value between 0 and 1, fraction complete.
///////////////////////////////////////////////////////////////////////////////
#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

void print_progress(double percentage);


struct ThreeVec {

  double v[3];

  double& operator[](int c) { 
    // assert(c < 3 && c >= 0); 
    return v[c]; 
  }
  const double& operator[](int c) const { 
    // assert(c < 3 && c >= 0); 
    return v[c]; 
  }

  // e.g. ThreeVec result = -Vec
  ThreeVec operator-() const {
    ThreeVec result {
      -v[0],
      -v[1],
      -v[2]
    };

    return result;
  }

  // e.g. ThreeVec result = Vec - second.
  ThreeVec operator-(const ThreeVec& second) const {
    ThreeVec result {
      v[0] - second.v[0],
      v[1] - second.v[1],
      v[2] - second.v[2]
    };

    return result;
  }

  ThreeVec operator/(double value) const {
    assert(value != 0.0);
    ThreeVec result {
      v[0] / value,
      v[1] / value,
      v[2] / value
    };

    return result;
  }


  ThreeVec& operator/=(double rhs) {
    assert(rhs != 0.0);
    v[0] /= rhs;
    v[1] /= rhs;
    v[2] /= rhs;

    return *this;
  }

  ThreeVec& operator*=(double rhs) {
    v[0] *= rhs;
    v[1] *= rhs;
    v[2] *= rhs;

    return *this;
  }

};

inline double dot(const ThreeVec& v1, const ThreeVec& v2) {

  double result { v1[0] * v2[0] +  v1[1] * v2[1] + v1[2] * v2[2] };
  return result;

}


/// \todo is a flattened t[9] preferred?
struct ThreeTensor { // By tensor we mean rank 2.
  double t[9];

  double& operator()(int c1, int c2) { 
    // assert(c1 < 3 && c1 >= 0); 
    // assert(c2 < 3 && c2 >= 0); 
    return t[c1 * 3 + c2];
  }
  const double& operator()(int c1, int c2) const { 
    assert(c1 < 3 && c1 >= 0); 
    // assert(c2 < 3 && c2 >= 0); 
    return t[c1 * 3 + c2];
  }

  // double t[3][3];

  // double& operator()(int c1, int c2) { 
  //   assert(c1 < 3 && c1 >= 0); 
  //   assert(c2 < 3 && c2 >= 0); 
  //   return t[c1][c2];
  // }
  // const double& operator()(int c1, int c2) const { 
  //   assert(c1 < 3 && c1 >= 0); 
  //   assert(c2 < 3 && c2 >= 0); 
  //   return t[c1][c2];
  // }

};

inline double dot(const ThreeTensor& t, const ThreeVec& v1, const ThreeVec& v2) {

  double result { 
    v1[0] * (t(0, 0) * v2[0] + t(0, 1) * v2[1] + t(0, 2) * v2[2]) +
    v1[1] * (t(1, 0) * v2[0] + t(1, 1) * v2[1] + t(1, 2) * v2[2]) +
    v1[2] * (t(2, 0) * v2[0] + t(2, 1) * v2[1] + t(2, 2) * v2[2])
  };

  return result;

}


}


#endif
