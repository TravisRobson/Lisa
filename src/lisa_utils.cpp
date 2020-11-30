

#include "lisa_utils.hpp"

#include <algorithm> 
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>


namespace lisa {


Timers g_timers; ///< Definition of global Timers object.


int find_num_lines(const std::string& filename)
{

  int result { 0 };

  std::string line{};
  std::ifstream file{ filename };

  while(std::getline(file, line)) {
    ++result;
  }

  return result;

}


// https://stackoverflow.com/questions/1719070/what-is-the-right-approach-when-using-stl-container-for-median-calculation/1719155#1719155
double median(std::vector<double>& v)
{

  if (v.empty()) {
    return 0.0; /// \todo should I throw an exception here?
  }

  size_t half_size { v.size() / 2 };

  std::nth_element(v.begin(), v.begin() + half_size, v.end());

  if (v.size() % 2 != 0) { // i.e. odd number of elements in vector
    return v[half_size];
  }
  else { // even number of elements in vector
    auto other_value { std::max_element(v.begin(), v.begin() + half_size) };
    return (v[half_size] + *other_value) / 2.0;
  }

}


void print_progress(double percentage)
{

  int val{ (int)(percentage * 100.0) };
  int lpad{ (int)(percentage * PBWIDTH) };
  int rpad{ PBWIDTH - lpad };

  printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
  fflush(stdout);

}


void Timer::stop() 
{

  assert(started); // Clock must have been started.

  std::chrono::nanoseconds elapsed_time  = 
    std::chrono::duration_cast<std::chrono::nanoseconds>(
      std::chrono::high_resolution_clock::now() - start_time);

  times.push_back(elapsed_time.count());
  started = false;

}


/// Function put here to avoid <numeric> include in every file which includes
/// lisa_utils.hpp
double Timer::total_time() 
{
  return std::accumulate(times.begin(), times.end(), 0.0);
}



void Timers::report_median_times() const
{
#ifdef ENABLE_LISA_PROFILE
  std::cout << "LISA Profiling Statistics\n";

  const std::string separator { " |" }; // Used to delimit values in table

  // Total width of a line of information
  const size_t total_width { 
    name_width + 2 * time_width + count_width + num_fields * separator.size() };

  // Construct a horizontal line of "------------------------"
  const std::string line = separator + std::string( total_width - 1, '-' ) + '|';

  std::cout << line << '\n' << separator
            << std::setw(name_width) << "name" << separator
            << std::setw(time_width) << "median time (ms)" << separator
            << std::setw(time_width) << "total time (ms)" << separator
            << std::setw(count_width) << "call count" << separator 
            << "\n" << line << "\n";

  /// \todo It would be nice if I could print this in alphabetical order
  for (auto& timer : timers) {

    std::string name { timer.first }; 
    double median_time { timer.second->median_time() / 1000000.0 };
    double total_time { timer.second->total_time() / 1000000.0 };
    size_t count = timer.second->count();

    std::cout << separator << std::setw(name_width) << name << separator
              << std::setw(time_width) << median_time << separator
              << std::setw(time_width) << total_time << separator
              << std::setw(count_width) << format_with_commans(count) << separator << "\n";

  }
  std::cout << line << '\n';

#endif
}


}

