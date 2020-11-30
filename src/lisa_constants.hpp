

#include <cmath>


namespace lisa {


namespace constants {


  constexpr double mod_freq { 3.168753575e-8 };   ///< LISA modulation frequency (Hz)
  constexpr double arm_length { 2.5e09 };         ///< Distance between spacecraft (meters)
  constexpr double delta_t { 15.0 };              ///< Detector cadence (s)
  constexpr double kappa { 0.0 };                 ///<
  constexpr double lambda { 0.0 };                ///<
  constexpr double freq_star { 0.0190853806 };    ///< Transfer frequency (Hz)
  constexpr double eccentricity { 0.0048241852 }; ///< LISA spacecraft eccentricity. 


  /// The following noise values were read from:
  /// https://iopscience.iop.org/article/10.1088/1742-6596/840/1/012024/pdf
  constexpr double position_noise { 8.9e-23 }; ///< (m^2 Hz^-1)
  constexpr double accel_noise { 9.0e-30 };    ///< (m^2 Hz^-1)


  constexpr double c_light { 299'792'458.0 }; ///< Speed of light vacuum (m/s)
  constexpr double year { 3.15581498e7 };     ///< sidereal year (s)
                                              ///< https://en.wikipedia.org/wiki/Sidereal_year
                                              ///< where sidereal year is 365.256 363 004 days.
  constexpr double au { 1.49597870660e11 };   ///< Astronomical Unit (m)
  constexpr double pi { M_PI };


}


}
