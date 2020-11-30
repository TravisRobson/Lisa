

#ifndef LISA_TYPES_H
#define LISA_TYPES_H


#include <stdint.h>


namespace lisa {


///////////////////////////////////////////////////////////////////////////////
/// \brief X, A, E channel values at a given frequency or time sample.
///
/// These channels are frequently used together in calculations.
/// \todo I might want other combinations of these data channels for certain
///       calculations.
/// \todo What source (paper) should I quote here which describes these data 
///       channels?
///////////////////////////////////////////////////////////////////////////////
struct XAE {
  double X; ///< X Michelson channel.
  double A; ///< A channel.
  double E; ///< E channel.
};


///////////////////////////////////////////////////////////////////////////////
/// \brief Parameters specifying LISA GW signal of slowly evolving GB.
///////////////////////////////////////////////////////////////////////////////
struct GB_params {
  double freq;            ///< Frequency (HZ)
  double freq_dot;        ///< Frequency time derivative, df/dt, ("f-dot")
  double theta;           ///< Polar angle on sky (solar system barycenter coord)
  double phi;             ///< Azimuthal angle on sky
  double amp;             ///< Amplitude
  double iota;            ///< Binary system inclination angle
  double psi;             ///< GW polarization angle
  double phase;           ///< Initial GW phase angle
  double freq_double_dot; ///< Second time derivative

  static constexpr int id_{ 1 };  ///< ID for this structure. 
  static constexpr int num_{ 9 }; ///< 9 params in this struct.
};


///////////////////////////////////////////////////////////////////////////////
/// \brief  X and AE noise PSD values for a given frequency or time sample.
///
/// Note: X, Y, Z channels have the same noise PSD. 
///
/// T channel is ignored in this struct because frequently that channel doesn't
/// need to be considered in a calculation.
///////////////////////////////////////////////////////////////////////////////
struct NoiseXAE {
  double XYZ; ///< XYZ Michelson channel noise PSD
  double AE;  ///< AE channels noise PSD
};


}


#endif
