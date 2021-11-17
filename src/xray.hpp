// Module for X-ray cluster fgas data
// (c) 2013 Adam Mantz, amantz@stanford.edu

#ifndef _XRAY_
#define _XRAY_

#include <iostream>
#include <string>
#include <vector>

#include "clusters.hpp"

using std::string;

typedef std::vector<double> Vector;

namespace Clusters {

  // converts to imaging-spectroscopy units from physical quantities
  class XrayConversions {
  public:
    static double G  // gravitational constant in keV Mpc/Msun^2
      ,           mp // proton mass in Msun
      ,      Mpc_cgs // 1 Mpc in cm
      ;              // These have default values, but can be changed for consistency with another code.
    double distance, density, mtot, mgas;
    void evaluate(const double z, const double rpp, const double ne_nH, const double ntot_ne, const double mu, const double dA, const double dL);
  };

  class XrayMeasurements {
  public:
    double rpp    // radians per fit radial unit
      ,    x1	  // inner radius of measurement (units of rDelta,measured)
      ,    x2	  // outer radius of measurement
      ,    fgas	  // Mgas(r1<r<r2) / Mtot(r1<r<r2)
      ,    ferr   // measurement error (natural log)
      ,    M      // Mtot(0<r<r2)
      ,    Merr   // measurement error (natural log)
      ,    errcor // correlation coefficient relating measurement errors
      ,    dA	  // angular diameter distance assumed in measuring rDelta (Mpc)
      , G_mp_mu_rhocr // other physical factor used in measuring rDelta (Mpc^-2)
      ;
    XrayConversions conv;
    Vector lens_cal_radius; // units consistent with rpp
    mcNFWmodel mass;
    inline bool has_mass_data() {return mass.size()>0;}
    bool load(const string&);
  };

}

std::istream& operator>>(std::istream &is, Clusters::XrayMeasurements&);
std::ostream& operator<<(std::ostream &os, Clusters::XrayMeasurements&);

#endif
