// Module for X-ray cluster fgas data
// (c) 2013, 2020 Anja von der Linden, Douglas Applegate,
//          Adam Mantz (amantz@stanford.edu)

#ifndef _LENSING_
#define _LENSING_

#include <iostream>
#include <string>
#include <vector>

#include "clusters.hpp"

using std::string;

typedef std::vector<double> Vector;

namespace Clusters {

  class ShearMeasurement {
  public:
    double r_arcmin
      ,    r_Mpc
      ,    g
      ,    sigma_g
      ,    sigma2_g
    ;
    inline void get_radius_Mpc(const double dA) {
      r_Mpc = r_arcmin * dA * (1./60.)*(pi/180.);
    }
  };
  typedef std::vector<ShearMeasurement> ShearProfile;

  class RedshiftPoint {
  public:
    double z, N;
  };
  typedef std::vector<RedshiftPoint> RedshiftHistogram;
  
  class LensingMeasurements { // color-cut method
  public:
    static double c2_over_G; // c^2/G in Msun/Mpc
    double betas, betas2, Sigmacrit; //rho_c_over_Sigma_c;
    ShearProfile shear;
    RedshiftHistogram zhist;
    double chisq(NFWmodel&) const;
    inline bool has_data() {return shear.size()>0;}
    bool load(const string &file);
    void redshift_dependent_calculations(ClusterModel&, const double zcluster, const double dAcluster);

  };

}

std::istream& operator>>(std::istream &is, Clusters::ShearMeasurement&);
std::ostream& operator<<(std::ostream &os, Clusters::ShearMeasurement&);
std::istream& operator>>(std::istream &is, Clusters::RedshiftPoint&);
std::ostream& operator<<(std::ostream &os, Clusters::RedshiftPoint&);
std::ostream& operator<<(std::ostream &os, Clusters::LensingMeasurements&);

#endif

