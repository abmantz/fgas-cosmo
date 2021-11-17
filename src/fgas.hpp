// Module for X-ray cluster fgas data
// (c) 2013, 2020 Adam Mantz, amantz@stanford.edu

#ifndef _FGAS_
#define _FGAS_

#include <iostream>
#include <string>
#include <vector>

#include "clusters.hpp"
#include "lensing.hpp"
#include "xray.hpp"

using std::string;

typedef std::vector<double> Vector;

namespace Clusters {

  class Cluster {
  public:
    bool do_lensing_calibration, use_fgas;
    string name, data_path;
    double z // redshift
      ,    last_lnP
      ,    ln_lensing_mass
    ;
    struct {
      NuisanceParameters *nuisance;
      double rho_c, lncalibration;
      NFWmodel mass;
    } for_integrals;
    XrayMeasurements xray;
    LensingMeasurements lensing;
    double lnP(ClusterModel&, NuisanceParameters&, const double dA, const double dL, const double rhocr);
    double lnP_lensing(ClusterModel&, NuisanceParameters&, const double dA, const double rhocr);
    double lnP_lensing_integrand(const double lnMwl);
    bool load_parameters(const double lnMwl);
  };
  typedef std::vector<Cluster*> ClusterVector;

  class Dataset {
  public:
    bool do_lensing_calibration, use_fgas;
    string data_path;
    Vector trial_dA, trial_dL, trial_rhocr;
    ClusterVector clusters;
    bool init(const string &file);
    double lnP(ClusterModel&, NuisanceParameters&);
    void load_trial(double *dAs, double *dLs, double *rhocrs);
    inline double redshift(const int i) const {return clusters[i]->z;}
    void simulate(ClusterModel&, NuisanceParameters&, const bool iscatter=false, const bool mscatter=false);
    inline size_t size() const {return clusters.size();}
  };

  double fgas_prediction(Cluster&, ClusterModel&, NuisanceParameters&, const double dA, const double dL, const double rhocr);

}

std::istream& operator>>(std::istream &is, Clusters::Cluster&);
std::ostream& operator<<(std::ostream &os, Clusters::Cluster&);
std::istream& operator>>(std::istream &is, Clusters::Dataset&);
std::ostream& operator<<(std::ostream &os, Clusters::Dataset&);

#endif
