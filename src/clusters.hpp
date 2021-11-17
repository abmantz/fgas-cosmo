// Module for X-ray cluster fgas data
// (c) 2013 Adam Mantz, amantz@stanford.edu

#ifndef _CLUSTERMODEL_
#define _CLUSTERMODEL_

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "util.hpp"

using std::pow;
using std::string;

typedef std::vector<double> Vector;
typedef double(*func2)(const double, const double);

namespace Clusters {

  class Profile {
    // simple power-law implementation; NB inputs and index are squared
  public:
    double index;
    Profile() : index(0.0) {}
    inline double operator()(const double x2_model, const double x2_measure) {return pow(x2_measure/x2_model, 0.5*index);}
  };

  class NFWmodel {
  public:
    double nfwa, nfwpot; // Nulsen NFWMASS parametrization!
    inline static double delta_c(const double c, const double Delta=200.0) {
      return Delta/3.0* c*c*c / ( log(1.0+c) - c/(1.0+c) );
    }
    static double f_inverse(const double); // inverse of x^3*(log(1+1/x)-1/(1+x)) from Hu & Kravtsov 2003 ApJ 584:702 Appendix C
    inline static double f_fcn(const double x) {
      // nb this is x^3*unn_mass(1/x)
      return pow3(x) * ( log(1.0+1.0/x) - 1.0/(1.0+x) );
    }
    inline double kappa(const double r, const double Sigmacrit) const {
      // predicted convergence as a function of radius
      return nfwpot/(2.0*pi * nfwa * Sigmacrit) * unn_kappa(r/nfwa);
    }
    inline double mass(const double r) const {
      return nfwa * nfwpot * unn_mass(r/nfwa);
    }
    inline void set_from_rs_c_rhocr(const double rs, const double c, const double rhoc, const double Delta=200.0) {
      nfwa = rs;
      nfwpot = 4.0*pi * pow2(rs) * rhoc * delta_c(c, Delta);
    }
    inline double shear(const double r, const double Sigmacrit) const{
      // predicted shear as a function of radius
      return nfwpot/(4.0*pi * nfwa * Sigmacrit) * unn_shear(r/nfwa);
    }
    static double unn_kappa(const double x); // un-normalized convergence as a function of scaled radius
    inline static double unn_mass(const double x) {
      // un-normalized mass as a function of scaled radius
      return log(1.0 + x) - x / (1.0 + x);
    }
    static double unn_shear(const double x); // un-normalized shear as a function of scaled radius
  };
  typedef std::vector<NFWmodel> mcNFWmodel;

  class ClusterModel {
  public:
    static double mHe_mp; // Ratio of alpha mass to proton mass - has a default value, but can be overwritten.
    func2 AngularDiameterDistance; // pointer to some external function
    double U_baryon_0   // cluster parameters
      ,    U_baryon_1
      ,    U_coldgas_0
      ,    U_coldgas_1
      ,    fgas_scatter // fgas will be short for f_hotgas
      ,    massScaling  // power-law dependence of fgas on Mtot
      // not a free parameter:
      ,    massPivot    // pivot for fgas(M)
      ;
    double fb_cosmic    // copied or derived from external cosmology class
      ,    ne_nH        // ratio of ICM number densities
      ,    ntot_ne      // ditto
      ,    mu           // mean particle mass of ICM (in proton masses)
    ;
    Profile fb_prof, mtot_prof;
    inline double f_baryon(const double z, const double x2_model, const double x2_measure) {
      // NB no check for negative values
      return U_baryon_0 * (1.0 + U_baryon_1*z) * fb_prof(x2_model, x2_measure) * fb_cosmic;
    }
    inline double f_coldgas(const double z) {
      // NB no check for negative values
      return U_coldgas_0 * (1.0 + U_coldgas_1*z) * fb_cosmic;
    }
    bool init(func2 &f, const string& filename="");
    bool load_parameters(const double Ub0, const double Ub1, const double Uc0, const double Uc1, const double scat, const double rslope, const double Mslope, const double Mscaling, const double fbcos, const double YHe=-1, const double nenh=-1, const double ntne=-1, const double myoo=-1);
    void set_from_YHe(const double YHe);
  };

  class NuisanceParameters {
  public:
    bool do_lensing_calibration, use_fgas;
    double nonthermal, calibration, calibration_evol, calibration_scatter, lensing_systematics, lensing_systematics_evol;
    Vector lensing_concentrations, lensing_concentrations500;
    NuisanceParameters();
    inline double calibration_z(const double z) const {
      return calibration * (1 + z*calibration_evol)*(lensing_systematics + (z-0.4)*lensing_systematics_evol);
      // z=0.4 is the appropriate pivot for simulations in Applegate 2016
    }
    inline void load_parameters(const double nth, const double cal, const double calev, const double calscat, const double lenssys, const double lenssysev) {
      nonthermal = nth;
      calibration = cal;
      calibration_evol = calev;
      calibration_scatter = calscat;
      lensing_systematics = lenssys;
      lensing_systematics_evol = lenssysev;
    }
  };

}

std::istream& operator>>(std::istream &is, Clusters::NFWmodel&);
std::ostream& operator<<(std::ostream &os, Clusters::NFWmodel&);
std::ostream& operator<<(std::ostream &os, Clusters::NuisanceParameters&);

#endif
