// Module for X-ray cluster fgas data
// (c) 2013, 2020 Anja von der Linden, Douglas Applegate,
//          Adam Mantz (amantz@stanford.edu)

#include "lensing.hpp"
#include "util.hpp"

#include <sstream>

#include "ConfigFile.h"

namespace Clusters {

  double LensingMeasurements::c2_over_G = 2.090e19; // c^2/G in Msun/Mpc

  bool LensingMeasurements::load(const string &file) {
    std::ifstream fin;
    fin.open(file.c_str());
    if ( !fin.good() ) return false;
    fin.close();

    ConfigFile config(file);
    string st;
    
    // read in shear profile
    if (!config.readInto(st, "gprof")) return false;
    else {
      std::istringstream iss(st);
      while (true) {
	shear.resize(shear.size() + 1);
	iss >> shear.back();
	if (!iss) break;
      }
      shear.pop_back();
    }
    if (shear.size() == 0) return false;

    // read in z histo
    if (!config.readInto(st, "zhisto")) return false;
    else {
      std::istringstream iss(st);
      while (true) {
	zhist.resize(zhist.size() + 1);
	iss >> zhist.back();
	if (!iss) break;
      }
      zhist.pop_back();
    }
    if (zhist.size() == 0) return false;

    //std::cout << "Read lensing data from " << file << std::endl;
    return true;
  }

  void LensingMeasurements::redshift_dependent_calculations(ClusterModel &mod, const double zcluster, const double dAcluster) {
    const double infinity = 1.0e12;

    for(int i=0; i<shear.size(); ++i) shear[i].get_radius_Mpc(dAcluster);

    double betainf = mod.AngularDiameterDistance(zcluster, infinity) / mod.AngularDiameterDistance(0.0, infinity);
    double sum=0.0, sum2=0.0, n=0.0, beta;
    for(int i=0; i<zhist.size(); ++i) {
      n += zhist[i].N;
      if (zhist[i].z > zcluster) {
	beta = mod.AngularDiameterDistance(zcluster, zhist[i].z) / mod.AngularDiameterDistance(0.0, zhist[i].z);
	sum += (zhist[i].N * beta);
	sum2 += (zhist[i].N * pow2(beta));
      }
    }
    betas = sum / (n * betainf);
    betas2 = sum2 / (n * pow2(betainf));

    Sigmacrit = c2_over_G / (4.0*pi * dAcluster * betainf);
  }

  double LensingMeasurements::chisq(NFWmodel &mass) const {
    double gamma_inf, kappa_inf, gmodel, chi2=0.0;
    for(int i=0; i<shear.size(); ++i) {
      gamma_inf = mass.shear(shear[i].r_Mpc, Sigmacrit);
      kappa_inf = mass.kappa(shear[i].r_Mpc, Sigmacrit);
      gmodel = (betas*gamma_inf)/(1.0 - betas2/betas*kappa_inf); // Seitz&Schneider 1997, A2.4
      chi2 += pow2(shear[i].g - gmodel) / shear[i].sigma2_g;
    }
    return chi2;
  }


}

std::istream& operator>>(std::istream &is, Clusters::ShearMeasurement &x) {
  is >> x.r_arcmin >> x.g >> x.sigma_g;
  x.sigma2_g = Clusters::pow2(x.sigma_g);
  return is;
}

std::ostream& operator<<(std::ostream &os, Clusters::ShearMeasurement &x) {
  const char s = ' ';
  return os << x.r_arcmin <<s<< x.g <<s<< x.sigma_g;
}

std::istream& operator>>(std::istream &is, Clusters::RedshiftPoint &x) {
  return is >> x.z >> x.N;;
}

std::ostream& operator<<(std::ostream &os, Clusters::RedshiftPoint &x) {
  const char s = ' ';
  return os << x.z <<s<< x.N;
}

std::ostream& operator<<(std::ostream &os, Clusters::LensingMeasurements &x) {
  const char s = ' ';
  return os << "LensShear:"<<x.shear.size() <<s<< "LensZ:"<<x.zhist.size();
}
