// Module for X-ray cluster fgas data
// (c) 2013 Adam Mantz, amantz@stanford.edu

#ifndef _FWRAPPER_
#define _FWRAPPER_

#include "wrapper.hpp"

namespace Clusters {
  Wrapper *theWrapper;
  double cangluardistance2(const double, const double);
}

extern "C" {

  typedef double fortran_dl;
  typedef float fortran_real;
  typedef long fortran_integer;
  typedef bool fortran_logical1;
  typedef char fortran_character;

  fortran_integer fclwrapperdatasetsize_(fortran_integer *i);
  fortran_dl fclwrappergetredshift_(fortran_integer *i, fortran_integer *j);
  fortran_integer fclwrapperinit_(fortran_character *a_char, fortran_integer a_len);
  fortran_dl fclwrapperlnp_(fortran_dl *datapars);
  fortran_integer fclwrapperloaddata_(fortran_character *a_char, fortran_integer a_len);
  fortran_integer fclwrapperloadparameters_(fortran_dl *modelpars);
  fortran_integer fclwrappernumdatapars_();
  fortran_integer fclwrappernumparameters_();
  fortran_integer fclwrappernumclusters_();
  fortran_integer fclwrappernumconst_();
  void fclwrappersetconst_(fortran_dl *physics);
  void fclwrappersetdataoptions_(fortran_integer *i, fortran_integer *j);
  void fclwrappersetmasspivot_(fortran_dl *M);
  void fclwrappersimulate_(fortran_dl *datapars, fortran_integer *iscatter, fortran_integer *mscatter);

  fortran_dl fangulardistance2_(fortran_dl*, fortran_dl*); // provided by cosmomc
  
}

#endif
