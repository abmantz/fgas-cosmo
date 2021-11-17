// Module for X-ray cluster fgas data
// (c) 2013 Adam Mantz, amantz@stanford.edu

#include "fwrapper.hpp"

namespace Clusters {

  double cangulardistance2(const double z1, const double z2) {
    double zz1=z1, zz2=z2;
    return fangulardistance2_(&zz1, &zz2);
  }

}

extern "C" {

  fortran_integer fclwrapperdatasetsize_(fortran_integer *i) {
    return ClWrapperDatasetSize(Clusters::theWrapper, *i);
  }

  fortran_dl fclwrappergetredshift_(fortran_integer *i, fortran_integer *j) {
    return ClWrapperGetRedshift(Clusters::theWrapper, *i, *j);
  }

  fortran_integer fclwrapperinit_(fortran_character *a_char, fortran_integer a_len) {
    Clusters::theWrapper = newClWrapper();
    return ClWrapperInit(Clusters::theWrapper, &Clusters::cangulardistance2, a_char);
  }

  fortran_dl fclwrapperlnp_(fortran_dl *datapars) {
    return ClWrapperLnP(Clusters::theWrapper, datapars);
  }

  fortran_integer fclwrapperloaddata_(fortran_character *a_char, fortran_integer a_len) {
    std::string s(a_char, a_len);
    return ClWrapperLoadData(Clusters::theWrapper, s.c_str());
  }

  fortran_integer fclwrapperloadparameters_(fortran_dl *modelpars) {
    return ClWrapperLoadParameters(Clusters::theWrapper, modelpars);
  }

  fortran_integer fclwrappernumdatapars_() {
    return ClWrapperNumDatapars(Clusters::theWrapper);
  }

  fortran_integer fclwrappernumparameters_() {
    return ClWrapperNumParameters(Clusters::theWrapper);
  }

  fortran_integer fclwrappernumclusters_() {
    return ClWrapperNumClusters(Clusters::theWrapper);
  }

  fortran_integer fclwrappernumconst_() {
    return ClWrapperNumConst(Clusters::theWrapper);
  }

  void fclwrappersetconst_(fortran_dl *physics) {
    ClWrapperSetConst(Clusters::theWrapper, physics);
  }

  void fclwrappersetdataoptions_(fortran_integer *i, fortran_integer *j) {
    ClWrapperSetDataOptions(Clusters::theWrapper, *i, *j);
  }

  void fclwrappersetmasspivot_(fortran_dl *M) {
    ClWrapperSetMassPivot(Clusters::theWrapper, *M);
  }

  void fclwrappersimulate_(fortran_dl *datapars, fortran_integer *iscatter, fortran_integer *mscatter) {
    ClWrapperSimulate(Clusters::theWrapper, datapars, *iscatter, *mscatter);
  }

}
