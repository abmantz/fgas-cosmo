// Module for X-ray cluster fgas data
// (c) 2013 Adam Mantz, amantz@stanford.edu

#ifndef _CLUSTERUTIL_
#define _CLUSTERUTIL_

namespace Clusters {

  template <class T> inline T const pow2(const T x) {return x*x;}
  template <class T> inline T const pow3(const T x) {return x*x*x;}
  template <class T> inline T const pow5(const T x) {T x2 = x*x; return x2*x2*x;}

  const double pi = 3.14159265358979323846264338328;

}

#endif
