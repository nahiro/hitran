#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <gsl/gsl_spline.h>

/// sigma: Gaussian std
/// gamma: Lorentzian half width
///    x0: X center
///  area: Integrated area
///    np: X array size
///    xp: X
///    yp: Y
int FFTVoigt(const double sigma,const double gamma,
             const double x0,const double area,
             const int np,const double* const xp,
             double* const yp);
int FTApproxVoigt(const double sigma,const double gamma,
                  const double x0,const double area,
                  const int np,const double* const xp,
                  double* const yp);

#ifdef __cplusplus
}
#endif
