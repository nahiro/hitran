#include <complex.h>
#include <fftw3.h>
#include <gsl/gsl_spline.h>

int FFTVoigt(const double sigma,const double gamma,
             const double x0,const double area,
             const int np,const double* const xp,double* const yp)
{
  static const double threshold=1.0e-50; // Y threshold
  static const double step=0.02; // X sampling interval in units of FWHM
  static const double width=1.0e4; // X width in units of FWHM
  const double xp_min = xp[0]; // xp must be aligned in ascending order
  const double xp_max = xp[np-1];
  const double xdif_min = (xp_min>=x0?xp_min-x0:(xp_max<=x0?x0-xp_max:0.0));
  const double xdif_max = (xp_min>=x0?xp_max-x0:(xp_max<=x0?x0-xp_min:MAX(xp_max-x0,x0-xp_min)));
  // Voigt profile width approximation by J.J.Olivero and R.L.Longbothum
  // JQSRT, Vol.17, pp.233-236 (1977)
  static const double w1 = 1.0692;
  static const double w2 = 0.86639;
  static const double w3 = 5.5451774444795623; // =log(2.0)*8.0
  const double fwhm = w1*gamma+sqrt(w2*gamma*gamma+w3*sigma*sigma);
  const double xstp = MAX(MIN(fwhm*step,xdif_max*step),MAX(xdif_max*2.0,fwhm*width)*1.0e-5);
  const int nbin = (int)(pow(2.0,ceil(log(MAX(xdif_max*2.0,fwhm*width)/xstp)/log(2.0)))+0.1);
  const int n_half = nbin/2;
  const int n_zero = nbin-n_half-1;
  const double twid = 2.0*M_PI/xstp;
  const double tstp = twid/nbin;
  const double norm = area/(xstp*nbin);
  double t;
  double fact;
  const double sign = (area<0.0?-1.0:1.0);
  static const int nmgn = 16;
  const int nx = n_half+1;
  const int nx_min = MAX((int)(floor(xdif_min/xstp)+0.1)-nmgn,0);
  const int nx_max = MIN((int)(ceil(xdif_max/xstp)+0.1)+nmgn,nx);
  const int nx_neg = MIN(nmgn,nx_max-nx_min-1);
  const int nx_dif = nx_max-nx_min+(nx_min==0?nx_neg:0);
  const int index = (xp_min>=x0?-1:(xp_max<=x0?np-1:gsl_interp_bsearch(xp,x0,0,np)));
  double yi;
  double *xn;
  if((xn=(double*)malloc(nx_dif*sizeof(double))) == NULL)
  {
    fprintf(stderr,"Error, failed in allocating memory.\n");
    return -1;
  }
  double *yn;
  if((yn=(double*)malloc(nx_dif*sizeof(double))) == NULL)
  {
    fprintf(stderr,"Error, failed in allocating memory.\n");
    free(xn);
    return -1;
  }
  double complex *p;
  if((p=(double complex*)malloc(nbin*sizeof(double complex))) == NULL)
  {
    fprintf(stderr,"Error, failed in allocating memory.\n");
    free(xn);
    free(yn);
    return -1;
  }
  double complex *q;
  if((q=(double complex*)malloc(nbin*sizeof(double complex))) == NULL)
  {
    fprintf(stderr,"Error, failed in allocating memory.\n");
    free(xn);
    free(yn);
    free(p);
    return -1;
  }
  int i,j;

  // fill p array
  p[0] = 1.0;
  for(i=1; i<=n_zero; i++)
  {
    t = tstp*i;
    fact = exp(-gamma*t-0.5*sigma*sigma*t*t);
    p[nbin-i] = p[i] = fact;
    if(fact < threshold) // No need to take fabs (fact > 0)
    {
      while(i < n_zero)
      {
        i++;
        p[nbin-i] = p[i] = 0.0;
      }
      break;
    }
  }
  if(n_half > n_zero)
  {
    if(fact < threshold)
    {
      p[n_half] = 0.0;
    }
    else
    {
      t = tstp*n_half;
      fact = exp(-gamma*t-0.5*sigma*sigma*t*t);
      p[n_half] = fact;
    }
  }

  // fill q array
  fftw_complex *inp = (fftw_complex*)p;
  fftw_complex *out = (fftw_complex*)q;
  fftw_plan fwd_plan = fftw_plan_dft_1d(nbin,inp,out,FFTW_FORWARD,FFTW_ESTIMATE);
  fftw_execute(fwd_plan);
  fftw_destroy_plan(fwd_plan);

  // fill x and y array
  if(nx_min == 0)
  {
    for(i=nx_min,j=nx_neg; i<nx_max; i++,j++)
    {
      xn[j] = xstp*i;
      yn[j] = creal(q[i])*norm;
    }
    for(i=nx_neg+1,j=nx_neg-1; j>=0; i++,j--)
    {
      xn[j] = -xn[i];
      yn[j] = yn[i];
    }
  }
  else
  {
    for(i=nx_min,j=0; i<nx_max; i++,j++)
    {
      xn[j] = xstp*i;
      yn[j] = creal(q[i])*norm;
    }
  }

  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline,nx_dif);
  gsl_spline_init(spline,xn,yn,nx_dif);
  for(i=index+1; i<np; i++)
  {
    yi = gsl_spline_eval(spline,xp[i]-x0,acc);
    if(yi*sign < threshold)
    {
      break;
    }
    yp[i] += yi;
  }
  for(i=index; i>=0; i--)
  {
    yi = gsl_spline_eval(spline,x0-xp[i],acc);
    if(yi*sign < threshold)
    {
      break;
    }
    yp[i] += yi;
  }
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  free(xn);
  free(yn);
  free(p);
  free(q);

  return 0;
}

int FTApproxVoigt(const double sigma,const double gamma,
                  const double x0,const double area,
                  const int np,const double* const xp,
                  double* const yp)
{
  static const int nmax = 20;
  static const double threshold = 1.0e-50;
  static const int nr = 150;
  static const double r = 1.1;
  static const double width = 5.0e3;
  static const int t1 = 12;
  static const int t2 = 144; // t1*t1
  int n;
  double n2p2[nmax];
  double n4p4[nmax];
  double pn[nmax];
  double en[nmax];
  double sn[nmax];
  static const double p1 = M_PI;
  static const double pp = 6.283185307179586; // p1+p1
  static const double p2 = 9.869604401089358; // p1*p1
  static const double p_2 = 1.1283791670955126; // 2/sqrt(pi)
  static const double s_2 = 0.70710678118654746; // 1/sqrt(2)
  static const double t_2 = 0.006944444444444444; // 1/t2
  const double s_s = s_2/sigma;
  const double y = s_s*gamma;
  const double y2 = y*y;
  const double ty = t1*y;
  const double t2y = t2*y;
  const double e_ty = exp(-ty);
  const double te_ty = 0.5*t1*e_ty;
  const double norm = area*p_2*s_s/t1;
  const double sign = (area<0.0?-1.0:1.0);
  double v;
  double x;
  double x2;
  double tx;
  double x2y2;
  double t2x2y2;
  double t2y2x2;
  double t4x4y4;
  double ptx;
  double n2p2t2x2y2;
  double stx;
  double ctx;
  double xstx;
  double yctx;
  double tyctx;

  const double xp_min = xp[0]; // xp must be aligned in ascending order
  const double xp_max = xp[np-1];
  const double xdif_min = (xp_min>=x0?xp_min-x0:(xp_max<=x0?x0-xp_max:0.0));
  const double xdif_max = (xp_min>=x0?xp_max-x0:(xp_max<=x0?x0-xp_min:MAX(xp_max-x0,x0-xp_min)));
  // Voigt profile width approximation by J.J.Olivero and R.L.Longbothum
  // JQSRT, Vol.17, pp.233-236 (1977)
  static const double w1 = (0.5*1.0692);
  static const double w2 = (0.25*0.86639);
  static const double w3 = 1.3862943611198906; // log(2)*2
  const double hwhm = w1*gamma+sqrt(w2*gamma*gamma+w3*sigma*sigma);
  const double x1 = MAX(hwhm*width,xdif_max)*(r-1.0)/(pow(r,(double)nr)-1.0);
  const double rx_1 = (r-1.0)/x1;
  const double lr_1 = 1.0/log(r);
  static const int nmgn = 16;
  const int nx = nr+1;
  const int nx_min = MAX((int)(floor(log(xdif_min*rx_1+1.0)*lr_1)+0.1)-nmgn,0);
  const int nx_max = MIN((int)(ceil(log(xdif_max*rx_1+1.0)*lr_1)+0.1)+nmgn,nx);
  const int nx_neg = MIN(nmgn,nx_max-nx_min-1);
  const int nx_dif = nx_max-nx_min+(nx_min==0?nx_neg:0);
  double xn[nx_dif];
  double yn[nx_dif];
  double dx;
  double sx;
  double yi;
  const int index = (xp_min>=x0?-1:(xp_max<=x0?np-1:gsl_interp_bsearch(xp,x0,0,np)));
  int i,j;

  for(n=0; n<nmax; n++)
  {
    n2p2[n] = n*n*p2;
    n4p4[n] = n2p2[n]*n2p2[n];
    pn[n] = p1*n;
    en[n] = exp(-n2p2[n]*t_2);
    sn[n] = (n%2==1?-1.0:1.0)*te_ty;
  }
  dx = x1;
  sx = 0.0;
  for(i=1; i<=nx_min; i++)
  {
    sx += dx;
    dx *= r;
  }
  if(nx_min == 0)
  {
    for(i=nx_min,j=nx_neg; i<nx_max; i++,j++)
    {
      xn[j] = sx;
      sx += dx;
      dx *= r;
      x = s_s*xn[j];
      x2 = x*x;
      tx = t1*x;
      x2y2 = x2+y2;
      t2x2y2 = t2*x2y2;
      t4x4y4 = t2x2y2*t2x2y2;
      t2y2x2 = 2.0*t2*(y2-x2);
      ptx = pp*tx;
      stx = sin(tx);
      ctx = cos(tx);
      xstx = x*stx;
      yctx = y*ctx;
      tyctx = ty*ctx;
      v = (-y+(yctx-xstx)*e_ty)/(2.0*x2y2);
      for(n=0; n<nmax; n++)
      {
        n2p2t2x2y2 = n2p2[n]+t2x2y2;
        v += en[n]*(t2y*n2p2t2x2y2/(n4p4[n]+n2p2[n]*t2y2x2+t4x4y4)+
             sn[n]*(((tx-pn[n])*stx-tyctx)/(n2p2t2x2y2-n*ptx)+
                    ((tx+pn[n])*stx-tyctx)/(n2p2t2x2y2+n*ptx)));
      }
      yn[j] = norm*v;
      if(yn[j]*sign < threshold)
      {
        for(i=i+1,j=j+1; i<nx_max; i++,j++)
        {
          xn[j] = sx;
          yn[j] = 0.0;
          sx += dx;
          dx *= r;
        }
        break;
      }
    }
    for(i=nx_neg+1,j=nx_neg-1; j>=0; i++,j--)
    {
      xn[j] = -xn[i];
      yn[j] = yn[i];
    }
  }
  else
  {
    for(i=nx_min,j=0; i<nx_max; i++,j++)
    {
      xn[j] = sx;
      sx += dx;
      dx *= r;
      x = s_s*xn[j];
      x2 = x*x;
      tx = t1*x;
      x2y2 = x2+y2;
      t2x2y2 = t2*x2y2;
      t4x4y4 = t2x2y2*t2x2y2;
      t2y2x2 = 2.0*t2*(y2-x2);
      ptx = pp*tx;
      stx = sin(tx);
      ctx = cos(tx);
      xstx = x*stx;
      yctx = y*ctx;
      tyctx = ty*ctx;
      v = (-y+(yctx-xstx)*e_ty)/(2.0*x2y2);
      for(n=0; n<nmax; n++)
      {
        n2p2t2x2y2 = n2p2[n]+t2x2y2;
        v += en[n]*(t2y*n2p2t2x2y2/(n4p4[n]+n2p2[n]*t2y2x2+t4x4y4)+
             sn[n]*(((tx-pn[n])*stx-tyctx)/(n2p2t2x2y2-n*ptx)+
                    ((tx+pn[n])*stx-tyctx)/(n2p2t2x2y2+n*ptx)));
      }
      yn[j] = norm*v;
      if(yn[j]*sign < threshold)
      {
        for(i=i+1,j=j+1; i<nx_max; i++,j++)
        {
          xn[j] = sx;
          yn[j] = 0.0;
          sx += dx;
          dx *= r;
        }
        break;
      }
    }
  }

  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline,nx_dif);
  gsl_spline_init(spline,xn,yn,nx_dif);
  for(i=index+1; i<np; i++)
  {
    yi = gsl_spline_eval(spline,xp[i]-x0,acc);
    if(yi*sign < threshold)
    {
      break;
    }
    yp[i] += yi;
  }
  for(i=index; i>=0; i--)
  {
    yi = gsl_spline_eval(spline,x0-xp[i],acc);
    if(yi*sign < threshold)
    {
      break;
    }
    yp[i] += yi;
  }
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);

  return 0;
}
