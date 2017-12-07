/***********************************************************/
/* HITRAN_COMMON     ... Calculate absorption cross        */
/*                       section using HITRAN database.    */
/* Author: N.Manago                                        */
/* $Revision: 187 $                                         */
/* $Date: 2014-12-06 12:40:51 +0900 (Sat, 06 Dec 2014) $   */
/***********************************************************/
#define		PI			M_PI			// PI
#define		PI2			6.283185307179586	// 2*PI
#define		PI_2			M_PI_2			// PI/2
#define		D_TO_R			1.7453292519943295e-02	// PI/180    deg  -> rad
#define		R_TO_D			57.295779513082323	// 180/PI    rad  -> deg
#define		GAUSS_FWHM		2.3548200450309493	// 2.0*sqrt(2.0*log(2.0))
#define		LOREN_FWHM		2.0			// 2.0
#define		GAUSS_NORM		0.39894228040143270	// 1.0/sqrt(2.0*pi)
#define		LOREN_NORM		M_1_PI			// 1.0/pi
#define		VOIGT_NORM1		M_SQRT1_2		// 1.0/sqrt(2.0)
#define		VOIGT_NORM2		0.56418958354775628	// 1.0/sqrt(pi)
#define		C_LIGHT			2.99792458e8		// Speed of light in m/s
#define		H_PLANK			6.62606896e-34		// Plank constant in J*s
#define		N_A			6.02214179e23		// Abogadoro constant in /mol
#define		K_B			1.3806504e-23		// Boltzmann constant in J/K
#define		C_RAD2			1.4387751601679206	// Second radiation constant (H_PLANK*C_LIGHT/K_B) in cm*K
#define		NMOL96			36
#define		NMOL04			38
#define		NMOL08			41
#define		NMOL12			51
#define		NMOL16			51
#define		MAXISO04		20
#define		MAXISO08		20
#define		MAXISO12		20
#define		MAXISO16		20
#define		NSPECI			85
#define		NRANGE			3
#define		NPOL			4
#define		HITRAN96_LINE_LENGTH	100
#define		HITRAN04_LINE_LENGTH	160
#define		TREF			296.0			// Reference temperature in K
#define		DNAM			"qt.dat"
#define		EPSILON			0.00001
#define		MAXCHAR			256
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <bits/nan.h>
#include <errno.h>
#include <gsl/gsl_spline.h>
#include "faddeeva.h"

static inline double MIN(double a,double b)
{
  return b < a ? b : a;
}

static inline double MAX(double a,double b)
{
  return b > a ? b : a;
}

struct hitran96
{
  int nmol;
  int niso;
  double freq;
  double sint;
  double mtrn;
  double wair;
  double wslf;
  double elow;
  double ctmp;
  double psft;
  int gqup;
  int gqlo;
  char lqup[MAXCHAR];
  char lqlo[MAXCHAR];
  int err1;
  int err2;
  int err3;
  int ref1;
  int ref2;
  int ref3;
};

struct hitran04
{
  int nmol;
  int niso;
  double freq;
  double sint;
  double aein;
  double wair;
  double wslf;
  double elow;
  double ctmp;
  double psft;
  char gqup[MAXCHAR];
  char gqlo[MAXCHAR];
  char lqup[MAXCHAR];
  char lqlo[MAXCHAR];
  int err1;
  int err2;
  int err3;
  int err4;
  int err5;
  int err6;
  int ref1;
  int ref2;
  int ref3;
  int ref4;
  int ref5;
  int ref6;
  char flag[MAXCHAR];
  double swup;
  double swlo;
};

double		tref			= TREF;
int		isonm[NMOL96];
double		qref[NSPECI];
double		qcoef[NRANGE][NSPECI][NPOL];
char		dnam[MAXCHAR]		= DNAM;

double tips96(int im,int ii,double t);
double tips04(int im,int ii,double t);
double tips08(int im,int ii,double t);
double tips12(int im,int ii,double t);
double tips16(int im,int ii,double t);
double (*tips)(int im,int ii,double t);
double intensity(double t,const struct hitran96 *h);
double width(double t,double p,double ps,const struct hitran96 *h);
double origin(double p,const struct hitran96 *h);
double lorentz(double g,double f0,double f);
double gauss(double g,double f0,double f);
int add_lorentz(double gamma,double xc,double norm,int n,const double *x,double *y);
int add_gauss(double sigma,double xc,double norm,int n,const double *x,double *y);
int add_voigt(double sigma,double gamma,double x0,double area,int np,const double *xp,double *yp);
int get_coeff(void);
int get_hitran96(char *s,struct hitran96 *h);
int get_hitran04(char *s,struct hitran04 *h);
int cnv_hitran04(char *s,struct hitran96 *h);
int get_chr(char *str,int start,int size,char *value);
int get_int(char *str,int start,int size,int *value);
int get_hex(char *str,int start,int size,int *value);
int get_dbl(char *str,int start,int size,double *value);
extern void bd_tips_2004_();
extern void bd_tips_2008_();
extern void bd_tips_2012_();
extern void bd_tips_2016_();

int get_coeff(void)
{
  int i,j,n;
  int err;
  int n1,n2;
  char line[MAXCHAR];
  char str1[MAXCHAR];
  char str2[MAXCHAR];
  char str3[MAXCHAR];
  char str4[MAXCHAR];
  char *endp;
  FILE *fp;

  if((fp=fopen(dnam,"r")) == NULL)
  {
    fprintf(stderr,"Error, cannot open %s\n",dnam);
    return -1;
  }
  err = 0;
  do
  {
    n = 1;
    if(fgets(line,MAXCHAR,fp) == NULL)
    {
      fprintf(stderr,"Error in reading line %d.\n",n);
      err = 1;
      break;
    }
    if(sscanf(line,"%s%s",str1,str2)!=2 || strcmp(str1,"NMOL")!=0 || strcmp(str2,"NSPECI")!=0)
    {
      fprintf(stderr,"Error in line %d >>> %s\n",n,line);
      err = 1;
      break;
    }
    n++;
    if(fgets(line,MAXCHAR,fp) == NULL)
    {
      fprintf(stderr,"Error in reading line %d.\n",n);
      err = 1;
      break;
    }
    if(sscanf(line,"%d%d",&n1,&n2)!=2 || n1!=NMOL96 || n2!=NSPECI)
    {
      fprintf(stderr,"Error in line %d >>> %s\n",n,line);
      err = 1;
      break;
    }
    n++;
    if(fgets(line,MAXCHAR,fp) == NULL)
    {
      fprintf(stderr,"Error in reading line %d.\n",n);
      err = 1;
      break;
    }
    if(sscanf(line,"%s",str1)!=1 || strcmp(str1,"ISONM")!=0)
    {
      fprintf(stderr,"Error in line %d >>> %s\n",n,line);
      err = 1;
      break;
    }
    for(i=0; i<NMOL96; i++)
    {
      n++;
      if(fgets(line,MAXCHAR,fp) == NULL)
      {
        fprintf(stderr,"Error in reading line %d.\n",n);
        err = 1;
        break;
      }
      if(sscanf(line,"%s",str1) != 1)
      {
        fprintf(stderr,"Error in line %d >>> %s\n",n,line);
        err = 1;
        break;
      }
      errno = 0;
      isonm[i] = strtol(str1,&endp,10);
      if(errno==ERANGE || *endp!='\0')
      {
        fprintf(stderr,"Convert error in line %d >>> %s\n",n,line);
        err = 1;
        break;
      }
    }
    if(err) break;
    n++;
    if(fgets(line,MAXCHAR,fp) == NULL)
    {
      fprintf(stderr,"Error in reading line %d.\n",n);
      err = 1;
      break;
    }
    if(sscanf(line,"%s",str1)!=1 || strcmp(str1,"Q296")!=0)
    {
      fprintf(stderr,"Error in line %d >>> %s\n",n,line);
      err = 1;
      break;
    }
    for(i=0; i<NSPECI; i++)
    {
      n++;
      if(fgets(line,MAXCHAR,fp) == NULL)
      {
        fprintf(stderr,"Error in reading line %d.\n",n);
        err = 1;
        break;
      }
      if(sscanf(line,"%s",str1) != 1)
      {
        fprintf(stderr,"Error in line %d >>> %s\n",n,line);
        err = 1;
        break;
      }
      errno = 0;
      qref[i] = strtod(str1,&endp);
      if(errno==ERANGE || *endp!='\0')
      {
        fprintf(stderr,"Convert error in line %d >>> %s\n",n,line);
        err = 1;
        break;
      }
    }
    if(err) break;
    for(j=0; j<NRANGE; j++)
    {
      n++;
      if(fgets(line,MAXCHAR,fp) == NULL)
      {
        fprintf(stderr,"Error in reading line %d.\n",n);
        err = 1;
        break;
      }
      if(sscanf(line,"%s%s%s%s",str1,str2,str3,str4)!=4 || strncmp(str1,"QCOEF",5)!=0
                                                        || strncmp(str2,"QCOEF",5)!=0
                                                        || strncmp(str3,"QCOEF",5)!=0
                                                        || strncmp(str4,"QCOEF",5)!=0)
      {
        fprintf(stderr,"Error in line %d >>> %s\n",n,line);
        err = 1;
        break;
      }
      for(i=0; i<NSPECI; i++)
      {
        n++;
        if(fgets(line,MAXCHAR,fp) == NULL)
        {
          fprintf(stderr,"Error in reading line %d.\n",n);
          err = 1;
          break;
        }
        if(sscanf(line,"%s%s%s%s",str1,str2,str3,str4) != 4)
        {
          fprintf(stderr,"Error in line %d >>> %s\n",n,line);
          err = 1;
          break;
        }
        errno = 0;
        qcoef[j][i][0] = strtod(str1,&endp);
        if(errno==ERANGE || *endp!='\0')
        {
          fprintf(stderr,"Convert error in line %d >>> %s\n",n,line);
          err = 1;
          break;
        }
        errno = 0;
        qcoef[j][i][1] = strtod(str2,&endp);
        if(errno==ERANGE || *endp!='\0')
        {
          fprintf(stderr,"Convert error in line %d >>> %s\n",n,line);
          err = 1;
          break;
        }
        errno = 0;
        qcoef[j][i][2] = strtod(str3,&endp);
        if(errno==ERANGE || *endp!='\0')
        {
          fprintf(stderr,"Convert error in line %d >>> %s\n",n,line);
          err = 1;
          break;
        }
        errno = 0;
        qcoef[j][i][3] = strtod(str4,&endp);
        if(errno==ERANGE || *endp!='\0')
        {
          fprintf(stderr,"Convert error in line %d >>> %s\n",n,line);
          err = 1;
          break;
        }
      }
      if(err) break;
    }
    if(err) break;
    n++;
    if(fgets(line,MAXCHAR,fp) != NULL)
    {
      fprintf(stderr,"Error, too many lines %d.\n",n);
      err = 1;
      break;
    }
  }
  while(0);
  fclose(fp);
  if(err)
  {
    return -1;
  }

  return 0;
}

double tips96(int im,int ii,double t)
{
  int i,n;
  int mol,iso;
  int irange;
  double qt;

  mol = im-1;
  iso = ii-1;
  if(ii > isonm[mol])
  {
    fprintf(stderr,"No such isotope(%d) for molecule %d.\n",ii,im);
    return NAN;
  }

  n = 0;
  for(i=0; i<mol; i++)
  {
    n += isonm[i];
  }
  n += iso;

  if(fabs(t-tref) < EPSILON)
  {
    qt = qref[n];
  }
  else
  {
    if(t<70.0 || t>3005.0)
    {
      fprintf(stderr,"Error, temperature out of range >>> %13.4e\n",t);
      return NAN;
    } else
    if(t <= 500.0)
    {
      irange = 0;
    } else
    if(t <= 1500.0)
    {
      irange = 1;
    }
    else
    {
      irange = 2;
    }
    qt  = qcoef[irange][n][0];
    qt += qcoef[irange][n][1]*t;
    qt += qcoef[irange][n][2]*t*t;
    qt += qcoef[irange][n][3]*t*t*t;
  }

  return qt;
}

double tips04(int im,int ii,double t)
{
  double qt,gi;

  if(im<1 || im>NMOL04 || ii<1 || ii>MAXISO04)
  {
    fprintf(stderr,"Invalid input >>> (%d,%d)\n",im,ii);
    return NAN;
  } else
  if(t<70.0 || t>3000.0)
  {
    fprintf(stderr,"Error, temperature out of range >>> %13.4e\n",t);
    return NAN;
  }
  bd_tips_2004_(&im,&t,&ii,&gi,&qt); // Note: no data for CH3OH (39)

  return qt;
}

double tips08(int im,int ii,double t)
{
  double qt,gi;

  if(im<1 || im>NMOL08 || ii<1 || ii>MAXISO08)
  {
    fprintf(stderr,"Invalid input >>> (%d,%d)\n",im,ii);
    return NAN;
  } else
  if(t<70.0 || t>3000.0)
  {
    fprintf(stderr,"Error, temperature out of range >>> %13.4e\n",t);
    return NAN;
  }
  bd_tips_2008_(&im,&t,&ii,&gi,&qt); // Note: no data for CO2_838, CO2_837, CH4_312,
                                     // C2H6_1231, CH3OH (39), CH3Br (40), CH3CN (41), CF4 (42)

  return qt;
}

double tips12(int im,int ii,double t)
{
  double qt,gi;

  if(im<1 || im>NMOL12 || ii<1 || ii>MAXISO12)
  {
    if(im==2 && ii==0)
    {
      ii = 10;
    }
    else
    {
      fprintf(stderr,"Invalid input >>> (%d,%d)\n",im,ii);
      return NAN;
    }
  }
  if(t<70.0 || t>3000.0)
  {
    fprintf(stderr,"Error, temperature out of range >>> %13.4e\n",t);
    return NAN;
  }
  bd_tips_2012_(&im,&t,&ii,&gi,&qt); // Note: CO2_827->CO2_728
                                     // no data for HNO3_156, HF_29, HCl_25, HCl_27,
                                     // HBr_29, HBr_21, HI_27, N2_45, COF2_369, SO3 (47)
                                     // H2 (45)->47
  if(fabs(qt) < 1.0e-50)
  {
    fprintf(stderr,"Warning, qt=%13.6e is modified to be 1.0 at t=%13.6e\n",qt,t);
    return 1.0;
  }

  return qt;
}

double tips16(int im,int ii,double t)
{
  double qt,gi;

  if(im<1 || im>NMOL16 || ii<1 || ii>MAXISO16)
  {
    if(im==2 && ii==0)
    {
      ii = 10;
    }
    else
    {
      fprintf(stderr,"Invalid input >>> (%d,%d)\n",im,ii);
      return NAN;
    }
  }
  if(t<70.0 || t>3000.0)
  {
    fprintf(stderr,"Error, temperature out of range >>> %13.4e\n",t);
    return NAN;
  }
  bd_tips_2016_(&im,&t,&ii,&gi,&qt);
  if(fabs(qt) < 1.0e-50)
  {
    fprintf(stderr,"Warning, qt=%13.6e is modified to be 1.0 at t=%13.6e\n",qt,t);
    return 1.0;
  }

  return qt;
}

double intensity(double t,const struct hitran96 *h)
{
  double s;

  s  = h->sint;
  s *= tips(h->nmol,h->niso,tref)/tips(h->nmol,h->niso,t);
  s *= exp(-C_RAD2*h->elow/t)/exp(-C_RAD2*h->elow/tref);
  s *= (1.0-exp(-C_RAD2*h->freq/t))/(1.0-exp(-C_RAD2*h->freq/tref));

  return s;
}

double width(double t,double p,double ps,const struct hitran96 *h)
{
  return (h->wair*(p-ps)+h->wslf*ps)*pow(tref/t,h->ctmp);
}

double origin(double p,const struct hitran96 *h)
{
  return h->freq+h->psft*p;
}

double lorentz(double g,double f0,double f)
{
  return LOREN_NORM*g/(g*g+(f-f0)*(f-f0));
}

double gauss(double s,double f0,double f)
{
  return GAUSS_NORM/s*exp(-0.5/(s*s)*(f-f0)*(f-f0));
}

double voigt(double s,double g,double f0,double f)
{
  const double fact = VOIGT_NORM1/s;
  const double x = (f-f0)*fact;
  const double y = g*fact;
  double w,z;
  int flag;
  wofz_(&x,&y,&w,&z,&flag);
  return w*fact*VOIGT_NORM2;
}

int add_lorentz(double gamma,double xc,double norm,int n,const double *x,double *y)
{
  int i;
  double alpha,beta;

  if(gamma < 1.0e-100) // Error
  {
    return -1;
  }
  else
  {
    alpha = LOREN_NORM*norm*gamma;
    beta = gamma*gamma;
    for(i=0; i<n; i++)
    {
      y[i] += alpha/((x[i]-xc)*(x[i]-xc)+beta);
    }
  }

  return 0;
}

int add_gauss(double sigma,double xc,double norm,int n,const double *x,double *y)
{
  int i;
  double alpha,beta;

  if(sigma < 1.0e-100) // Error
  {
    return -1;
  }
  else
  {
    alpha = GAUSS_NORM*norm/sigma;
    beta = -0.5/(sigma*sigma);
    for(i=0; i<n; i++)
    {
      y[i] += alpha*exp(beta*(x[i]-xc)*(x[i]-xc));
    }
  }

  return 0;
}

int add_voigt(double sigma,double gamma,double x0,double area,int np,const double *xp,double *yp)
{
  static const double threshold = 1.0e-50; // Y threshold
  static const int nr = 150;
  static const double r = 1.1;
  static const double width = 5.0e3;
  const double xp_min = xp[0]; // xp must be aligned in ascending order
  const double xp_max = xp[np-1];
  const double xdif_min = (xp_min>=x0?xp_min-x0:(xp_max<=x0?x0-xp_max:0.0));
  const double xdif_max = (xp_min>=x0?xp_max-x0:(xp_max<=x0?x0-xp_min:MAX(xp_max-x0,x0-xp_min)));
  // Voigt profile width approximation by J.J.Olivero and R.L.Longbothum
  // JQSRT, Vol.17, pp.233-236 (1977)
  static const double w1 = (0.5*1.0692);
  static const double w2 = (0.25*0.86639);
  static const double w3 = 1.3862943611198906; // log(2.0)*2.0
  const double hwhm = w1*gamma+sqrt(w2*gamma*gamma+w3*sigma*sigma);
  const double x1 = MAX(hwhm*width,xdif_max)*(r-1.0)/(pow(r,(double)nr)-1.0); // MAX(hwhm*width,xdif_max)*6.2e-8
  const double rx_1 = (r-1.0)/x1;
  const double lr_1 = 1.0/log(r);
  const double sign = (area<0.0?-1.0:1.0);
  static const int nmgn = 16;
  const int nx = nr+1;
  const int nx_min = MAX((int)(floor(log(xdif_min*rx_1+1.0)*lr_1)+0.1)-nmgn,0);
  const int nx_max = MIN((int)(ceil(log(xdif_max*rx_1+1.0)*lr_1)+0.1)+nmgn,nx);
  const int nx_neg = MIN(nmgn,nx_max-nx_min-1);
  const int nx_dif = nx_max-nx_min+(nx_min==0?nx_neg:0);
  const int index = (xp_min>=x0?-1:(xp_max<=x0?np-1:gsl_interp_bsearch(xp,x0,0,np)));
  double xn[nx_dif];
  double yn[nx_dif];
  double dx;
  double sx;
  double alpha,beta;
  double fact1,fact2;
  double yi;
  double zr,zi;
  double wr,wi;
  int i,j;
  int flag;

  dx = x1;
  sx = 0.0;
  for(i=1; i<=nx_min; i++)
  {
    sx += dx;
    dx *= r;
  }
  if(nx_min == 0)
  {
    if(sigma<1.0e-100 && gamma<1.0e-100) // Error
    {
      return -1;
    } else
    if(gamma>1.0e-100 && sigma/gamma<1.0e-10) // Lorentz
    {
      alpha = LOREN_NORM*area*gamma;
      beta = gamma*gamma;
      for(i=nx_min,j=nx_neg; i<nx_max; i++,j++)
      {
        xn[j] = sx;
        sx += dx;
        dx *= r;
        yn[j] = alpha/(xn[j]*xn[j]+beta);
        if(yn[j] < threshold)
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
    } else
    if(sigma>1.0e-100 && gamma/sigma<1.0e-10) // Gauss
    {
      alpha = GAUSS_NORM*area/sigma;
      beta = -0.5/(sigma*sigma);
      for(i=nx_min,j=nx_neg; i<nx_max; i++,j++)
      {
        xn[j] = sx;
        sx += dx;
        dx *= r;
        yn[j] = alpha*exp(beta*xn[j]*xn[j]);
        if(yn[j] < threshold)
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
    else // Voigt
    {
      fact1 = VOIGT_NORM1/sigma;
      fact2 = VOIGT_NORM2*area*fact1;
      zi = gamma*fact1;
      for(i=nx_min,j=nx_neg; i<nx_max; i++,j++)
      {
        xn[j] = sx;
        sx += dx;
        dx *= r;
        zr = xn[j]*fact1;
        wofz_(&zr,&zi,&wr,&wi,&flag);
        yn[j] = wr*fact2;
        if(yn[j] < threshold)
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
    for(i=nx_neg+1,j=nx_neg-1; j>=0; i++,j--)
    {
      xn[j] = -xn[i];
      yn[j] = yn[i];
    }
  }
  else
  {
    if(sigma<1.0e-100 && gamma<1.0e-100) // Error
    {
      return -1;
    } else
    if(gamma>1.0e-100 && sigma/gamma<1.0e-10) // Lorentz
    {
      alpha = LOREN_NORM*area*gamma;
      beta = gamma*gamma;
      for(i=nx_min,j=0; i<nx_max; i++,j++)
      {
        xn[j] = sx;
        sx += dx;
        dx *= r;
        yn[j] = alpha/(xn[j]*xn[j]+beta);
        if(yn[j] < threshold)
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
    } else
    if(sigma>1.0e-100 && gamma/sigma<1.0e-10) // Gauss
    {
      alpha = GAUSS_NORM*area/sigma;
      beta = -0.5/(sigma*sigma);
      for(i=nx_min,j=0; i<nx_max; i++,j++)
      {
        xn[j] = sx;
        sx += dx;
        dx *= r;
        yn[j] = alpha*exp(beta*xn[j]*xn[j]);
        if(yn[j] < threshold)
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
    else // Voigt
    {
      fact1 = VOIGT_NORM1/sigma;
      fact2 = VOIGT_NORM2*area*fact1;
      zi = gamma*fact1;
      for(i=nx_min,j=0; i<nx_max; i++,j++)
      {
        xn[j] = sx;
        sx += dx;
        dx *= r;
        zr = xn[j]*fact1;
        wofz_(&zr,&zi,&wr,&wi,&flag);
        yn[j] = wr*fact2;
        if(yn[j] < threshold)
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

int get_hitran96(char *s,struct hitran96 *h)
{
  int n;
  int err;
  char temp[MAXCHAR];
  char *p;

  strncpy(temp,s,MAXCHAR);
  if((p=strchr(temp,'\r')) != NULL)
  {
    *p = '\0';
  }
  if((p=strchr(temp,'\n')) != NULL)
  {
    *p = '\0';
  }
  if(strlen(temp) != HITRAN96_LINE_LENGTH)
  {
    if(strcmp(s,"\x1a") != 0)
    {
      fprintf(stderr,"Error, length=%ld >>> %s\n",strlen(temp),temp);
    }
    else
    {
      fprintf(stderr,"EOF\n");
    }
    return -1;
  }
  // read values
  err = 0;
  n = 0;
  if(get_int(temp,n, 2,&h->nmol) < 0) err = 1; n +=  2;
  if(get_hex(temp,n, 1,&h->niso) < 0) err = 1; n +=  1;
  if(get_dbl(temp,n,12,&h->freq) < 0) err = 1; n += 12;
  if(get_dbl(temp,n,10,&h->sint) < 0) err = 1; n += 10;
  if(get_dbl(temp,n,10,&h->mtrn) < 0) err = 1; n += 10;
  if(get_dbl(temp,n, 5,&h->wair) < 0) err = 1; n +=  5;
  if(get_dbl(temp,n, 5,&h->wslf) < 0) err = 1; n +=  5;
  if(get_dbl(temp,n,10,&h->elow) < 0) err = 1; n += 10;
  if(get_dbl(temp,n, 4,&h->ctmp) < 0) err = 1; n +=  4;
  if(get_dbl(temp,n, 8,&h->psft) < 0) err = 1; n +=  8;
  if(get_int(temp,n, 3,&h->gqup) < 0) err = 1; n +=  3;
  if(get_int(temp,n, 3,&h->gqlo) < 0) err = 1; n +=  3;
  if(get_chr(temp,n, 9,h->lqup)  < 0) err = 1; n +=  9;
  if(get_chr(temp,n, 9,h->lqlo)  < 0) err = 1; n +=  9;
  if(get_int(temp,n, 1,&h->err1) < 0) err = 1; n +=  1;
  if(get_int(temp,n, 1,&h->err2) < 0) err = 1; n +=  1;
  if(get_int(temp,n, 1,&h->err3) < 0) err = 1; n +=  1;
  if(get_int(temp,n, 2,&h->ref1) < 0) err = 1; n +=  2;
  if(get_int(temp,n, 2,&h->ref2) < 0) err = 1; n +=  2;
  if(get_int(temp,n, 2,&h->ref3) < 0) err = 1; n +=  2;
  if(err)
  {
    fprintf(stderr,"Read error >>> %s\n",temp);
    return -1;
  }

  return 0;
}

int get_hitran04(char *s,struct hitran04 *h)
{
  int n;
  int err;
  char temp[MAXCHAR];
  char *p;

  strncpy(temp,s,MAXCHAR);
  if((p=strchr(temp,'\r')) != NULL)
  {
    *p = '\0';
  }
  if((p=strchr(temp,'\n')) != NULL)
  {
    *p = '\0';
  }
  if(strlen(temp) != HITRAN04_LINE_LENGTH)
  {
    fprintf(stderr,"Error, length=%ld >>> %s\n",strlen(temp),temp);
    return -1;
  }
  // read values
  err = 0;
  n = 0;
  if(get_int(temp,n, 2,&h->nmol) < 0) err = 1; n +=  2;
  if(get_hex(temp,n, 1,&h->niso) < 0) err = 1; n +=  1;
  if(get_dbl(temp,n,12,&h->freq) < 0) err = 1; n += 12;
  if(get_dbl(temp,n,10,&h->sint) < 0) err = 1; n += 10;
  if(get_dbl(temp,n,10,&h->aein) < 0) err = 1; n += 10;
  if(get_dbl(temp,n, 5,&h->wair) < 0) err = 1; n +=  5;
  if(get_dbl(temp,n, 5,&h->wslf) < 0) err = 1; n +=  5;
  if(get_dbl(temp,n,10,&h->elow) < 0) err = 1; n += 10;
  if(get_dbl(temp,n, 4,&h->ctmp) < 0) err = 1; n +=  4;
  if(get_dbl(temp,n, 8,&h->psft) < 0) err = 1; n +=  8;
  if(get_chr(temp,n,15,h->gqup)  < 0) err = 1; n += 15;
  if(get_chr(temp,n,15,h->gqlo)  < 0) err = 1; n += 15;
  if(get_chr(temp,n,15,h->lqup)  < 0) err = 1; n += 15;
  if(get_chr(temp,n,15,h->lqlo)  < 0) err = 1; n += 15;
  if(get_int(temp,n, 1,&h->err1) < 0) err = 1; n +=  1;
  if(get_int(temp,n, 1,&h->err2) < 0) err = 1; n +=  1;
  if(get_int(temp,n, 1,&h->err3) < 0) err = 1; n +=  1;
  if(get_int(temp,n, 1,&h->err4) < 0) err = 1; n +=  1;
  if(get_int(temp,n, 1,&h->err5) < 0) err = 1; n +=  1;
  if(get_int(temp,n, 1,&h->err6) < 0) err = 1; n +=  1;
  if(get_int(temp,n, 2,&h->ref1) < 0) err = 1; n +=  2;
  if(get_int(temp,n, 2,&h->ref2) < 0) err = 1; n +=  2;
  if(get_int(temp,n, 2,&h->ref3) < 0) err = 1; n +=  2;
  if(get_int(temp,n, 2,&h->ref4) < 0) err = 1; n +=  2;
  if(get_int(temp,n, 2,&h->ref5) < 0) err = 1; n +=  2;
  if(get_int(temp,n, 2,&h->ref6) < 0) err = 1; n +=  2;
  if(get_chr(temp,n, 1,h->flag)  < 0) err = 1; n +=  1;
  if(get_dbl(temp,n, 7,&h->swup) < 0) err = 1; n +=  7;
  if(get_dbl(temp,n, 7,&h->swlo) < 0) err = 1; n +=  7;
  if(err)
  {
    fprintf(stderr,"Read error >>> %s\n",temp);
    return -1;
  }

  return 0;
}

int cnv_hitran04(char *s,struct hitran96 *h)
{
  int n;
  int err;
  int itmp;
  double dtmp;
  char temp[MAXCHAR];
  char stmp[MAXCHAR];
  char *p;

  strncpy(temp,s,MAXCHAR);
  if((p=strchr(temp,'\r')) != NULL)
  {
    *p = '\0';
  }
  if((p=strchr(temp,'\n')) != NULL)
  {
    *p = '\0';
  }
  if(strlen(temp) != HITRAN04_LINE_LENGTH)
  {
    fprintf(stderr,"Error, length=%ld >>> %s\n",strlen(temp),temp);
    return -1;
  }
  // read values
  err = 0;
  n = 0;
  if(get_int(temp,n, 2,&h->nmol) < 0) err = 1; n +=  2;
  if(get_hex(temp,n, 1,&h->niso) < 0) err = 1; n +=  1;
  if(get_dbl(temp,n,12,&h->freq) < 0) err = 1; n += 12;
  if(get_dbl(temp,n,10,&h->sint) < 0) err = 1; n += 10;
  if(get_dbl(temp,n,10,&dtmp)    < 0) err = 1; n += 10;
  if(get_dbl(temp,n, 5,&h->wair) < 0) err = 1; n +=  5;
  if(get_dbl(temp,n, 5,&h->wslf) < 0) err = 1; n +=  5;
  if(get_dbl(temp,n,10,&h->elow) < 0) err = 1; n += 10;
  if(get_dbl(temp,n, 4,&h->ctmp) < 0) err = 1; n +=  4;
  if(get_dbl(temp,n, 8,&h->psft) < 0) err = 1; n +=  8;
  if(get_chr(temp,n,15,stmp)     < 0) err = 1; n += 15;
  if(get_chr(temp,n,15,stmp)     < 0) err = 1; n += 15;
  if(get_chr(temp,n,15,h->lqup)  < 0) err = 1; n += 15;
  if(get_chr(temp,n,15,h->lqlo)  < 0) err = 1; n += 15;
  if(get_int(temp,n, 1,&h->err1) < 0) err = 1; n +=  1;
  if(get_int(temp,n, 1,&h->err2) < 0) err = 1; n +=  1;
  if(get_int(temp,n, 1,&h->err3) < 0) err = 1; n +=  1;
  if(get_int(temp,n, 1,&itmp)    < 0) err = 1; n +=  1;
  if(get_int(temp,n, 1,&itmp)    < 0) err = 1; n +=  1;
  if(get_int(temp,n, 1,&itmp)    < 0) err = 1; n +=  1;
  if(get_int(temp,n, 2,&h->ref1) < 0) err = 1; n +=  2;
  if(get_int(temp,n, 2,&h->ref2) < 0) err = 1; n +=  2;
  if(get_int(temp,n, 2,&h->ref3) < 0) err = 1; n +=  2;
  if(get_int(temp,n, 2,&itmp)    < 0) err = 1; n +=  2;
  if(get_int(temp,n, 2,&itmp)    < 0) err = 1; n +=  2;
  if(get_int(temp,n, 2,&itmp)    < 0) err = 1; n +=  2;
  if(get_chr(temp,n, 1,stmp)     < 0) err = 1; n +=  1;
  if(get_dbl(temp,n, 7,&dtmp)    < 0) err = 1; n +=  7;
  if(get_dbl(temp,n, 7,&dtmp)    < 0) err = 1; n +=  7;
  if(err)
  {
    fprintf(stderr,"Read error >>> %s\n",temp);
    return -1;
  }

  return 0;
}

int get_chr(char *str,int start,int size,char *value)
{
  int i,j;

  if(start<0 || size<0)
  {
    fprintf(stderr,"get_chr: invalid input >>> %d,%d\n",start,size);
    return -1;
  } else
  if(start+size >= MAXCHAR)
  {
    fprintf(stderr,"get_chr: out of range >>> %d,%d\n",start,size);
    return -1;
  }

  for(i=start,j=0; i<start+size; i++,j++)
  {
    *(value+j) = *(str+i);
  }
  *(value+j) = '\0';

  return 0;
}

int get_int(char *str,int start,int size,int *value)
{
  int i,j;
  char temp[MAXCHAR];
  char strm[MAXCHAR];
  char *endp;

  if(start<0 || size<0)
  {
    fprintf(stderr,"get_int: invalid input >>> %d,%d\n",start,size);
    return -1;
  } else
  if(start+size >= MAXCHAR)
  {
    fprintf(stderr,"get_int: out of range >>> %d,%d\n",start,size);
    return -1;
  }

  for(i=start,j=0; i<start+size; i++,j++)
  {
    temp[j] = *(str+i);
  }
  temp[j] = '\0';
  strm[0] = '\0';
  sscanf(temp,"%s",strm);
  if(strlen(strm) < 1)
  {
    *value = -1;
  }
  else
  {
    errno = 0;
    *value = strtol(strm,&endp,10);
    if(errno==ERANGE || *endp!='\0')
    {
      fprintf(stderr,"get_int: convert error >>> %s\n",strm);
      return -1;
    }
  }

  return 0;
}

int get_hex(char *str,int start,int size,int *value)
{
  int i,j;
  char temp[MAXCHAR];
  char strm[MAXCHAR];
  char *endp;

  if(start<0 || size<0)
  {
    fprintf(stderr,"get_hex: invalid input >>> %d,%d\n",start,size);
    return -1;
  } else
  if(start+size >= MAXCHAR)
  {
    fprintf(stderr,"get_hex: out of range >>> %d,%d\n",start,size);
    return -1;
  }

  for(i=start,j=0; i<start+size; i++,j++)
  {
    temp[j] = *(str+i);
  }
  temp[j] = '\0';
  strm[0] = '\0';
  sscanf(temp,"%s",strm);
  if(strlen(strm) < 1)
  {
    *value = -1;
  }
  else
  {
    errno = 0;
    *value = strtol(strm,&endp,16);
    if(errno==ERANGE || *endp!='\0')
    {
      fprintf(stderr,"get_hex: convert error >>> %s\n",strm);
      return -1;
    }
  }

  return 0;
}

int get_dbl(char *str,int start,int size,double *value)
{
  int i,j;
  char temp[MAXCHAR];
  char strm[MAXCHAR];
  char *endp;

  if(start<0 || size<0)
  {
    fprintf(stderr,"get_dbl: invalid input >>> %d,%d\n",start,size);
    return -1;
  } else
  if(start+size >= MAXCHAR)
  {
    fprintf(stderr,"get_dbl: out of range >>> %d,%d\n",start,size);
    return -1;
  }

  for(i=start,j=0; i<start+size; i++,j++)
  {
    temp[j] = *(str+i);
  }
  temp[j] = '\0';
  strm[0] = '\0';
  sscanf(temp,"%s",strm);
  if(strlen(strm) < 1)
  {
    *value = NAN;
  }
  else
  {
    errno = 0;
    *value = strtod(strm,&endp);
    if(errno==ERANGE || *endp!='\0')
    {
      fprintf(stderr,"get_dbl: convert error >>> %s\n",strm);
      return -1;
    }
  }

  return 0;
}
