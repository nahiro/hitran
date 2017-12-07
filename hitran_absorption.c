/***********************************************************/
/* HITRAN_ABSORPTION ... Calculate absorption cross        */
/*                       section using HITRAN database.    */
/* Author: N.Manago                                        */
/* $Revision: 171 $                                         */
/* $Date: 2013-12-02 17:50:45 +0900 (Mon, 02 Dec 2013) $   */
/***********************************************************/
#include "hitran_common.c"
#include <getopt.h>
#include "strutil.h"
#include "numutil.h"
#ifndef	NAN
#define	NAN				(nanval())
#endif

#define		TOUT			300.0			// Output temperature in K
#define		POUT			1.0			// Output pressure in atm
#define		ZOUT			NAN			// Output altitude in km
#define		RMIX			NAN			// Volume mixing ratio in mol/mol
#define		RISO			1.0			// Abundance ratio in mol/mol
#define		MASS			NAN			// Mass in g/mol
#define		XMIN			350.0			// X min
#define		XMAX			1050.0			// X max
#define		XSTP			0.01			// X step
#define		XSGM			NAN			// X sigma
#define		WSGM			20.0			// Gaussian line width in sigma
#define		WGAM			1.0e4			// Lorentzian line width in gamma
#define		FAIR			1.0			// Gamma-air factor
#define		NISO			-1			// Isotope ID
#define		XUNI_NM			1			// Wavelength in nm
#define		XUNI_WN			2			// Wavenumber in cm-1
#define		XUNI_HZ			3			// Frequency in Hz
#define		XUNI_GHZ		4			// Frequency in GHz
#define		XUNI			XUNI_NM			// X Unit
#define		HVER			2012			// HITRAN version

double		tout			= TOUT;			// Output temperature in K
double		pout			= POUT;			// Output pressure in atm
double		zout			= ZOUT;			// Output altitude in km
double		rmix			= RMIX;			// Volume mixing ratio in mol/mol
double		riso			= RISO;			// Abundance ratio in mol/mol
double		mass			= MASS;			// Mass in g/mol
double		xmin			= XMIN;			// X min
double		xmax			= XMAX;			// X max
double		xstp			= XSTP;			// X step
double		xsgm			= XSGM;			// X sigma
double		wsgm			= WSGM;			// Gaussian line width in sigma
double		wgam			= WGAM;			// Lorentzian line width in gamma
double		fair			= FAIR;			// Gamma-air factor
double		dwid			= 0.0;			// Doppler width factor
double		*xval			= NULL;			// X array
double		*xabs			= NULL;			// X array
double		*kabs			= NULL;			// Absorption cross-section
int		imax			= 0;			// Buffer size
int		niso			= NISO;			// Isotope ID
int		xuni			= XUNI;			// X Unit
int		hver			= HVER;			// HITRAN version
int		vb			= 0;			// Verbose mode
int		db			= 0;			// Debug   mode
int		hp			= 0;			// Help    mode

int Init(void);
int GetOpt(int argn,char **args);
int Usage(void);
int (*get_hitran)(char *s,struct hitran96 *h);
double pressure(double z);
double temperature(double z);
double air(double z);
double h2o(double z);

int main(int argc,char **argv)
{
  int i,i1,i2;
  double wl,wg;
  double fs,fw,fo;
  double xs,xw,xo,xf;
  char line[MAXCHAR];
  struct hitran96 h;

  if(GetOpt(argc,argv) < 0) return -1;
  if(hp) {Usage(); return 0;}
  if(Init() < 0) return -1;

  if(xuni == XUNI_NM)
  {
    while(fgets(line,MAXCHAR,stdin) != NULL)
    {
      if(get_hitran(line,&h) < 0) continue;
      if(niso > 0)
      {
        if(h.niso != niso) continue;
        h.sint /= riso;
      }
      h.wair *= fair;
      if(vb == 2)
      {
        fprintf(stderr,"%2d %d %13.6f %10.3e %10.3e %6.4f %6.4f "
                       "%11.4f %5.2f %9.6f %3d %3d %9s %9s %d %d %d %2d %2d %2d\n",
                h.nmol,h.niso,h.freq,h.sint,h.mtrn,h.wair,h.wslf,h.elow,h.ctmp,h.psft,
                h.gqup,h.gqlo,h.lqup,h.lqlo,h.err1,h.err2,h.err3,h.ref1,h.ref2,h.ref3);
      }
      fo = origin(pout,&h);
      xo = 1.0e7/fo; // wavelength in nm
      fs = intensity(tout,&h);
      if(vb == 3)
      {
        fprintf(stderr,"%13.6f %13.6e\n",fo,fs);
      }
      if(!isnan(xsgm))
      {
        xs = fs*1.0e-7*xo*xo;
        if(vb == 1)
        {
          if(xo>=xmin && xo<=xmax)
          {
            fprintf(stderr,"%13.6e %13.6e\n",xo,xs);
          }
        }
        wg = GAUSS_FWHM*xsgm*wsgm;
        i1 = (int)((xo-wg-xmin)/xstp);
        i2 = (int)((xo+wg-xmin)/xstp)+1;
        if(i1 < 0) i1 = 0;
        if(i2 > imax) i2 = imax;
        if(i1>=imax || i2<=0) // Out of range
        {
          continue;
        }
        add_gauss(xsgm,xo,xs,i2-i1,&xabs[i1],&kabs[i1]);
      } else
      if(!isnan(mass))
      {
        fw = width(tout,pout,pout*rmix,&h);
        xw = fw*1.0e-7*xo*xo;
        wl = LOREN_FWHM*xw*wgam;
        wg = GAUSS_FWHM*dwid*xo*wsgm;
        if(wg > wl)
        {
          i1 = imax-(int)((xo+wg-xmin)/xstp)-1;
          i2 = imax-(int)((xo-wg-xmin)/xstp);
        }
        else
        {
          i1 = imax-(int)((xo+wl-xmin)/xstp)-1;
          i2 = imax-(int)((xo-wl-xmin)/xstp);
        }
        if(i1 < 0) i1 = 0;
        if(i2 > imax) i2 = imax;
        if(i1>=imax || i2<=0) // Out of range
        {
          continue;
        }
        add_voigt(dwid*fo,fw,fo,fs,i2-i1,&xabs[i1],&kabs[i1]);
      }
      else
      {
        fw = width(tout,pout,pout*rmix,&h);
        xw = fw*1.0e-7*xo*xo;
        wl = LOREN_FWHM*xw*wgam;
        if(db)
        {
          fprintf(stderr,"%13.6e %13.6e %13.6e\n",fo,fs,fw);
        }
        if(vb == 1)
        {
          if(xo>=xmin && xo<=xmax)
          {
            fprintf(stderr,"%13.6e %13.6e %13.6e\n",xo,fs,xw);
          }
        }
        i1 = imax-(int)((xo+wl-xmin)/xstp)-1;
        i2 = imax-(int)((xo-wl-xmin)/xstp);
        if(i1 < 0) i1 = 0;
        if(i2 > imax) i2 = imax;
        if(i1>=imax || i2<=0) // Out of range
        {
          continue;
        }
        add_lorentz(fw,fo,fs,i2-i1,&xabs[i1],&kabs[i1]);
      }
    }
  }
  else
  {
    switch(xuni)
    {
      case XUNI_WN:
        xf = 1.0;
        break;
      case XUNI_HZ:
        xf = C_LIGHT*1.0e2;
        break;
      case XUNI_GHZ:
        xf = C_LIGHT*1.0e-7;
        break;
      default:
        fprintf(stderr,"Error, unknown x unit >>> %d\n",xuni);
        return -1;
        break;
    }
    while(fgets(line,MAXCHAR,stdin) != NULL)
    {
      if(get_hitran(line,&h) < 0) continue;
      if(niso > 0)
      {
        if(h.niso != niso) continue;
        h.sint /= riso;
      }
      h.wair *= fair;
      if(vb == 2)
      {
        fprintf(stderr,"%2d %d %13.6f %10.3e %10.3e %6.4f %6.4f "
                       "%11.4f %5.2f %9.6f %3d %3d %9s %9s %d %d %d %2d %2d %2d\n",
                h.nmol,h.niso,h.freq,h.sint,h.mtrn,h.wair,h.wslf,h.elow,h.ctmp,h.psft,
                h.gqup,h.gqlo,h.lqup,h.lqlo,h.err1,h.err2,h.err3,h.ref1,h.ref2,h.ref3);
      }
      fo = origin(pout,&h);
      xo = fo*xf;
      fs = intensity(tout,&h);
      xs = fs*xf;
      if(vb == 3)
      {
        fprintf(stderr,"%13.6f %13.6e\n",fo,fs);
      }
      if(!isnan(xsgm))
      {
        if(vb == 1)
        {
          if(xo>=xmin && xo<=xmax)
          {
            fprintf(stderr,"%13.6e %13.6e\n",xo,xs);
          }
        }
        wg = GAUSS_FWHM*xsgm*wsgm;
        i1 = (int)((xo-wg-xmin)/xstp);
        i2 = (int)((xo+wg-xmin)/xstp)+1;
        if(i1 < 0) i1 = 0;
        if(i2 > imax) i2 = imax;
        if(i1>=imax || i2<=0) // Out of range
        {
          continue;
        }
        add_gauss(xsgm,xo,xs,i2-i1,&xabs[i1],&kabs[i1]);
      } else
      if(!isnan(mass))
      {
        fw = width(tout,pout,pout*rmix,&h);
        xw = fw*xf;
        wl = LOREN_FWHM*xw*wgam;
        wg = GAUSS_FWHM*dwid*xo*wsgm;
        if(wg > wl)
        {
          i1 = (int)((xo-wg-xmin)/xstp);
          i2 = (int)((xo+wg-xmin)/xstp)+1;
        }
        else
        {
          i1 = (int)((xo-wl-xmin)/xstp);
          i2 = (int)((xo+wl-xmin)/xstp)+1;
        }
        if(i1 < 0) i1 = 0;
        if(i2 > imax) i2 = imax;
        if(i1>=imax || i2<=0) // Out of range
        {
          continue;
        }
        add_voigt(dwid*xo,xw,xo,xs,i2-i1,&xabs[i1],&kabs[i1]);
      }
      else
      {
        fw = width(tout,pout,pout*rmix,&h);
        xw = fw*xf;
        wl = LOREN_FWHM*xw*wgam;
        if(db)
        {
          fprintf(stderr,"%13.6e %13.6e %13.6e\n",fo,fs,fw);
        }
        if(vb == 1)
        {
          if(xo>=xmin && xo<=xmax)
          {
            fprintf(stderr,"%13.6e %13.6e %13.6e\n",xo,fs,xw);
          }
        }
        i1 = (int)((xo-wl-xmin)/xstp);
        i2 = (int)((xo+wl-xmin)/xstp)+1;
        if(i1 < 0) i1 = 0;
        if(i2 > imax) i2 = imax;
        if(i1>=imax || i2<=0) // Out of range
        {
          continue;
        }
        add_lorentz(xw,xo,xs,i2-i1,&xabs[i1],&kabs[i1]);
      }
    }
  }
  if(xuni==XUNI_NM && isnan(xsgm))
  {
    for(i=imax-1; i>=0; i--)
    {
      printf("%11.6f %20.13e\n",xval[i],kabs[i]);
    }
  }
  else
  {
    for(i=0; i<imax; i++)
    {
      printf("%11.6f %20.13e\n",xval[i],kabs[i]);
    }
  }
  free(xval);
  xval = NULL;
  if(xuni==XUNI_NM && isnan(xsgm))
  {
    free(xabs);
    xabs = NULL;
  }
  free(kabs);
  kabs = NULL;

  return 0;
}

int Init(void)
{
  int i,j;

  if(hver >= 2016)
  {
    tips = tips16;
    get_hitran = cnv_hitran04;
  } else
  if(hver >= 2012)
  {
    tips = tips12;
    get_hitran = cnv_hitran04;
  } else
  if(hver >= 2008)
  {
    tips = tips08;
    get_hitran = cnv_hitran04;
  } else
  if(hver >= 2004)
  {
    tips = tips04;
    get_hitran = cnv_hitran04;
  }
  else
  {
    tips = tips96;
    get_hitran = get_hitran96;
    if(get_coeff() < 0)
    {
      return -1;
    }
  }

  if(!isnan(mass))
  {
    dwid = sqrt(K_B*tout/(mass*1.0e-3/N_A*C_LIGHT*C_LIGHT));
  }

  imax = (int)((xmax-xmin)/xstp+0.5);
  if((xval=(double *)malloc(imax*sizeof(double))) == NULL)
  {
    fprintf(stderr,"Error in allocating memory.\n");
    return -1;
  }
  if(xuni==XUNI_NM && isnan(xsgm))
  {
    if((xabs=(double *)malloc(imax*sizeof(double))) == NULL)
    {
      fprintf(stderr,"Error in allocating memory.\n");
      free(xval);
      xval = NULL;
      return -1;
    }
    for(i=0,j=imax-1; i<imax; i++,j--)
    {
      xval[j] = xmin+xstp*i;
      xabs[j] = 1.0e7/xval[j];
    }
  }
  else
  {
    for(i=0; i<imax; i++)
    {
      xval[i] = xmin+xstp*i;
    }
    xabs = xval;
  }
  if((kabs=(double *)malloc(imax*sizeof(double))) == NULL)
  {
    fprintf(stderr,"Error in allocating memory.\n");
    free(xval);
    xval = NULL;
    if(xuni==XUNI_NM && isnan(xsgm))
    {
      free(xabs);
      xabs = NULL;
    }
    return -1;
  }
  for(i=0; i<imax; i++)
  {
    kabs[i] = 0.0;
  }

  return 0;
}

double pressure(double z)
{
  return 1.00527*exp(-z/7.8364); // atm
}

double temperature(double z)
{
  if(z < 12.25)
  {
    return 297.2-6.105*z;
  }
  else
  {
    return 310.38-12.88*z+0.56183*z*z-8.3034e-3*z*z*z+3.7879e-5*z*z*z*z; // K
  }
}

double air(double z)
{
  return 2.55e19*exp(-z/8.5098); // cm-3
}

double h2o(double z)
{
  return 5.3703e17*exp(-z/1.8325); // cm-3
}

int GetOpt(int argn,char **args)
{
  int c,rt;
  int option_index = 0;
  int this_option_optind;
  int ntmp;
  double xtmp;
  char *p;
  struct option long_options[] =
  {
    {"dnam",1,0,'f'},
    {"tout",1,0,'o'},
    {"pout",1,0,'p'},
    {"zout",1,0,'z'},
    {"rmix",1,0,'R'},
    {"riso",1,0,'r'},
    {"niso",1,0,'N'},
    {"mass",1,0,'M'},
    {"xmin",1,0,'x'},
    {"xmax",1,0,'X'},
    {"xstp",1,0,'D'},
    {"xsgm",1,0,'S'},
    {"wsgm",1,0,'w'},
    {"wgam",1,0,'W'},
    {"fair",1,0,'F'},
    {"xuni",1,0,'U'},
    {"hver",1,0,'n'},
    {"debug",0,0,'d'},
    {"verbose",0,0,'v'},
    {"help",0,0,'h'}
  };

  rt = 0;
  while(1)
  {
    this_option_optind = optind?optind:1;
    c = getopt_long(argn,args,":f:o:p:z:R:r:N:M:x:X:D:S:w:W:F:U:n:dvh",long_options,&option_index);
    if(c == -1) break;

    switch(c)
    {
      case 'f':
        strncpy(dnam,optarg,MAXCHAR);
        break;
      case 'o':
        errno = 0;
        xtmp = strtod(optarg,&p);
        if(errno!=ERANGE && *p=='\0' && xtmp>=70.0 && xtmp<=3005.0) tout = xtmp;
        else
        {
          fprintf(stderr,"Temperature -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'p':
        errno = 0;
        xtmp = strtod(optarg,&p);
        if(errno!=ERANGE && *p=='\0' && xtmp>0.0 && xtmp<10.0) pout = xtmp;
        else
        {
          fprintf(stderr,"Pressure -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'z':
        errno = 0;
        xtmp = strtod(optarg,&p);
        if(errno!=ERANGE && *p=='\0' && xtmp>=0.0) zout = xtmp;
        else
        {
          fprintf(stderr,"Altitude -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'R':
        errno = 0;
        xtmp = strtod(optarg,&p);
        if(errno!=ERANGE && *p=='\0' && xtmp>0.0) rmix = xtmp;
        else
        {
          fprintf(stderr,"Mixing ratio -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'r':
        errno = 0;
        xtmp = strtod(optarg,&p);
        if(errno!=ERANGE && *p=='\0' && xtmp>0.0) riso = xtmp;
        else
        {
          fprintf(stderr,"Abundance ratio -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'N':
        errno = 0;
        ntmp = strtol(optarg,&p,10);
        if(errno!=ERANGE && *p=='\0' && ntmp>0) niso = ntmp;
        else
        {
          fprintf(stderr,"Isotope ID -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'M':
        errno = 0;
        xtmp = strtod(optarg,&p);
        if(errno!=ERANGE && *p=='\0' && xtmp>0.0) mass = xtmp;
        else
        {
          fprintf(stderr,"Mass -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'x':
        errno = 0;
        xtmp = strtod(optarg,&p);
        if(errno!=ERANGE && *p=='\0' && xtmp>=0.0) xmin = xtmp;
        else
        {
          fprintf(stderr,"X min -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'X':
        errno = 0;
        xtmp = strtod(optarg,&p);
        if(errno!=ERANGE && *p=='\0' && xtmp>0.0) xmax = xtmp;
        else
        {
          fprintf(stderr,"X max -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'D':
        errno = 0;
        xtmp = strtod(optarg,&p);
        if(errno!=ERANGE && *p=='\0' && xtmp>0.0) xstp = xtmp;
        else
        {
          fprintf(stderr,"X step -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'S':
        errno = 0;
        xtmp = strtod(optarg,&p);
        if(errno!=ERANGE && *p=='\0' && xtmp>0.0) xsgm = xtmp;
        else
        {
          fprintf(stderr,"X sigma -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'w':
        errno = 0;
        xtmp = strtod(optarg,&p);
        if(errno!=ERANGE && *p=='\0' && xtmp>0.0 && xtmp<=1.0e10) wsgm = xtmp;
        else
        {
          fprintf(stderr,"Gaussian line width -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'W':
        errno = 0;
        xtmp = strtod(optarg,&p);
        if(errno!=ERANGE && *p=='\0' && xtmp>0.0 && xtmp<=1.0e10) wgam = xtmp;
        else
        {
          fprintf(stderr,"Lorentzian line width -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'F':
        errno = 0;
        xtmp = strtod(optarg,&p);
        if(errno!=ERANGE && *p=='\0' && xtmp>0.0 && xtmp<=1.0e10) fair = xtmp;
        else
        {
          fprintf(stderr,"Gamma-air factor -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'U':
        errno = 0;
        ntmp = strtol(optarg,&p,10);
        if(errno!=ERANGE && *p=='\0' && ntmp>=1 && ntmp<=4) xuni = ntmp;
        else
        {
          fprintf(stderr,"X unit -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'n':
        errno = 0;
        ntmp = strtol(optarg,&p,10);
        if(errno!=ERANGE && *p=='\0' && ntmp>=1996 && ntmp<=2016) hver = ntmp;
        else
        {
          fprintf(stderr,"HITRAN version -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'd':
        db++;
        break;
      case 'v':
        vb++;
        break;
      case 'h':
        hp = 1;
        break;
      case '?':
        if(optopt == '\0')
        {
          fprintf(stderr,"Invalid option : %s\n",args[this_option_optind]);
        }
        else
        {
          fprintf(stderr,"Invalid option : %c\n",optopt);
        }
        rt = -1;
        break;
      case ':':
        if(optopt == '\0')
        {
          fprintf(stderr,"Option requires an argument : %s\n",args[this_option_optind]);
        }
        else
        {
          fprintf(stderr,"Option requires an argument : %c\n",optopt);
        }
        rt = -1;
        break;
      default:
        fprintf(stderr,"?? getopt returned character code 0%o ??\n",c);
        rt = -1;
        break;
    }
  }

  if(!isnan(zout))
  {
    tout = temperature(zout);
    pout = pressure(zout);
  }
  if(isnan(rmix))
  {
    if(isnan(zout))
    {
      zout = 0.0;
    }
    rmix = h2o(zout)/air(zout);
  }
  if(vb>1 || hp)
  {
    if(!isnan(zout))
    {
      fprintf(stderr,"Altitude    : %13.2f km\n",zout);
    }
    fprintf(stderr,"Temperature : %13.2f K\n",tout);
    fprintf(stderr,"Pressure    : %13.5f atm\n",pout);
    fprintf(stderr,"Mixing ratio: %13.6e mol/mol\n",rmix);
  }

  if(optind < argn)
  {
    fprintf(stderr,"non-option ARGV-elements:\n");
    while(optind < argn)
    {
      fprintf(stderr,"%s\n",args[optind++]);
    }
    rt = -1;
  }

  if(hp) rt = 1;
  return rt;
}

int Usage(void)
{
  int n = 15;
  char e[MAXCHAR];
  char a[MAXCHAR];
  char d[MAXCHAR];

  fprintf(stderr,"Usage:\n");
  fprintf(stderr,"hitran_absorption -[option] (argument) -[option] (argument) ...\n");
  fprintf(stderr,"-----------------------------------------------------------------------------\n");
  fprintf(stderr,"   option   |%s|%s|%s| current\n",As(e,"",n),         As(a,"argument",n),   As(d,"default",n));
  fprintf(stderr," f -dnam    |%s|%s|%s| %s\n",As(e,"Data file",n),     As(a,"name",n),       As(d,DNAM,n),dnam);
  fprintf(stderr," o -tout    |%s|%s|%s| %e\n",As(e,"Temperature",n),   As(a,"K",n),          Af(d,TOUT,n),tout);
  fprintf(stderr," p -pout    |%s|%s|%s| %e\n",As(e,"Pressure",n),      As(a,"atm",n),        Af(d,POUT,n),pout);
  fprintf(stderr," z -zout    |%s|%s|%s| %e\n",As(e,"Altitude",n),      As(a,"km",n),         Af(d,0.0,n),zout);
  fprintf(stderr," R -rmix    |%s|%s|%s| %e\n",As(e,"Mixing ratio ",n), As(a,"mol/mol",n),    Ae(d,h2o(0)/air(0),n),rmix);
  fprintf(stderr," r -riso    |%s|%s|%s| %e\n",As(e,"Abundance ratio ",n),As(a,"mol/mol",n),  Ae(d,RISO,n),riso);
  fprintf(stderr," N -niso    |%s|%s|%s| %d\n",As(e,"Isotope ID",n),    As(a,"#",n),          Ad(d,NISO,n),niso);
  fprintf(stderr," M -mass    |%s|%s|%s| %e\n",As(e,"Mass",n),          As(a,"g/mol",n),      Af(d,MASS,n),mass);
  fprintf(stderr," x -xmin    |%s|%s|%s| %f\n",As(e,"X min  ",n),       As(a,"X unit   ",n),  Af(d,XMIN,n),xmin);
  fprintf(stderr," X -xmax    |%s|%s|%s| %f\n",As(e,"X max  ",n),       As(a,"X unit   ",n),  Af(d,XMAX,n),xmax);
  fprintf(stderr," D -xstp    |%s|%s|%s| %f\n",As(e,"X step ",n),       As(a,"X unit(*)",n),  Af(d,XSTP,n),xstp);
  fprintf(stderr," S -xsgm    |%s|%s|%s| %f\n",As(e,"X sigma",n),       As(a,"X unit   ",n),  Af(d,XSGM,n),xsgm);
  fprintf(stderr," w -wsgm    |%s|%s|%s| %f\n",As(e,"Gauss width",n),   As(a,"sigma",n),      Af(d,WSGM,n),wsgm);
  fprintf(stderr," W -wgam    |%s|%s|%s| %f\n",As(e,"Lorentz width",n), As(a,"gamma",n),      Af(d,WGAM,n),wgam);
  fprintf(stderr," F -fair    |%s|%s|%s| %f\n",As(e,"G-air factor",n),  As(a,"factor",n),     Af(d,FAIR,n),fair);
  fprintf(stderr," U -xuni    |%s|%s|%s| %d\n",As(e,"X Unit",n),        As(a,"Unit#",n),      Ad(d,XUNI,n),xuni);
  fprintf(stderr," n -hver    |%s|%s|%s| %d\n",As(e,"HITRAN version",n),As(a,"version",n),    Ad(d,HVER,n),hver);
  fprintf(stderr," d -debug   |%s|%s|%s| %d\n",As(e,"Debug   mode",n),  As(a,"nothing",n),    Ad(d,0,n),db);
  fprintf(stderr," v -verbose |%s|%s|%s| %d\n",As(e,"Verbose mode",n),  As(a,"nothing",n),    Ad(d,0,n),vb);
  fprintf(stderr," h -help    |%s|%s|%s| %d\n",As(e,"Help    mode",n),  As(a,"nothing",n),    Ad(d,0,n),1);
  fprintf(stderr,"-----------------------------------------------------------------------------\n");
  fprintf(stderr,"The zout option is used to calculate temperature, pressure, and partial pressure of H2O.\n");
  fprintf(stderr,"They can be set individually using the tout, pout, and rmix options.\n");
  fprintf(stderr,"Note: use the rmix option for molecules other than H2O.\n");
  fprintf(stderr,"The xsgm option is used to calculate Gauss profiles (default: Lorentz profile).\n");
  fprintf(stderr,"X Unit 1:nm 2:cm-1 3:Hz 4:GHz\n");

  return 0;
}
