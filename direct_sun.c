/***********************************************************/
/* DIRECT_SUN ... direct solar irradiance simulation       */
/* Author: N.Manago                                        */
/* $Revision: 162 $                                         */
/* $Date: 2013-12-02 13:22:33 +0900 (Mon, 02 Dec 2013) $   */
/***********************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <getopt.h>
#include <math.h>
//#include <bits/nan.h>
#include "strutil.h"
#include "numutil.h"
#ifndef	NAN
#define	NAN				(nanval())
#endif
#include "hitran_common.c"

#define		SOLAR_CONST		1370.0			// Solar constant in W/m2
#define		SOLAR_TEMP		5780.0			// Solar temperature in K
#define		KOSHM_COEFF		3.912023		// Koshmieder coefficient
#define		SIGRAY_550		4.30e-31		// Rayleigh cross section in m2
#define		NCOL_AIR		2.17e29			// Column air density at ground in m-2
#define		NVOL_AIR		2.55e25			// Volume air density at ground in m-3
#define		HITRAN_NMOL_H2O		1			// H2O
#define		HITRAN_NMOL_O2		7			// O2
#define		HVER			2012			// HITRAN version
#define		MIXR_O2			0.21			// O2 Mixing ratio
#define		HRAY			(NCOL_AIR/NVOL_AIR)	// Rayleigh scale height in m
#define		VMIE			20.0e3			// Visibility in m
#define		HMIE			2.0e3			// Mie scale height in m
#define		PMIE			1.0			// Angstrom exponent
#define		WSCL			1.0			// H2O scale
#define		XMIN			350.0			// Min wavelength in nm
#define		XMAX			1050.0			// Max wavelength in nm
#define		XSTP			0.01			// Wavelength step in nm
#define		SSTP			1.0			// Wavelength step in nm
#define		LWID			1.0			// Line width in nm
#define		RWID			50.0			// Resolution in nm
#define		SALT			(70.0*D_TO_R)		// Solar altitude in radian
#define		MAXSTEP			100000
#define		MAXSDAT			1000
#define		MAXDATA			100000
#define		DELTA			1.0e-3

double		xmin			= XMIN;
double		xmax			= XMAX;
double		xstp			= XSTP;
double		sstp			= SSTP;
double		lwid			= LWID;
double		rwid			= RWID;
double		vmie			= VMIE;
double		hmie			= HMIE;
double		pmie			= PMIE;
double		wscl			= WSCL;
double		salt			= SALT;
double		th_s;
double		cth_s;
int		imax			= 0;
int		ndat			= 0;
int		hver			= HVER;			// HITRAN version
int		vb			= 0;			// Verbose mode
int		db			= 0;			// Debug   mode
int		hp			= 0;			// Help    mode
char		snam[MAXCHAR]		= "";
struct	hitran96 data[MAXDATA];

// Optional data
double		sun_flux[MAXSDAT];
int		nsun;

int Init(void);
int GetOpt(int argn,char **args);
int Usage(void);
int ReadHITRAN(void);
int ReadSUN(void);
double pressure(double z);
double temperature(double z);
double air(double z);
double h2o(double z);
double alp_ray(double l,double z);
double alp_mie(double hs,double vis,double p,double l,double z);
double bbr(double l,double t);
double sun(double l);

int main(int argc,char **argv)
{
  int i,n;
  double v;
  double x,xo;
  double f,fs,fw,fo;
  double zmin = 0.0;
  double zmax = 10.0e3;
  double zstp = 10.0;
  double zout,tout,pout;
  double dens;
  double tau[MAXSTEP];
  double val[MAXSTEP];

  if(GetOpt(argc,argv) < 0) exit(-1);
  if(hp) {Usage(); return 0;}
  if(Init() < 0) exit(-1);

  for(i=0; i<imax; i++)
  {
    x = xmin+xstp*i;
    if(nsun > 0)
    {
      val[i] = sun(x);
    }
    else
    {
      val[i] = bbr(x,SOLAR_TEMP);
    }
  }

  for(zout=zmax; zout>zmin; zout-=zstp)
  {
    fprintf(stderr,"z=%13.4e\n",zout);
    for(i=0; i<imax; i++)
    {
      tau[i] = 0.0;
    }
    tout = temperature(zout);
    pout = pressure(zout);
    for(n=0; n<ndat; n++)
    {
      switch(data[n].nmol)
      {
        case HITRAN_NMOL_H2O:
          v = wscl*h2o(zout);
          fw = width(tout,pout,pout*v/air(zout),&data[n]);
          dens = v*1.0e2; // molecule/m/cm2
          break;
        case HITRAN_NMOL_O2:
          fw = width(tout,pout,pout*MIXR_O2,&data[n]);
          dens = air(zout)*MIXR_O2*1.0e2; // molecule/m/cm2
          break;
        default:
          fprintf(stderr,"Unsupported molecule >>> %d\n",data[n].nmol);
          return -1;
          break;
      }
      fs = intensity(tout,&data[n])*dens;
      fo = origin(pout,&data[n]);
      xo = 1.0e7/fo;
      if(xo>=xmin && xo<=xmax)
      {
        for(x=xo-lwid; x<=xo+lwid; x+=xstp)
        {
          i = (int)((x+xstp*DELTA-xmin)/xstp);
          if(i<0 || i>=MAXSTEP) continue;
          f = 1.0e7/x;
          tau[i] += fs*lorentz(fw,fo,f);
        }
      }
    }
    for(i=0; i<imax; i++)
    {
      x = xmin+xstp*i;
      val[i] -= val[i]*(alp_ray(x,zout)+alp_mie(hmie,vmie,pmie,x,zout)+tau[i])*zstp/cth_s;
    }
  }

  for(i=0; i<imax; i++)
  {
    x = xmin+xstp*i;
    printf("%8.2f %13.4e\n",x,val[i]);
  }

  return 0;
}

int Init(void)
{
  th_s = PI_2-salt;
  cth_s = cos(th_s);

  if(hver >= 2004)
  {
    tips = tips04;
  }
  else
  {
    tips = tips96;
    if(get_coeff() < 0)
    {
      return -1;
    }
  }

  imax = (int)((xmax-xmin)/xstp);
  if(imax >= MAXSTEP)
  {
    fprintf(stderr,"Error, #step exceed the limit >>> %d\n",imax);
    return -1;
  }

  if(ReadHITRAN() < 0) return -1;
  if(strcmp(snam,"") != 0)
  {
    if(ReadSUN() < 0) return -1;
  }

  return 0;
}

int ReadHITRAN(void)
{
  int i;
  double xo;
  char line[MAXCHAR];
  struct hitran96 h;

  i = 0;
  while(fgets(line,MAXCHAR,stdin) != NULL)
  {
    if(hver >= 2004)
    {
      if(cnv_hitran04(line,&h) < 0) return -1;
    }
    else
    {
      if(get_hitran96(line,&h) < 0) return -1;
    }
    if(vb > 1)
    {
      fprintf(stderr,"%2d %d %13.6f %10.3e %10.3e %6.4f %6.4f "
                     "%11.4f %5.2f %9.6f %3d %3d %9s %9s %d %d %d %2d %2d %2d\n",
              h.nmol,h.niso,h.freq,h.sint,h.mtrn,h.wair,h.wslf,h.elow,h.ctmp,h.psft,
              h.gqup,h.gqlo,h.lqup,h.lqlo,h.err1,h.err2,h.err3,h.ref1,h.ref2,h.ref3);
    }
    xo = 1.0e7/h.freq;
    if(xo>=xmin && xo<=xmax)
    {
      if(i >= MAXDATA)
      {
        fprintf(stderr,"Error, #data exceed the limit %d.\n",i);
        return -1;
      }
      data[i] = h;
      i++;
    }
  }
  ndat = i;

  return 0;
}

double pressure(double z)
{
  return 1.00527*exp(-z/7.8364e+3); // atm
}

double temperature(double z)
{
  if(z < 12.25e3)
  {
    return 297.2-6.105e-3*z;
  }
  else
  {
    return 310.38-12.88e-3*z+0.56183e-6*z*z-8.3034e-12*z*z*z+3.7879e-17*z*z*z*z; // K
  }
}

double air(double z)
{
  return 2.55e19*exp(-z/8.5098e+3); // cm-3
}

double h2o(double z)
{
  return 5.3703e17*exp(-z/1.8325e+3); // cm-3
}

// Rayleigh extinction coefficient
double alp_ray(double l,double z)
{
  return SIGRAY_550*pow(l/550.0,-4.0)*NVOL_AIR*exp(-z/HRAY);
}

// Mie extinction coefficient
double alp_mie(double hs,double vis,double p,double l,double z)
{
  return KOSHM_COEFF/vis*pow(l/550.0,-p)*exp(-z/hs);
}

// Blackbody radiation
double bbr(double l,double t)
{
  double f;
  double x;
  double s;

  x = H_PLANK*C_LIGHT/(l*1.0e-9*K_B*t);
  s = x/PI;
  f = 15.0*SOLAR_CONST*s*s*s*s/(l*(exp(x)-1.0));

  return f;
}

int ReadSUN(void)
{
  int i;
  int nerr;
  double w,v;
  char line[MAXCHAR];
  char str1[MAXCHAR];
  char str2[MAXCHAR];
  char *p;
  FILE *fp;

  if((fp=fopen(snam,"r")) == NULL)
  {
    fprintf(stderr,"Error, cannot open %s\n",snam);
    return -1;
  }
  i = 0;
  nerr = 0;
  while(fgets(line,MAXCHAR,fp) != NULL)
  {
    if(sscanf(line,"%s%s",str1,str2) != 2)
    {
      fprintf(stderr,"Read error >>> %s\n",line);
      nerr = 1;
      break;
    }
    errno = 0;
    w = strtod(str1,&p);
    if(errno==ERANGE || *p!='\0')
    {
      fprintf(stderr,"Convert error >>> %s\n",line);
      nerr = 1;
      break;
    }
    if(w<xmin || w>xmax)
    {
      continue;
    }
    if(fabs(xmin+sstp*i-w) > EPSILON)
    {
      continue;
    }
    errno = 0;
    v = strtod(str2,&p);
    if(errno==ERANGE || *p!='\0')
    {
      fprintf(stderr,"Convert error >>> %s\n",line);
      nerr = 1;
      break;
    }
    if(i >= MAXSDAT)
    {
      fprintf(stderr,"#data exceed the limit >>> %d (%s)\n",i,snam);
      nerr = 1;
      break;
    }
    sun_flux[i] = v*1.0e-3;
    i++;
  }
  fclose(fp);
  nsun = i;
  if(nerr > 0)
  {
    fprintf(stderr,"Failed in reading wavelength.\n");
    return -1;
  }
  if(vb)
  {
    fprintf(stderr,"Solar flux:\n");
    for(i=0; i<nsun; i++)
    {
      fprintf(stderr,"%3d %8.2f %8.2f\n",i+1,xmin+sstp*i,sun_flux[i]);
    }
  }

  return 0;
}

double sun(double l)
{
  int i1,i2;
  double x1;

  if(l<xmin || l>xmax) return 0.0;
  i1 = (int)((l+xstp*DELTA-xmin)/sstp); // assuming dl = xstp
  i2 = i1+1;
  if(i2 >= nsun)
  {
    i2 = nsun-1;
    i1 = i2-1;
  }
  x1 = xmin+sstp*i1;

  return sun_flux[i1]+(sun_flux[i2]-sun_flux[i1])*(l-x1)/(sstp);
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
    {"snam",1,0,'S'},
    {"xmin",1,0,'x'},
    {"xmax",1,0,'X'},
    {"xstp",1,0,'D'},
    {"hmie",1,0,'H'},
    {"pmie",1,0,'P'},
    {"wscl",1,0,'s'},
    {"lwid",1,0,'w'},
    {"rwid",1,0,'W'},
    {"salt",1,0,'a'},
    {"hver",1,0,'n'},
    {"debug",0,0,'d'},
    {"verbose",0,0,'v'},
    {"help",0,0,'h'}
  };

  rt = 0;
  while(1)
  {
    this_option_optind = optind?optind:1;
    c = getopt_long(argn,args,":f:S:x:X:D:H:P:s:w:W:a:n:dvh",long_options,&option_index);
    if(c == -1) break;

    switch(c)
    {
      case 'f':
        strncpy(dnam,optarg,MAXCHAR);
        break;
      case 'S':
        strncpy(snam,optarg,MAXCHAR);
        break;
      case 'x':
        errno = 0;
        xtmp = strtod(optarg,&p);
        if(errno!=ERANGE && *p=='\0' && xtmp>=100.0 && xtmp<=10000.0) xmin = xtmp;
        else
        {
          fprintf(stderr,"Min wavelength -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'X':
        errno = 0;
        xtmp = strtod(optarg,&p);
        if(errno!=ERANGE && *p=='\0' && xtmp>=100.0 && xtmp<=10000.0) xmax = xtmp;
        else
        {
          fprintf(stderr,"Max wavelength -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'D':
        errno = 0;
        xtmp = strtod(optarg,&p);
        if(errno!=ERANGE && *p=='\0' && xtmp>0.0 && xtmp<=10000.0) xstp = xtmp;
        else
        {
          fprintf(stderr,"Wavelength step -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'H':
        errno = 0;
        xtmp = strtod(optarg,&p);
        if(errno!=ERANGE && *p=='\0' && xtmp>0.0) hmie = xtmp;
        else
        {
          fprintf(stderr,"Mie scale hight -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'P':
        errno = 0;
        xtmp = strtod(optarg,&p);
        if(errno!=ERANGE && *p=='\0') pmie = xtmp;
        else
        {
          fprintf(stderr,"Angstrom parameter -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 's':
        errno = 0;
        xtmp = strtod(optarg,&p);
        if(errno!=ERANGE && *p=='\0') wscl = xtmp;
        else
        {
          fprintf(stderr,"H2O scale -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'w':
        errno = 0;
        xtmp = strtod(optarg,&p);
        if(errno!=ERANGE && *p=='\0' && xtmp>0.0 && xtmp<=10000.0) lwid = xtmp;
        else
        {
          fprintf(stderr,"Line width -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'W':
        errno = 0;
        xtmp = strtod(optarg,&p);
        if(errno!=ERANGE && *p=='\0' && xtmp>0.0 && xtmp<=10000.0) rwid = xtmp;
        else
        {
          fprintf(stderr,"Resolution -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'a':
        errno = 0;
        xtmp = strtod(optarg,&p);
        if(errno!=ERANGE && *p=='\0' && xtmp>=0.0 && xtmp<=90.0) salt = xtmp*D_TO_R;
        else
        {
          fprintf(stderr,"Solar altitude -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'n':
        errno = 0;
        ntmp = strtol(optarg,&p,10);
        if(errno!=ERANGE && *p=='\0' && ntmp>=1996 && ntmp<=2012) hver = ntmp;
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
  fprintf(stderr,"direct_sun -[option] (argument) -[option] (argument) ...\n");
  fprintf(stderr,"-----------------------------------------------------------------------------\n");
  fprintf(stderr,"   option   |%s|%s|%s| current\n",As(e,"",n),         As(a,"argument",n),   As(d,"default",n));
  fprintf(stderr," f -dnam    |%s|%s|%s| %s\n",As(e,"Data file",n),     As(a,"name",n),       As(d,DNAM,n),dnam);
  fprintf(stderr," S -snam    |%s|%s|%s| %s\n",As(e,"Solar file",n),    As(a,"name",n),       As(d,"",n),snam);
  fprintf(stderr," x -xmin    |%s|%s|%s| %f\n",As(e,"Min wavelength",n),As(a,"nm",n),         Af(d,XMIN,n),xmin);
  fprintf(stderr," X -xmax    |%s|%s|%s| %f\n",As(e,"Max wavelength",n),As(a,"nm",n),         Af(d,XMAX,n),xmax);
  fprintf(stderr," D -xstp    |%s|%s|%s| %f\n",As(e,"Step wavelen",n),  As(a,"nm",n),         Af(d,XSTP,n),xstp);
  fprintf(stderr," H -hmie    |%s|%s|%s| %e\n",As(e,"Mie scale h",n),   As(a,"m",n),          Ae(d,HMIE,n),hmie);
  fprintf(stderr," P -pmie    |%s|%s|%s| %f\n",As(e,"Angstrom par",n),  As(a,"value",n),      Af(d,PMIE,n),pmie);
  fprintf(stderr," s -wscl    |%s|%s|%s| %e\n",As(e,"H2O scale",n),     As(a,"ratio",n),      Ae(d,WSCL,n),wscl);
  fprintf(stderr," w -lwid    |%s|%s|%s| %f\n",As(e,"Line width",n),    As(a,"nm",n),         Af(d,LWID,n),lwid);
  fprintf(stderr," W -rwid    |%s|%s|%s| %f\n",As(e,"Resolution",n),    As(a,"nm",n),         Af(d,RWID,n),rwid);
  fprintf(stderr," a -salt    |%s|%s|%s| %e\n",As(e,"Solar altitude",n),As(a,"degree",n),     Af(d,SALT*R_TO_D,n),salt*R_TO_D);
  fprintf(stderr," n -hver    |%s|%s|%s| %d\n",As(e,"HITRAN version",n),As(a,"version",n),    Ad(d,HVER,n),hver);
  fprintf(stderr," d -debug   |%s|%s|%s| %d\n",As(e,"Debug   mode",n),  As(a,"nothing",n),    Ad(d,0,n),db);
  fprintf(stderr," v -verbose |%s|%s|%s| %d\n",As(e,"Verbose mode",n),  As(a,"nothing",n),    Ad(d,0,n),vb);
  fprintf(stderr," h -help    |%s|%s|%s| %d\n",As(e,"Help    mode",n),  As(a,"nothing",n),    Ad(d,0,n),1);
  fprintf(stderr,"-----------------------------------------------------------------------------\n");

  return 0;
}
