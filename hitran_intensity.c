/***********************************************************/
/* HITRAN_INTENSITY  ... Calculate line intensity          */
/*                       using HITRAN database.            */
/* Author: N.Manago                                        */
/* $Revision: 45 $                                         */
/* $Date: 2009-12-04 17:11:03 +0900 (Fri, 04 Dec 2009) $   */
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

#define		IMOL			-1			// Mol#
#define		NISO			-1			// Iso#
#define		XUNI_NM			1			// Wavelength in nm
#define		XUNI_WN			2			// Wavenumber in cm-1
#define		XUNI_HZ			3			// Frequency in Hz
#define		XUNI_GHZ		4			// Frequency in GHz
#define		XUNI			XUNI_NM			// X Unit
#define		YUNI_H			1			// HITRAN unit
#define		YUNI_J			2			// JPL unit
#define		YUNI_LOG		3			// JPL unit (log10)
#define		YUNI			YUNI_H			// Y Unit
#define		YFAC_H2J		2.99792458e18		// HITRAN->JPL conversion factor
#define		YNRM			1.0			// Y Weight
#define		RISO			1.0			// Abundance ratio in mol/mol
#define		XMIN			350.0			// X Min
#define		XMAX			1050.0			// X Max
#define		TMIN			200.0			// T Min in K
#define		TMAX			300.0			// T Max in K
#define		TSTP			0.1			// T Step in K
#define		HVER			2012			// HITRAN version

int		nmol			= IMOL;			// Mol#
int		niso			= NISO;			// Iso#
double		riso			= RISO;			// Abundance ratio in mol/mol
double		xmin			= XMIN;			// X Min
double		xmax			= XMAX;			// X Max
double		xfac			= 1.0;			// X Factor
double		yfac			= 1.0;			// Y Factor
double		ynrm			= YNRM;			// Y Weight
double		tmin			= TMIN;			// T Min in K
double		tmax			= TMAX;			// T Max in K
double		tstp			= TSTP;			// T Step in K
int		xuni			= XUNI;			// X Unit
int		yuni			= YUNI;			// Y Unit
int		hver			= HVER;			// HITRAN version
int		lm			= 0;			// Long    mode
int		db			= 0;			// Debug   mode
int		vb			= 0;			// Verbose mode
int		hp			= 0;			// Help    mode

int Init(void);
int GetOpt(int argn,char **args);
int Usage(void);

int main(int argc,char **argv)
{
  double x,t;
  char line[MAXCHAR];
  struct hitran96 h;
  int (*get_hitran)(char *s,struct hitran96 *h);

  if(GetOpt(argc,argv) < 0) return -1;
  if(hp) {Usage(); return 0;}
  if(Init() < 0) return -1;

  if(hver >= 2004)
  {
    get_hitran = cnv_hitran04;
  }
  else
  {
    get_hitran = get_hitran96;
  }
  while(fgets(line,MAXCHAR,stdin) != NULL)
  {
    if(get_hitran(line,&h) < 0) continue;
    if(nmol>0 && h.nmol!=nmol)
    {
      continue;
    }
    if(niso > 0)
    {
      if(h.niso != niso) continue;
      h.sint /= riso;
    }
    if(xuni == XUNI_NM)
    {
      x = 1.0e7/h.freq;
    }
    else
    {
      x = h.freq*xfac;
    }
    if(x<xmin || x>xmax)
    {
      continue;
    }
    if(vb == 1)
    {
      fprintf(stderr,"%3d %2d %13.4f %13.6e\n",h.nmol,h.niso,x,h.sint*yfac);
    } else
    if(vb == 2)
    {
      fprintf(stderr,"%2d %d %13.6f %10.3e %10.3e %6.4f %6.4f "
                     "%11.4f %5.2f %9.6f %3d %3d %9s %9s %d %d %d %2d %2d %2d\n",
              h.nmol,h.niso,h.freq,h.sint,h.mtrn,h.wair,h.wslf,h.elow,h.ctmp,h.psft,
              h.gqup,h.gqlo,h.lqup,h.lqlo,h.err1,h.err2,h.err3,h.ref1,h.ref2,h.ref3);
    }
    if(lm)
    {
      if(tmax < tmin)
      {
        if(yuni == YUNI_LOG)
        {
          printf("%17.8f %13.6e %4d\n",x,log10(intensity(tmin,&h)*yfac),h.niso);
        }
        else
        {
          printf("%17.8f %13.6e %4d\n",x,intensity(tmin,&h)*yfac,h.niso);
        }
      }
      else
      {
        if(yuni == YUNI_LOG)
        {
          for(t=tmin; t<=tmax; t+=tstp)
          {
            printf("%17.8f %7.2f %13.6e %4d\n",x,t,log10(intensity(t,&h)*yfac),h.niso);
          }
        }
        else
        {
          for(t=tmin; t<=tmax; t+=tstp)
          {
            printf("%17.8f %7.2f %13.6e %4d\n",x,t,intensity(t,&h)*yfac,h.niso);
          }
        }
      }
    }
    else
    {
      if(tmax < tmin)
      {
        if(yuni == YUNI_LOG)
        {
          printf("%17.8f %13.6e\n",x,log10(intensity(tmin,&h)*yfac));
        }
        else
        {
          printf("%17.8f %13.6e\n",x,intensity(tmin,&h)*yfac);
        }
      }
      else
      {
        if(yuni == YUNI_LOG)
        {
          for(t=tmin; t<=tmax; t+=tstp)
          {
            printf("%17.8f %7.2f %13.6e\n",x,t,log10(intensity(t,&h)*yfac));
          }
        }
        else
        {
          for(t=tmin; t<=tmax; t+=tstp)
          {
            printf("%17.8f %7.2f %13.6e\n",x,t,intensity(t,&h)*yfac);
          }
        }
      }
    }
  }

  return 0;
}

int Init(void)
{
  switch(xuni)
  {
    case XUNI_HZ:
      xfac = C_LIGHT*1.0e2;
      break;
    case XUNI_GHZ:
      xfac = C_LIGHT*1.0e-7;
      break;
    default:
      xfac = 1.0;
      break;
  }
  switch(yuni)
  {
    case YUNI_J:
    case YUNI_LOG:
      yfac = YFAC_H2J;
      break;
    default:
      yfac = 1.0;
      break;
  }
  yfac /= ynrm;

  if(hver >= 2012)
  {
    tips = tips12;
  } else
  if(hver >= 2008)
  {
    tips = tips08;
  } else
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

  return 0;
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
    {"nmol",1,0,'m'},
    {"niso",1,0,'M'},
    {"riso",1,0,'r'},
    {"xmin",1,0,'x'},
    {"xmax",1,0,'X'},
    {"xuni",1,0,'U'},
    {"yuni",1,0,'V'},
    {"ynrm",1,0,'W'},
    {"tmin",1,0,'t'},
    {"tmax",1,0,'T'},
    {"tstp",1,0,'S'},
    {"hver",1,0,'n'},
    {"long",0,0,'l'},
    {"debug",0,0,'d'},
    {"verbose",0,0,'v'},
    {"help",0,0,'h'}
  };

  rt = 0;
  while(1)
  {
    this_option_optind = optind?optind:1;
    c = getopt_long(argn,args,":f:m:M:r:x:X:U:V:W:t:T:S:n:ldvh",long_options,&option_index);
    if(c == -1) break;

    switch(c)
    {
      case 'f':
        strncpy(dnam,optarg,MAXCHAR);
        break;
      case 'm':
        errno = 0;
        ntmp = strtol(optarg,&p,10);
        if(errno!=ERANGE && *p=='\0' && ntmp>0 && ntmp<100) nmol = ntmp;
        else
        {
          fprintf(stderr,"Mol# -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'M':
        errno = 0;
        ntmp = strtol(optarg,&p,10);
        if(errno!=ERANGE && *p=='\0' && ntmp>0 && ntmp<100) niso = ntmp;
        else
        {
          fprintf(stderr,"Iso# -> out of range %s\n",optarg);
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
      case 'V':
        errno = 0;
        ntmp = strtol(optarg,&p,10);
        if(errno!=ERANGE && *p=='\0' && ntmp>=1 && ntmp<=3) yuni = ntmp;
        else
        {
          fprintf(stderr,"Y unit -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'W':
        errno = 0;
        xtmp = strtod(optarg,&p);
        if(errno!=ERANGE && *p=='\0' && xtmp>0.0) ynrm = xtmp;
        else
        {
          fprintf(stderr,"Y weight -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 't':
        errno = 0;
        xtmp = strtod(optarg,&p);
        if(errno!=ERANGE && *p=='\0' && xtmp>0.0) tmin = xtmp;
        else
        {
          fprintf(stderr,"T min -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'T':
        errno = 0;
        xtmp = strtod(optarg,&p);
        if(errno!=ERANGE && *p=='\0' && xtmp>0.0) tmax = xtmp;
        else
        {
          fprintf(stderr,"T max -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'S':
        errno = 0;
        xtmp = strtod(optarg,&p);
        if(errno!=ERANGE && *p=='\0' && xtmp>0.0) tstp = xtmp;
        else
        {
          fprintf(stderr,"T step -> out of range %s\n",optarg);
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
      case 'l':
        lm++;
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
  fprintf(stderr,"hitran_intensity -[option] (argument) -[option] (argument) ...\n");
  fprintf(stderr,"-----------------------------------------------------------------------------\n");
  fprintf(stderr,"   option   |%s|%s|%s| current\n",As(e,"",n),         As(a,"argument",n),   As(d,"default",n));
  fprintf(stderr," f -dnam    |%s|%s|%s| %s\n",As(e,"Data file",n),     As(a,"name",n),       As(d,DNAM,n),dnam);
  fprintf(stderr," m -nmol    |%s|%s|%s| %d\n",As(e,"Mol#",n),          As(a,"#",n),          Ad(d,IMOL,n),nmol);
  fprintf(stderr," M -niso    |%s|%s|%s| %d\n",As(e,"Iso#",n),          As(a,"#",n),          Ad(d,NISO,n),niso);
  fprintf(stderr," r -riso    |%s|%s|%s| %e\n",As(e,"Abundance ratio ",n),As(a,"mol/mol",n),  Ae(d,RISO,n),riso);
  fprintf(stderr," x -xmin    |%s|%s|%s| %e\n",As(e,"X Min",n),         As(a,"X unit",n),     Af(d,XMIN,n),xmin);
  fprintf(stderr," X -xmax    |%s|%s|%s| %e\n",As(e,"X Max",n),         As(a,"X unit",n),     Af(d,XMAX,n),xmax);
  fprintf(stderr," U -xuni    |%s|%s|%s| %d\n",As(e,"X Unit",n),        As(a,"Unit#",n),      Ad(d,XUNI,n),xuni);
  fprintf(stderr," V -yuni    |%s|%s|%s| %d\n",As(e,"Y Unit",n),        As(a,"Unit#",n),      Ad(d,YUNI,n),yuni);
  fprintf(stderr," W -ynrm    |%s|%s|%s| %e\n",As(e,"Y Weight",n),      As(a,"Value",n),      Ad(d,YNRM,n),ynrm);
  fprintf(stderr," t -tmin    |%s|%s|%s| %f\n",As(e,"T Min",n),         As(a,"K",n),          Af(d,TMIN,n),tmin);
  fprintf(stderr," T -tmax    |%s|%s|%s| %f\n",As(e,"T Max",n),         As(a,"K",n),          Af(d,TMAX,n),tmax);
  fprintf(stderr," S -tstp    |%s|%s|%s| %f\n",As(e,"T Step",n),        As(a,"K",n),          Af(d,TSTP,n),tstp);
  fprintf(stderr," n -hver    |%s|%s|%s| %d\n",As(e,"HITRAN version",n),As(a,"version",n),    Ad(d,HVER,n),hver);
  fprintf(stderr," l -long    |%s|%s|%s| %d\n",As(e,"Long    mode",n),  As(a,"nothing",n),    Ad(d,0,n),lm);
  fprintf(stderr," d -debug   |%s|%s|%s| %d\n",As(e,"Debug   mode",n),  As(a,"nothing",n),    Ad(d,0,n),db);
  fprintf(stderr," v -verbose |%s|%s|%s| %d\n",As(e,"Verbose mode",n),  As(a,"nothing",n),    Ad(d,0,n),vb);
  fprintf(stderr," h -help    |%s|%s|%s| %d\n",As(e,"Help    mode",n),  As(a,"nothing",n),    Ad(d,0,n),1);
  fprintf(stderr,"-----------------------------------------------------------------------------\n");
  fprintf(stderr,"X Unit 1:nm 2:cm-1 3:Hz 4:GHz\n");
  fprintf(stderr,"Y Unit 1:HITRAN 2:JPL 3:JPL(log10)\n");

  return 0;
}
