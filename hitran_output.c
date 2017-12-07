/***********************************************************/
/* HITRAN_OUTPUT     ... Output contents of HITRAN         */
/*                       database.                         */
/* Author: N.Manago                                        */
/* $Revision: 39 $                                         */
/* $Date: 2009-10-03 21:36:26 +0900 (Sat, 03 Oct 2009) $   */
/***********************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <getopt.h>
#include <math.h>
#include <bits/nan.h>
#include "strutil.h"
#include "numutil.h"
#ifndef	NAN
#define	NAN				(nanval())
#endif
#include "hitran_common.c"

#define		NISO			-1			// Isotope ID
#define		XUNI_NM			1			// Wavelength in nm
#define		XUNI_WN			2			// Wavenumber in cm-1
#define		XUNI_HZ			3			// Frequency in Hz
#define		XUNI_GHZ		4			// Frequency in GHz
#define		XUNI			XUNI_NM			// X Unit
#define		YUNI_H			1			// HITRAN unit
#define		YUNI_J			2			// JPL unit
#define		YUNI			YUNI_H			// Y Unit
#define		YFAC_H2J		2.99792458e18		// HITRAN->JPL conversion factor
#define		XMIN			350.0			// X Min
#define		XMAX			1050.0			// X Max
#define		HVER			2012			// HITRAN version
#define		SPAC			'_'			// Space character

double		xmin			= XMIN;
double		xmax			= XMAX;
double		xfac			= 1.0;			// X Factor
double		yfac			= 1.0;			// Y Factor
int		niso			= NISO;			// Isotope ID
int		xuni			= XUNI;			// X Unit
int		yuni			= YUNI;			// Y Unit
char		spac			= SPAC;			// Space character
int		hver			= HVER;			// HITRAN version
int		vb			= 0;			// Verbose mode
int		db			= 0;			// Debug   mode
int		hp			= 0;			// Help    mode

int Init(void);
int GetOpt(int argn,char **args);
int Usage(void);

int main(int argc,char **argv)
{
  double x,y;
  char line[MAXCHAR];
  char *p,*q;
  struct hitran96 h96;
  struct hitran04 h04;

  if(GetOpt(argc,argv) < 0) return -1;
  if(hp) {Usage(); return 0;}
  if(Init() < 0) return -1;

  while(fgets(line,MAXCHAR,stdin) != NULL)
  {
    if(hver >= 2004)
    {
      if(get_hitran04(line,&h04) < 0) continue;
      if(niso > 0)
      {
        if(h04.niso != niso) continue;
      }
      x = h04.freq*xfac;
      y = h04.sint*yfac;
    }
    else
    {
      if(get_hitran96(line,&h96) < 0) continue;
      if(niso > 0)
      {
        if(h96.niso != niso) continue;
      }
      x = h96.freq*xfac;
      y = h96.sint*yfac;
    }
    if(xuni == XUNI_NM)
    {
      x = 1.0e7/x;
    }
    if(x>xmin && x<xmax)
    {
      if(hver >= 2004)
      {
        // strip spaces
        q = h04.gqup;
        for(p=strchr(q,'\0')-1; *p==' '&&p>=q; p--) *p = '\0';
        for(p=q; *p==' '; p++);
        while(1)
        {
          *q = (*p==' '?spac:*p);
          if(!*p) break;
          p++;
          q++;
        }
        if(strlen(h04.gqup) == 0) sprintf(h04.gqup,"%c",spac);
        q = h04.gqlo;
        for(p=strchr(q,'\0')-1; *p==' '&&p>=q; p--) *p = '\0';
        for(p=q; *p==' '; p++);
        while(1)
        {
          *q = (*p==' '?spac:*p);
          if(!*p) break;
          p++;
          q++;
        }
        if(strlen(h04.gqlo) == 0) sprintf(h04.gqlo,"%c",spac);
        q = h04.lqup;
        for(p=strchr(q,'\0')-1; *p==' '&&p>=q; p--) *p = '\0';
        for(p=q; *p==' '; p++);
        while(1)
        {
          *q = (*p==' '?spac:*p);
          if(!*p) break;
          p++;
          q++;
        }
        if(strlen(h04.lqup) == 0) sprintf(h04.lqup,"%c",spac);
        q = h04.lqlo;
        for(p=strchr(q,'\0')-1; *p==' '&&p>=q; p--) *p = '\0';
        for(p=q; *p==' '; p++);
        while(1)
        {
          *q = (*p==' '?spac:*p);
          if(!*p) break;
          p++;
          q++;
        }
        if(strlen(h04.lqlo) == 0) sprintf(h04.lqlo,"%c",spac);
        q = h04.flag;
        for(p=strchr(q,'\0')-1; *p==' '&&p>=q; p--) *p = '\0';
        for(p=q; *p==' '; p++);
        while(1)
        {
          *q = (*p==' '?spac:*p);
          if(!*p) break;
          p++;
          q++;
        }
        if(strlen(h04.flag) == 0) sprintf(h04.flag,"%c",spac);
        printf("%2d %d %15.8e %13.6e %10.3e %5.4f %5.4f "
               "%10.4f %4.2f %8.6f %15s %15s %15s %15s %d %d %d %d %d %d %2d %2d %2d %2d %2d %2d %s %7.1f %7.1f\n",
                h04.nmol,h04.niso,x,y,h04.aein,h04.wair,h04.wslf,h04.elow,h04.ctmp,h04.psft,
                h04.gqup,h04.gqlo,h04.lqup,h04.lqlo,h04.err1,h04.err2,h04.err3,h04.err4,h04.err5,h04.err6,
                h04.ref1,h04.ref2,h04.ref3,h04.ref4,h04.ref5,h04.ref6,h04.flag,h04.swup,h04.swlo);
      }
      else
      {
        q = h96.lqup;
        for(p=strchr(q,'\0')-1; *p==' '&&p>=q; p--) *p = '\0';
        for(p=q; *p==' '; p++);
        while(1)
        {
          *q = (*p==' '?spac:*p);
          if(!*p) break;
          p++;
          q++;
        }
        if(strlen(h96.lqup) == 0) sprintf(h96.lqup,"%c",spac);
        q = h96.lqlo;
        for(p=strchr(q,'\0')-1; *p==' '&&p>=q; p--) *p = '\0';
        for(p=q; *p==' '; p++);
        while(1)
        {
          *q = (*p==' '?spac:*p);
          if(!*p) break;
          p++;
          q++;
        }
        if(strlen(h96.lqlo) == 0) sprintf(h96.lqlo,"%c",spac);
        printf("%2d %d %13.6e %13.6e %10.3e %5.4f %5.4f "
               "%10.4f %4.2f %8.6f %3d %3d %9s %9s %d %d %d %2d %2d %2d\n",
                h96.nmol,h96.niso,x,y,h96.mtrn,h96.wair,h96.wslf,h96.elow,h96.ctmp,h96.psft,
                h96.gqup,h96.gqlo,h96.lqup,h96.lqlo,h96.err1,h96.err2,h96.err3,h96.ref1,h96.ref2,h96.ref3);
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
      yfac = YFAC_H2J;
      break;
    default:
      yfac = 1.0;
      break;
  }

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
    {"niso",1,0,'N'},
    {"xmin",1,0,'x'},
    {"xmax",1,0,'X'},
    {"xuni",1,0,'U'},
    {"yuni",1,0,'V'},
    {"space",1,0,'S'},
    {"hver",1,0,'n'},
    {"debug",0,0,'d'},
    {"verbose",0,0,'v'},
    {"help",0,0,'h'}
  };

  rt = 0;
  while(1)
  {
    this_option_optind = optind?optind:1;
    c = getopt_long(argn,args,":f:N:x:X:U:V:S:n:dvh",long_options,&option_index);
    if(c == -1) break;

    switch(c)
    {
      case 'f':
        strncpy(dnam,optarg,MAXCHAR);
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
        if(errno!=ERANGE && *p=='\0' && ntmp>=1 && ntmp<=2) yuni = ntmp;
        else
        {
          fprintf(stderr,"Y unit -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'S':
        spac = optarg[0];
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
  fprintf(stderr,"hitran_output -[option] (argument) -[option] (argument) ...\n");
  fprintf(stderr,"-----------------------------------------------------------------------------\n");
  fprintf(stderr,"   option   |%s|%s|%s| current\n",As(e,"",n),         As(a,"argument",n),   As(d,"default",n));
  fprintf(stderr," f -dnam    |%s|%s|%s| %s\n",As(e,"Data file",n),     As(a,"name",n),       As(d,DNAM,n),dnam);
  fprintf(stderr," N -niso    |%s|%s|%s| %d\n",As(e,"Isotope ID",n),    As(a,"#",n),          Ad(d,NISO,n),niso);
  fprintf(stderr," x -xmin    |%s|%s|%s| %f\n",As(e,"X Min",n),         As(a,"X unit",n),     Af(d,XMIN,n),xmin);
  fprintf(stderr," X -xmax    |%s|%s|%s| %f\n",As(e,"X Max",n),         As(a,"X unit",n),     Af(d,XMAX,n),xmax);
  fprintf(stderr," U -xuni    |%s|%s|%s| %d\n",As(e,"X Unit",n),        As(a,"Unit#",n),      Ad(d,XUNI,n),xuni);
  fprintf(stderr," V -yuni    |%s|%s|%s| %d\n",As(e,"Y Unit",n),        As(a,"Unit#",n),      Ad(d,YUNI,n),yuni);
  fprintf(stderr," S -space   |%s|%s|%s| %c\n",As(e,"Space caracter",n),As(a,"character",n),  Ac(d,SPAC,n),spac);
  fprintf(stderr," n -hver    |%s|%s|%s| %d\n",As(e,"HITRAN version",n),As(a,"version",n),    Ad(d,HVER,n),hver);
  fprintf(stderr," d -debug   |%s|%s|%s| %d\n",As(e,"Debug   mode",n),  As(a,"nothing",n),    Ad(d,0,n),db);
  fprintf(stderr," v -verbose |%s|%s|%s| %d\n",As(e,"Verbose mode",n),  As(a,"nothing",n),    Ad(d,0,n),vb);
  fprintf(stderr," h -help    |%s|%s|%s| %d\n",As(e,"Help    mode",n),  As(a,"nothing",n),    Ad(d,0,n),1);
  fprintf(stderr,"-----------------------------------------------------------------------------\n");
  fprintf(stderr,"X Unit 1:nm 2:cm-1 3:Hz 4:GHz\n");
  fprintf(stderr,"Y Unit 1:HITRAN 2:JPL\n");

  return 0;
}
