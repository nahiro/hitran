#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <getopt.h>
#include <math.h>
#include <bits/nan.h>
#include "strutil.h"
#include "hitran_common.c"

#define		IMOL			1
#define		IISO			1
#define		TOUT			300.0
#define		HVER			2012			// HITRAN version

int		imol			= IMOL;
int		iiso			= IISO;
double		tout			= TOUT;
int		hver			= HVER;				// HITRAN04
int		vb			= 0;				// Verbose mode
int		db			= 0;				// Debug   mode
int		hp			= 0;				// Help    mode

int Init(void);
int GetOpt(int argn,char **args);
int Usage(void);

int main(int argc,char **argv)
{
  if(GetOpt(argc,argv) < 0) exit(-1);
  if(hp) {Usage(); return 0;}
  if(Init() < 0) exit(-1);

  printf("%2d %2d %10.3f %13.6f\n",imol,iiso,tout,tips(imol,iiso,tout));

  return 0;
}

int Init(void)
{
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
    {"dnam",1,0,'D'},
    {"imol",1,0,'m'},
    {"iiso",1,0,'i'},
    {"tout",1,0,'o'},
    {"tref",1,0,'r'},
    {"hver",1,0,'n'},
    {"debug",0,0,'d'},
    {"verbose",0,0,'v'},
    {"help",0,0,'h'}
  };

  rt = 0;
  while(1)
  {
    this_option_optind = optind?optind:1;
    c = getopt_long(argn,args,":D:m:i:o:r:n:dvh",long_options,&option_index);
    if(c == -1) break;

    switch(c)
    {
      case 'D':
        strncpy(dnam,optarg,MAXCHAR);
        break;
      case 'm':
        errno = 0;
        ntmp = strtol(optarg,&p,10);
        if(errno!=ERANGE && *p=='\0' && ntmp>0) imol = ntmp;
        else
        {
          fprintf(stderr,"Molecule# -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'i':
        errno = 0;
        ntmp = strtol(optarg,&p,10);
        if(errno!=ERANGE && *p=='\0' && ntmp>0) iiso = ntmp;
        else
        {
          fprintf(stderr,"Isotope# -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'o':
        errno = 0;
        xtmp = strtod(optarg,&p);
        if(errno!=ERANGE && *p=='\0' && xtmp>=70.0 && xtmp<=3005.0) tout = xtmp;
        else
        {
          fprintf(stderr,"Output temperature -> out of range %s\n",optarg);
          rt = -1;
        }
        break;
      case 'r':
        errno = 0;
        xtmp = strtod(optarg,&p);
        if(errno!=ERANGE && *p=='\0' && xtmp>=70.0 && xtmp<=3005.0) tref = xtmp;
        else
        {
          fprintf(stderr,"Reference temperature -> out of range %s\n",optarg);
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
  fprintf(stderr,"hitran_tips -[option] (argument) -[option] (argument) ...\n");
  fprintf(stderr,"-----------------------------------------------------------------------------\n");
  fprintf(stderr,"   option   |%s|%s|%s| current\n",As(e,"",n),         As(a,"argument",n),   As(d,"default",n));
  fprintf(stderr," D -dnam    |%s|%s|%s| %s\n",As(e,"Data file",n),     As(a,"name",n),       As(d,DNAM,n),dnam);
  fprintf(stderr," m -imol    |%s|%s|%s| %d\n",As(e,"Mol#",n),          As(a,"#",n),          Ad(d,IMOL,n),imol);
  fprintf(stderr," i -iiso    |%s|%s|%s| %d\n",As(e,"Iso#",n),          As(a,"#",n),          Ad(d,IISO,n),iiso);
  fprintf(stderr," o -tout    |%s|%s|%s| %e\n",As(e,"Out. temp",n),     As(a,"K",n),          Af(d,TOUT,n),tout);
  fprintf(stderr," r -tref    |%s|%s|%s| %f\n",As(e,"Ref. temp",n),     As(a,"K",n),          Af(d,TREF,n),tref);
  fprintf(stderr," n -hver    |%s|%s|%s| %d\n",As(e,"HITRAN version",n),As(a,"version",n),    Ad(d,HVER,n),hver);
  fprintf(stderr," d -debug   |%s|%s|%s| %d\n",As(e,"Debug   mode",n),  As(a,"nothing",n),    Ad(d,0,n),db);
  fprintf(stderr," v -verbose |%s|%s|%s| %d\n",As(e,"Verbose mode",n),  As(a,"nothing",n),    Ad(d,0,n),vb);
  fprintf(stderr," h -help    |%s|%s|%s| %d\n",As(e,"Help    mode",n),  As(a,"nothing",n),    Ad(d,0,n),1);
  fprintf(stderr,"-----------------------------------------------------------------------------\n");

  return 0;
}
