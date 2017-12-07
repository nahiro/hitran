#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <bits/nan.h>

#define		HITRAN_LINE_LENGTH		100
#define		HITRAN_NITEM			17
#define		MAXCHAR				256

int get_chr(char *str,int start,int size,char *value);
int get_int(char *str,int start,int size,int *value);
int get_dbl(char *str,int start,int size,double *value);

int main(void)
{
  int n;
  int err;
  char line[MAXCHAR];
  char *p;
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

  while(fgets(line,MAXCHAR,stdin) != NULL)
  {
    if((p=strchr(line,'\r')) != NULL)
    {
      *p = '\0';
    }
    if((p=strchr(line,'\n')) != NULL)
    {
      *p = '\0';
    }
    if(strlen(line) != HITRAN_LINE_LENGTH)
    {
      fprintf(stderr,"Error, length=%ld >>> %s\n",strlen(line),line);
      return -1;
    }
    // read values
    err = 0;
    n = 0;
    if(get_int(line,n, 2,&nmol) < 0) err = 1; n +=  2;
    if(get_int(line,n, 1,&niso) < 0) err = 1; n +=  1;
    if(get_dbl(line,n,12,&freq) < 0) err = 1; n += 12;
    if(get_dbl(line,n,10,&sint) < 0) err = 1; n += 10;
    if(get_dbl(line,n,10,&mtrn) < 0) err = 1; n += 10;
    if(get_dbl(line,n, 5,&wair) < 0) err = 1; n +=  5;
    if(get_dbl(line,n, 5,&wslf) < 0) err = 1; n +=  5;
    if(get_dbl(line,n,10,&elow) < 0) err = 1; n += 10;
    if(get_dbl(line,n, 4,&ctmp) < 0) err = 1; n +=  4;
    if(get_dbl(line,n, 8,&psft) < 0) err = 1; n +=  8;
    if(get_int(line,n, 3,&gqup) < 0) err = 1; n +=  3;
    if(get_int(line,n, 3,&gqlo) < 0) err = 1; n +=  3;
    if(get_chr(line,n, 9,lqup)  < 0) err = 1; n +=  9;
    if(get_chr(line,n, 9,lqlo)  < 0) err = 1; n +=  9;
    if(get_int(line,n, 1,&err1) < 0) err = 1; n +=  1;
    if(get_int(line,n, 1,&err2) < 0) err = 1; n +=  1;
    if(get_int(line,n, 1,&err3) < 0) err = 1; n +=  1;
    if(get_int(line,n, 2,&ref1) < 0) err = 1; n +=  2;
    if(get_int(line,n, 2,&ref2) < 0) err = 1; n +=  2;
    if(get_int(line,n, 2,&ref3) < 0) err = 1; n +=  2;
    if(err)
    {
      fprintf(stderr,"Read error >>> %s\n",line);
      return -1;
    }
    printf("%2d %d %13.6f %10.3e %10.3e %6.4f %6.4f %11.4f %5.2f %9.6f %3d %3d %9s %9s %d %d %d %2d %2d %2d\n",
            nmol,niso,freq,sint,mtrn,wair,wslf,elow,ctmp,psft,
            gqup,gqlo,lqup,lqlo,err1,err2,err3,ref1,ref2,ref3);
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
