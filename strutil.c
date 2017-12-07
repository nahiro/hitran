#include <stdio.h>
#include <string.h>
#include "strutil.h"

char *As(char *dst,char *src,int size)
{
  int i,n,p,q;

  if (size > MAX_STR_LEN) {
    fprintf(stderr,"Align center: size too big %d\n",size);
    return (char *)NULL;
  }

  n = strlen(src);
  if (n > size) {
    p = 0;
    q = size;
  } else {
    p = (size-n)/2;
    q = p+n;
  }

  for (i=0; i<p; i++) {
    *(dst+i) = ' ';
  }
  for (i=p; i<q; i++) {
    *(dst+i) = *(src+i-p);
  }
  for (i=q; i<size; i++) {
    *(dst+i) = ' ';
  }
  *(dst+size) = '\0';

  return dst;
}

char *Ad(char *dst,int d,int size)
{
  int i,n,p,q;
  char src[MAX_STR_LEN];

  if (size > MAX_STR_LEN) {
    fprintf(stderr,"Align center: size too big %d\n",size);
    return (char *)NULL;
  }

  sprintf(src,"%d",d);
  n = strlen(src);
  if (n > size) {
    p = 0;
    q = size;
  } else {
    p = (size-n)/2;
    q = p+n;
  }

  for (i=0; i<p; i++) {
    *(dst+i) = ' ';
  }
  for (i=p; i<q; i++) {
    *(dst+i) = *(src+i-p);
  }
  for (i=q; i<size; i++) {
    *(dst+i) = ' ';
  }
  *(dst+size) = '\0';

  return dst;
}

char *Af(char *dst,double f,int size)
{
  int i,n,p,q;
  char src[MAX_STR_LEN];

  if (size > MAX_STR_LEN) {
    fprintf(stderr,"Align center: size too big %d\n",size);
    return (char *)NULL;
  }

  sprintf(src,"%.2f",f);
  n = strlen(src);
  if (n > size) {
    p = 0;
    q = size;
  } else {
    p = (size-n)/2;
    q = p+n;
  }

  for (i=0; i<p; i++) {
    *(dst+i) = ' ';
  }
  for (i=p; i<q; i++) {
    *(dst+i) = *(src+i-p);
  }
  for (i=q; i<size; i++) {
    *(dst+i) = ' ';
  }
  *(dst+size) = '\0';

  return dst;
}

char *Ae(char *dst,double e,int size)
{
  int i,n,p,q;
  char src[MAX_STR_LEN];

  if (size > MAX_STR_LEN) {
    fprintf(stderr,"Align center: size too big %d\n",size);
    return (char *)NULL;
  }

  sprintf(src,"%.2e",e);
  n = strlen(src);
  if (n > size) {
    p = 0;
    q = size;
  } else {
    p = (size-n)/2;
    q = p+n;
  }

  for (i=0; i<p; i++) {
    *(dst+i) = ' ';
  }
  for (i=p; i<q; i++) {
    *(dst+i) = *(src+i-p);
  }
  for (i=q; i<size; i++) {
    *(dst+i) = ' ';
  }
  *(dst+size) = '\0';

  return dst;
}

char *Ac(char *dst,char c,int size)
{
  int i,n,p,q;
  char src[MAX_STR_LEN];

  if (size > MAX_STR_LEN) {
    fprintf(stderr,"Align center: size too big %d\n",size);
    return (char *)NULL;
  }

  sprintf(src,"%c",c);
  n = strlen(src);
  if (n > size) {
    p = 0;
    q = size;
  } else {
    p = (size-n)/2;
    q = p+n;
  }

  for (i=0; i<p; i++) {
    *(dst+i) = ' ';
  }
  for (i=p; i<q; i++) {
    *(dst+i) = *(src+i-p);
  }
  for (i=q; i<size; i++) {
    *(dst+i) = ' ';
  }
  *(dst+size) = '\0';

  return dst;
}

int ungets(char *s, int size, FILE *stream)
{
  int i,end;
  char *c;

  end = strlen(s)-1;
  c = s+end;

  for (i=end; i>=0; i--,c--)
  {
    if (ungetc((int)*c,stream) != (int)*c)
    {
      fprintf(stderr,"Error in ungets >>> %s\n",s);
      return -1;
    }
  }

  return 0;
}

int getarg(char *src,int size,char *dst,char *ptr[])
{
  int i,j,k;
  int flg;
  char c;

  flg = 0;
  for(i=j=k=0; i<size; i++)
  {
    if(i == size-1)
    {
      *(dst+j) = '\0';
      break;
    }
    c = *(src+i);
    switch(c)
    {
      case '\n':
      case '\0':
        *(dst+j) = '\0';
        i = size;
        break;
      case ' ':
      case '\t':
        if(flg == 1)
        {
          *(dst+j) = c;
          if(j==0 || (j>0 &&(*(dst+j-1)=='\0')))
          {
            ptr[k++] = dst+j;
          }
          j++;
        }
        else
        {
          *(dst+j) = '\0';
          j++;
        }
        break;
      case '\\':
        if(i<size-1 && *(src+i+1)!='"')
        {
          *(dst+j) = c;
          if(j==0 || (j>0 &&(*(dst+j-1)=='\0')))
          {
            ptr[k++] = dst+j;
          }
          j++;
        }
        break;
      case '"':
        if(i>0 && *(src+i-1)=='\\')
        {
          *(dst+j) = c;
          if(j==0 || (j>0 &&(*(dst+j-1)=='\0')))
          {
            ptr[k++] = dst+j;
          }
          j++;
        }
        else
        {
          flg = 1-flg;
        }
        break;
      default:
        *(dst+j) = c;
        if(j==0 || (j>0 &&(*(dst+j-1)=='\0')))
        {
          ptr[k++] = dst+j;
        }
        j++;
        break;
    }
  }

  return k;
}
