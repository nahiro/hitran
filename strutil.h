#ifdef __cplusplus
extern "C" {
#endif

#define MAX_STR_LEN		256

char *As(char *dst,char *src,int size);
char *Ad(char *dst,int d,int size);
char *Af(char *dst,double f,int size);
char *Ae(char *dst,double e,int size);
char *Ac(char *dst,char c,int size);
int ungets(char *s, int size, FILE *stream);
int getarg(char *src,int size,char *dst,char *ptr[]);

#ifdef __cplusplus
}
#endif
