/*

  utils.c

  some handy routines


  (c) 2003 Thomas Neff

*/

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <sys/times.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "md5.h"
#include "utils.h"


// why is this not defined in stdio.h as the man page says ?
char* cuserid(char* string);


static char infostring[4096];

void createinfo(int argc, char* argv[])
{
  int i; 
  char *c = infostring;

  sprintf(c, "# FMD git: %s \n#\n", git_version()); 
  sprintf(c+strlen(infostring), "# %s @ %s, %s\n", cuserid(NULL), hostname(), curdate());
  sprintf(c+strlen(infostring), "# ");
  for (i=0; i<argc; i++)
    sprintf(c+strlen(infostring), "%s ", argv[i]);
  sprintf(c+strlen(infostring), "\n");
}


void fprintinfo(FILE* fp)
{
  fputs(infostring, fp);
}

void gzprintinfo(gzFile fp)
{
  gzputs(fp, infostring);
}


char* hostname(void)
{
  char* host = (char*) malloc(255*sizeof(char));

  gethostname(host, 255);
  return host;
}


char* curdate(void)
{
  time_t curtime;
  struct tm *curtimep;

  curtime = time(NULL);
  curtimep = localtime(&curtime);

  return (stripstr(asctime(curtimep),"\n"));
}


int usertime(void)
{
  struct tms t;

  times(&t);
  return (t.tms_utime / sysconf(_SC_CLK_TCK));
}


char* stripstr(char* str, const char* strip)
{
  char* cp;

  if ((cp = strstr(str, strip)))
    *cp = '\0';

  return str;
}


char* strjoin(const char* stra, const char* strb)
{
  char* c = malloc(strlen(stra)+strlen(strb));
  sprintf(c, "%s%s", stra, strb);
  return c;
}


int fileexists(const char* fname)
{
  struct stat filestat;

  if (!stat(fname, &filestat))
    if (S_ISREG(filestat.st_mode))
      return 0;

  return -1;
}


void ensuredir(const char* dirname)
{
  struct stat dirstat;

  if (!stat(dirname, &dirstat)) {
    if (!S_ISDIR(dirstat.st_mode)) {
      fprintf(stderr, "%s is not a directory\n", dirname);
      exit(-1);
    }
  } else {
    if (mkdir(dirname, 0755)) {
      fprintf(stderr, "couldn't create directory %s\n", dirname);
      exit(-1);
    }
  }
}
	      

int backup(const char* fname)
{
  char backupname[1024];

  if (!access(fname, F_OK)) {
    ensuredir("BAK");
    sprintf(backupname, "BAK/%s-%d", fname, (int) time(NULL));
    return rename(fname, backupname);
  }
  return 0;
}


char* filepart(const char* fullname)
{
  char *p;
  char *filep;

  if ((p = strrchr(fullname, '/'))) {
    filep = (char*) malloc(strlen(p)*sizeof(char));
    strcpy(filep, p+1);
  } else {
    filep = (char*) malloc((strlen(fullname)+1)*sizeof(char));
    strcpy(filep, fullname);
  }
  return filep;
}

char* pathpart(const char* fullname)
{
  char *p, *pp;
  char* pathp;
  const char *i;

  if ((p = strrchr(fullname, '/'))) {
    pathp = (char*) malloc((strlen(fullname)-strlen(p)+2)*sizeof(char));
    for (pp = pathp, i = fullname; i<=p; i++)
      *pp++ = *i;
    *pp = '\0';
  } else {
    pathp = (char*) malloc(sizeof(char));
    *pathp = '\0';
  }
  return pathp;
}

#define BUFSIZE 1024
char* md5hash(const char* fname)
{
  md5_state_t state;
  md5_byte_t digest[16];
  char* hex_output = malloc(33*sizeof(char));
  int di;
  FILE* fp;
  char buf[BUFSIZE];

  if (!(fp = fopen(fname, "r"))) {
    fprintf(stderr, "couldn't open %s for reading\n", fname);
    exit(-1);
  }

  md5_init(&state);
  while(1) {
    fgets(buf, BUFSIZE, fp);
    if (feof(fp)) break;
    md5_append(&state, (const md5_byte_t *) buf, strlen(buf));
  }
  md5_finish(&state, digest);
  fclose(fp);
  
  for (di=0; di<16; di++)
    sprintf(hex_output+di*2, "%02x", digest[di]);

  return hex_output;
}
  

char* nucleusIDLformat(const char* name)
{
  int i=0;
  char* idlstr;

  idlstr = (char*) malloc(9*sizeof(char));
  idlstr[i++] = '!'; idlstr[i++] = 'U';
  if (isdigit(name[1]))
    idlstr[i++] = name[1];
  if (isdigit(name[2]))
    idlstr[i++] = name[2];
  if (isdigit(name[3]))
    idlstr[i++] = name[3];
  idlstr[i++] = '!'; idlstr[i++] = 'N';
  if (isalpha(name[0])) 
    idlstr[i++] = name[0];
  if (isalpha(name[1])) 
    idlstr[i++] = name[1];
  idlstr[i] = '\0';

  return idlstr; 
}


#define STRLEN 255
int readstringsfromfile(const char* fname, int* n, char* string[])
{
  FILE* fp;

  if (!(fp = fopen(fname, "r"))) {
    fprintf(stderr, "couldn't open %s for reading\n", fname);
    exit(-1);
  }

  int i=0; 
  int res;
  do {
    string[i] = malloc(STRLEN*sizeof(char));
    res = fscanf(fp, " %s", string[i]);
    i++;
  } while (res != EOF);

  fclose(fp);

  *n = i-1;

  return 0;
}
