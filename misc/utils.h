/*

  utils.h

  some handy routines


  (c) 2003 Thomas Neff

*/


#ifndef _UTILS_H
#define _UTILS_H

#include <stdio.h>
#include <zlib.h>

inline static int min(int a, int b) { return(a<b ? a : b); }
inline static int max(int a, int b) { return(a>b ? a : b); }


/// get current SVN version
const char* svn_version(void);

/// save info about running process
void createinfo(int argc, char* argv[]);

/// print info about running process
void fprintinfo(FILE* fp);

/// print info about running process
void gzprintinfo(gzFile fp);

/// returns string with hostname
char* hostname(void);

/// returns string with current time and date
char* curdate(void);

/// returns process user time in seconds
int usertime(void);

/// strips strip from str
char* stripstr(char* str, const char* strip);

/// join two strings
char* strjoin(const char* stra, const char* strb);

/// check if file exists and is readable (returns 0 for success)
int fileexists(const char* fname);

/// check for existence, create if not of directory
void ensuredir(const char* dirname);

/// renames fname to fname.bak
int backup(const char* fname);

/// get the file component of fullname
char* filepart(const char* fullname);

/// get the path component of fullname
char* pathpart(const char* fullname);

/// create md5 hash for file fname
char* md5hash(const char* fname);

/// convert name e.g. "Si28" into IDL format 
char* nucleusIDLformat(const char* name);

/// read list of strings from file
int readstringsfromfile(const char* fname, int* n, char* string[]);


#endif
