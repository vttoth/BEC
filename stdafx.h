#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex>
#include <math.h>

#define _TCHAR char
#define _tmain main
#define GetShortPathNameA(src,dst,len) strncpy(dst,src,len)
#define _getcwd getcwd
#define MAX_PATH PATH_MAX
