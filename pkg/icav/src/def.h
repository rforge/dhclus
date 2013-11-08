#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Boolean.h>

#include <float.h>
#include <string.h>
#include <omp.h>
#include <limits.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>


#define xor(x,y) (x||y)&&(!(x&&y))
#define MINDOUBLE DBL_MIN
#define MAXDOUBLE DBL_MAX
