/*
 * trap floating point exceptions
 *
 * Thomas Neff 2003
*/


#define _GNU_SOURCE 1
#include <fenv.h>

static void __attribute__ ((constructor)) trapfpe(void)
{
	  /* Enable some exceptions.  At startup all exceptions are masked. */
	  feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
}
