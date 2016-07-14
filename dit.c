/*
   Decimation in Time form of FFT.

   Power of 2

   Bit reversed in transform space

   Tony Skjellum

   intermediate form of 10/30/13 19:41
*/

#include <math.h>

/* preliminaries: */

#define PI 3.1415926535897932385 /* 20 digits of pi */

typedef struct _dcomplex
{
  double r;
  double i;
} dcomplex;

dcomplex cadd(dcomplex a, dcomplex b) 
{
  dcomplex c;
  c.r = a.r + b.r;
  c.i = a.i + b.i;
  return c;
}

dcomplex csub(dcomplex a, dcomplex b) 
{
  dcomplex c;
  c.r = a.r - b.r;
  c.i = a.i - b.i;

  return c;
}

dcomplex cmul(dcomplex a, dcomplex b) 
{
  dcomplex c;
  /* c = (a.r + i*a.i)*(b.r + i*b.i) = (a.r*b.r - a.i*b.i) + i*(a.r*b.i + a.i*.b.r) */
  c.r = a.r*b.r - a.i*b.i;
  c.i = a.r*b.i + a.i*b.r;

  return c;
}

dcomplex cdiv(dcomplex a, dcomplex b) 
{
  dcomplex c;
  double denom;
  /* c = (a.r + i*a.i)/(b.r + i*b.i) = (a.r + i*a.i)*(b.r - i*b.i)/(b.r*b.r + b.i*b.i)  */
  denom = (b.r*b.r + b.i*b.i);
  c.r =( a.r*b.r + a.i*b.i)/denom;
  c.i =(-a.r*b.i + a.i*b.r)/denom;

  return c;
}

/* twiddle: cos(-2*pi*k/N) + i*sin(-2*pi*k/N) = cos(2*pi*k/N) - i*sin(2*pi*k/N) */
dcomplex compute_twiddle(int k, int N)
{
  dcomplex tmp;
  double argument = (2.0*PI*(double)k)/(double)N;
  tmp.r = cos(argument);
  tmp.i = -sin(argument);
  return tmp;
}


void dcomplex_fft_dit(dcomplex *in, dcomplex *out, int N, int s) /* assume N is a power of 2 for now */
{
  int k;

  /* use a recursive formulation */

  if(N==1)
  {
    out[0] = in[0];
  }   
  else
  {
    dcomplex_fft_dit(in, out, N/2, 2*s);
    dcomplex_fft_dit(in+s, out+N/2, N/2, 2*s);
    
    /* combine */
    for(k = 0; k < N/2 ; ++k)
    {
       dcomplex t = out[k];
       dcomplex twiddle; /* compute twiddle [should not be done dynamically] */
       dcomplex term;
       twiddle = compute_twiddle(k,N);
       term    = cmul(twiddle, out[k+N/2]);
       out[k]     = cadd(t, term);
       out[k+N/2] = csub(t, term);
    }    
  }

}




