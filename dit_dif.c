/*
 dit_more.c: Decimation in Time form of FFT.

 Power of 2

 Bit reversed in transform space

 Tony Skjellum
 */

#define MAIN_DRIVER

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* preliminaries: */

#define PI 3.1415926535897932385 /* 20 digits of pi */

typedef struct _dcomplex {
  double r;
  double i;
} dcomplex;

inline dcomplex cadd(dcomplex a, dcomplex b) {
  dcomplex c;
  c.r = a.r + b.r;
  c.i = a.i + b.i;
  return c;
}

inline dcomplex csub(dcomplex a, dcomplex b) {
  dcomplex c;
  c.r = a.r - b.r;
  c.i = a.i - b.i;

  return c;
}

inline dcomplex cmul(dcomplex a, dcomplex b) {
  dcomplex c;
  /* c = (a.r + i*a.i)*(b.r + i*b.i) = (a.r*b.r - a.i*b.i) + i*(a.r*b.i + a.i*.b.r) */
  c.r = a.r * b.r - a.i * b.i;
  c.i = a.r * b.i + a.i * b.r;

  return c;
}

inline dcomplex cdiv(dcomplex a, dcomplex b) {
  dcomplex c;
  double denom;
  /* c = (a.r + i*a.i)/(b.r + i*b.i) = (a.r + i*a.i)*(b.r - i*b.i)/(b.r*b.r + b.i*b.i)  */
  denom = (b.r * b.r + b.i * b.i);
  c.r = (a.r * b.r + a.i * b.i) / denom;
  c.i = (-a.r * b.i + a.i * b.r) / denom;

  return c;
}

/* twiddle: cos(-2*pi*k/N) + i*sin(-2*pi*k/N) = cos(2*pi*k/N) - i*sin(2*pi*k/N) */
dcomplex compute_twiddle(int k, int N) /* NOTE: Not efficient, optimized, should precompute and simpify */
{
  dcomplex tmp;
  double argument = (2.0 * PI * (double) k) / (double) N;
  tmp.r = cos(argument);
  tmp.i = -sin(argument);
  return tmp;
}

void dcomplex_fft_dit(dcomplex *in, dcomplex *out, int N, int s) /* assume N is a power of 2 for now */
{
  int k;

  /* use a recursive formulation(we would not do recursively for performance) */

  if (N == 1) {
    out[0] = in[0];
  } else {
    dcomplex_fft_dit(in, out, N / 2, 2 * s);
    dcomplex_fft_dit(in + s, out + N / 2, N / 2, 2 * s);

    /* combine */
    for (k = 0; k < N / 2; ++k) {
      dcomplex t = out[k];
      dcomplex twiddle; /* compute twiddle [should not be done dynamically] */
      dcomplex term;
      twiddle = compute_twiddle(k, N);
      term = cmul(twiddle, out[k + N / 2]);
      out[k] = cadd(t, term);
      out[k + N / 2] = csub(t, term);
    }
  }
}

/*************************************/
/* REf::http://dsp-book.narod.ru/fxtbook.pdf
 /*
 /************************************/

void dcomplex_fft_dif(dcomplex *in, dcomplex *out, int N) /* assume N is a power of 2 for now */
{

  dcomplex fe[N / 2];
  dcomplex fo[N / 2];
  /* use a recursive formulation(we would not do recursively for performance) */

  if (N == 1) {
    out[0] = in[0];
  } else {
    int nh = N / 2;
    dcomplex b[N / 2];
    dcomplex c[N / 2];
    for (int k = 0; k < nh; k++) {
      dcomplex twiddle;
      twiddle = compute_twiddle(k, N);
      fe[k] = cadd(in[k], in[k + nh]);
      fo[k] = cmul(csub(in[k], in[k + nh]), twiddle);
    }
    dcomplex_fft_dif(fe, b, N / 2);
    dcomplex_fft_dif(fo, c, N / 2);
    int j = 0;
    for (int k = 0; k < nh; k++) {
      out[j] = b[k];
      out[j + 1] = c[k];
      j = j + 2;
    }
  }
}

void dcomplex_inverse_fft(dcomplex *in, dcomplex *out, int N) {
  //swap real and imaginary
  for (int i = 0; i < N; i++) {
    double temp = in[i].r;
    in[i].r = in[i].i;
    in[i].i = temp;
  }
  //fft
  dcomplex_fft_dif(in, out, N);
  //swap and divide by N
  for (int i = 0; i < N; i++) {
    double temp = out[i].r;
    out[i].r = out[i].i / N;
    out[i].i = temp / N;
  }

}

#ifdef MAIN_DRIVER

#define NEWCLEAR(_type, _N) (_type *)calloc((_N),sizeof(_type));
#define NEW(_type, N)       (_type *)malloc((_N) * sizeof(_type));
#define FREE(_p)            {free(_p); _p = 0;}

void get_dcomplex_input(dcomplex *buf, int N, char *str) {
  int i;
  for (i = 0; i < N; i++) {
    printf("Enter %s[%5d] = ", str, i);
    scanf("%lf %f", &buf[i].r, &buf[i].i);
  }
}

void print_dcomplex_output(dcomplex *buf, int N, char *str) {
  int i;
  for (i = 0; i < N; i++) {
    printf("%s[%5d] = %g + %g*i\n", str, i, buf[i].r, buf[i].i);
  }
}

void fill_test_function(int N, dcomplex *in, double (*fn)(double)) {
  int i;
  double delta = 2 * PI / (double) N, t = 0.0;
  for (i = 0; i < N; i++, t += delta) {
#define FILL_REAL
#ifdef FILL_REAL
    in[i].r = (*fn)(t);
    in[i].i = 0.0; /* in this simple test, just a real test input */
#else
    in[i].i = (*fn)(t);
    in[i].r = 0.0; /* in this simple test, just a imag test input */

#endif
  }
}

void do_test_function(int N, double (*fn)(double)) {
  dcomplex *in = NEWCLEAR(dcomplex, N)
  ;
  dcomplex *out = NEWCLEAR(dcomplex, N)
  ;

  fill_test_function(N, in, fn);
  dcomplex_fft_dit(in, out, N, 1);

  printf("\n");
  print_dcomplex_output(out, N, "output");

  FREE(out);
  FREE(in);
}
/********************************************/

void do_test_function_dif(int N, double (*fn)(double)) {
  dcomplex *in = NEWCLEAR(dcomplex, N)
  ;
  dcomplex *out = NEWCLEAR(dcomplex, N)
  ;

  fill_test_function(N, in, fn);
  dcomplex_fft_dif(in, out, N);

  printf("\n");
  print_dcomplex_output(out, N, "output");

  FREE(out);
  FREE(in);
}

void do_main(int N) {
  dcomplex *in = NEWCLEAR(dcomplex, N)
  ;
  dcomplex *out = NEWCLEAR(dcomplex, N)
  ;

  if (in != 0 || out != 0) {
    get_dcomplex_input(in, N, "input");

    /* do the FFT: */
    dcomplex_fft_dif(in, out, N);

    printf("\n");
    print_dcomplex_output(out, N, "output");
  }

  FREE(out);
  FREE(in);
}

void do_mainloop(int N) {
  dcomplex *in = NEWCLEAR(dcomplex, N)
  ;
  dcomplex *out = NEWCLEAR(dcomplex, N)
  ;

  if (in != 0 || out != 0) {
    get_dcomplex_input(in, N, "input");

    /* do the FFT: */
    dcomplex_fft_dit(in, out, N, 1);

    printf("\n");
    print_dcomplex_output(out, N, "output");
  }

  FREE(out);
  FREE(in);
}

double my_function(double t) {
  return (sin(t));
}

double test_function(double t) {
  return (sin(t) + sin(2 * t) + sin(3 * t));
}

int main(int argc, char **argv) {
  int i;
  int N;

//  do {
//    printf("Enter the length of the vector (0 <= N, 0 to exit): ");
//    scanf("%i", &N);
//  } while (N < 0);
//  if (N > 0) {
//    do_mainloop(N);
//  }

  do {
    printf("Enter the length of the vector (0 <= N, 0 to exit): ");
    scanf("%i", &N);
  } while (N < 0);
  if (N > 0) {
    do_main(N);
  }

//sin(t)+sin(2t)+sin(3t) function test;
  printf("\n Function: sint+sin2t+sin3t testing");
#define MINLEN 4
#define MAXLEN 8
  for (i = MINLEN; i <= MAXLEN; i *= 2) {
    printf("\nTest function for %d: \n", i);
    printf("\n dit: \n");
    do_test_function(i, test_function);
    printf("\n dif: \n");
    do_test_function_dif(i, test_function);
  }
  printf("\n");

  printf("Inverse DFT test:");

  dcomplex x[4],  output[4];
  x[0].r = 10;
  x[0].i = 0;
  x[1].r = -2;
  x[1].i = 2;
  x[2].r = -2;
  x[2].i = 0;
  x[3].r = -2;
  x[3].i = -2;
  printf("\n");
  print_dcomplex_output(x, 4, "input");
  dcomplex_inverse_fft(x,output,4);
  printf("\n");
  print_dcomplex_output(output, 4, "output");

  return 0;
}

#endif

