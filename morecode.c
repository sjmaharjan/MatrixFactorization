/* Note this is the slow test code: */
void dcomplex_dft(dcomplex *in, dcomplex *out, int N, int s) /* assume N is a power of 2 for now */
{
   int k,m;

   for(k = 0; k < N; k++)
   {
      dcomplex tmp;  tmp.r = tmp.i = 0.0;

      for(m = 0; m < N; m++)
      {
	if(k == 0 || m == 0) 
          tmp = cadd(tmp, in[m]);
        else
        {
	  /* w^(k*m) : */
	  double arg = 2.0*PI*(double)k*(double)m/(double)N;
          dcomplex w; w.r= cos(arg); w.i = -sin(arg);
	  dcomplex val; 
	  val = cmul(w, in[m]);
	  tmp = cadd(tmp, val);
        }
      }  
      out[k] = tmp;
   }
}

void print_dcomplex_comparison(dcomplex *buf1, dcomplex *buf2, int N, char *str1, char *str2)
{
  int i;
  for(i = 0; i < N; i++)
  {
    double res;
    dcomplex tmp1 = csub(buf1[i],buf2[i]);
    res = sqrt(tmp1.r*tmp1.r + tmp1.i*tmp1.i);

    printf("[%5d] = %s (%g + %g*i) vs. (%g + %g*i) %s, residual =%g\n", i, str1, buf1[i].r, buf1[i].i,
	   buf2[i].r, buf2[i].i, str2, res);

  }
}

