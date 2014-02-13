#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "stat.h"

/**
   This code is a C version of Homer's Statistics.pm
   http://biowhat.ucsd.edu/homer/
*/

double betacf(double a, double b, double x){
  double m;
  double m2;
  double aa;
  double c = 1.0;
  double d;
  double del;
  double h;
  double qab;
  double qam;
  double qap;

  qab = a + b;
  qap = a + 1.0;
  qam = a - 1.0;
  d = 1.0 - qab * x / qap;
  if(fabs(d) < FPMIN) d = FPMIN;
  d = 1.0 / d;
  h = d;

  for(m = 1.0; m < MAXITER; m++){
     m2 = 2.0 * m;
     aa = m * (b - m) * x / ( (qam + m2) * (a+m2) );
     d = 1.0 + aa * d;
     if(fabs(d) < FPMIN) d= FPMIN;
     c = 1.0 + aa / c;
     if(fabs(c) < FPMIN) c = FPMIN;
     d = 1.0 / d;
     h *= d * c;
     aa = -(a + m) * (qab + m) * x / ( (a+m2) * (qap+m2) );
     d = 1.0 + aa * d;
     if(fabs(d) < FPMIN) d = FPMIN;
     c = 1.0 + aa / c;
     if(fabs(c) < FPMIN) c = FPMIN;
     d = 1.0 / d;
     del = d * c;
     h *= del;
     if(fabs(del - 1.0) < EPS) break;
  }

  return(h);
}

double gammaln(double xx){
   double cof[] = {76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5396239384953e-5};
   double x, y, tmp, ser;
   int j;

   x   = xx;
   y   = xx;
   tmp = x + 5.5;
   tmp -= (x + 0.5) * log(tmp);
   ser = 1.000000000190015;

   for(j = 0; j < 6; j++) ser += cof[j] / (++y);

   return(-1.0 * tmp + log(2.5066282746310005 * ser / x));
}

double logbetai(double a, double b, double x){
   double bt;

   if(x <= 0.0 || x >= 1.0){
      printf("Invalid x in function logbetai\n");
      return(0.0);
   }

   bt = gammaln(a + b) - gammaln(a) - gammaln(b) + a * log(x) + b * log(1.0 - x);
   if(x < (a+1.0)/(a+b+2.0)){
      return bt + log(betacf(a,b,x)) - log(a);
   }
   else{
      return log(1.0 - exp(bt) * betacf(b,a,1.0 - x) / b);
   }
}

int correlation(double *x, double *y, int count, double *r, double *logp){
   double xysum = 0.0;
   double xsum  = 0.0;
   double ysum  = 0.0;
   double x2sum = 0.0;
   double y2sum = 0.0;
   double df;
   double den;
   double t;
   double numerator;
   double denomerator;
   int i;

   for(i = 0; i < count; i++){
      xysum += x[i] * y[i];
      xsum  += x[i];
      ysum  += y[i];
      x2sum += x[i] * x[i];
      y2sum += y[i] * y[i];
   }

   numerator = xysum - xsum * ysum / (double)count;
   denomerator = (x2sum - xsum * xsum / (double)count) * (y2sum - ysum * ysum / (double)count);

   if(denomerator <= 0.0) return(NA);

   *r  = numerator / sqrt(denomerator);
   df  = count - 2;
   den = ((1.0 - *r) * (1.0 + *r));
   if(den < 1e-10) den = 1e-10;
   t = *r * sqrt(df / den);
   *logp = 1.0;
   if(fabs(t) > 1e-5){
      *logp = logbetai(0.5*df, 0.5, df / (df + t * t));
   }

   return(OK);
}
