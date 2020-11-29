//
// Libraries for table generation
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>

#include "defs.h"
#include "gentable.h"

#include "chebyshev.h"

extern int B_MDG_CHECK;
extern int MDG_DEBUG_FLAG;

#ifdef DEBUG
#define DEBUG_PRINTF(...) { if (MDG_DEBUG_FLAG) printf(__VA_ARGS__); }
#else
#define DEBUG_PRINTF(...) 
#endif

#ifdef DEBUG
static double get_min(int n, double dmin, int ndiv_man) {
  int k = n>>ndiv_man;
  double xmin = 0.0;
  for(int i=0;i<k;++i) {
    xmin += dmin;
    if (i>0) dmin*=2.0;
  }
  n -= k*(1<<ndiv_man);
  xmin += n*dmin/(1<<ndiv_man);
  return xmin;
}

static int get_n(double x, double dmin, int ndiv_man) {
  int n=0;
  while (x>=dmin) {
    x    -= dmin;
    if (n>0) dmin*=2.0;
    n += (1<<ndiv_man);
  }
  n += (1<<ndiv_man)*x/dmin;
  return n;
}

static double get_w(int n, double dmin, int ndiv_man) {
  int k = n>>ndiv_man;
  for(int i=1;i<k;++i) {
    dmin*=2.0;
  }
  return dmin/(1<<ndiv_man);
}
#endif

void
mdg_table_init(mdg_func_table *t, int norder, int nmax, const char *name, int allow_carryup) {
  t->norder = norder;
  t->max_entries = nmax;
  t->nentries=0;
  t->nprec=0;
  t->floating=0;
  t->exporg=0;
  t->xmin = (double *)malloc(sizeof(double)*nmax);
  t->w    = (double *)malloc(sizeof(double)*nmax);
  t->c    = (double *)malloc(sizeof(double)*nmax*(norder+1));
  t->c1   = (double *)malloc(sizeof(double)*nmax*(norder+2));
  t->cp   = (double *)malloc(sizeof(double)*nmax*(norder+1));
  t->ci   = (int64_t *)malloc(sizeof(int64_t)*nmax*(norder+1));
  t->exp  = (int *)malloc(sizeof(int)*nmax);
  memset(t->exp, 0, sizeof(int)*nmax);
  t->cbit = (int *)malloc(sizeof(int)*(norder+1));
  t->sign = (int *)malloc(sizeof(int)*(norder+1));
  t->hsign = (int *)malloc(sizeof(int)*(norder+1)); // Sign for Horner's method
  t->name = (char *)malloc(sizeof(char)*strlen(name));
  t->allow_carryup = allow_carryup;
  strcpy(t->name, name);
}

// Add n intervals of the same width w starting from xmin
void
mdg_table_add_intervals
(
 double (*f)(double),
 int n,
 double xmin,
 double w,
 mdg_func_table *t) {
  if (t->nentries+n>t->max_entries) {
    fprintf(stderr, "Error: table overflow.\n");
    exit(1);
  } else {
    double x=xmin;
    double *c=t->c + (t->norder+1)*t->nentries;
    double *c1=t->c1 + (t->norder+2)*t->nentries;
    double *cp=t->cp + (t->norder+1)*t->nentries;
    double *xminp = t->xmin + t->nentries;
    double *wp    = t->w + t->nentries;
    for(int i=0;i<n;++i) {
      mdg_cheb_fit(f, x, x+w, t->norder, c);
      mdg_cheb_fit(f, x, x+w, t->norder+1, c1);
      // Conversion to polynomial coefficients
      mdg_cheb_conv(t->norder,c,cp);
      *xminp = x;
      *wp = w;
      x+=w;
      c+=t->norder+1;
      cp+=t->norder+1;
      c1+=t->norder+2;
      ++xminp;
      ++wp;
    }
    t->nentries += n;
  }
}

static void
check_sign
(
 mdg_func_table *t
 )
{
  int n=t->norder+1;
  for(int i=0;i<n;++i) {
    int s = 0;
    for(int j=0;j<t->nentries;++j) {
      double c = t->cp[j*n+i];
      if (s==0) {
	s =
	  ( c > 0.0) ? 1 :
	  ( ( c < 0.0 ) ? -1 : 0);
      } else if (s*c < 0.0) {
	s = 0; break;
      }
    }
    t->sign[i]=s;
  }

  t->hsign[0] = t->sign[0];
  for(int i=1;i<=t->norder;++i) {
    t->hsign[i] = t->sign[i] * t->sign[i-1];
  }
}
  
void
mdg_table_integer_conversion 
(
 int nprec,     // precision for fractions
 mdg_func_table *t
 ) {
  t->nprec=nprec;
  t->floating=0;
  //  int cbit_total=0;
  check_sign(t);
  for(int i=0;i<=t->norder;++i) {
    double cmin=t->cp[i];
    double cmax=t->cp[i];
    double cabsmax;
    for(int j=1;j<t->nentries;++j) {
      int k=j*(t->norder+1);
      if (t->cp[k+i]>cmax) cmax=t->cp[k+i];
      if (t->cp[k+i]<cmin) cmin=t->cp[k+i];
    }
    cabsmax = (fabs(cmax)>fabs(cmin)) ? fabs(cmax) : fabs(cmin);

    // ilogb returns 0 when [1,2)
    // c in [1/2,1)   -> Nprecision - Ntable bit
    // c in [1,2)     -> Nprecision - Ntable + 1bit
    // c in [1/4,1/2) -> Nprecision - Ntable - 1bit
    // Here, we always include sign bit, so plus 1.
    t->cbit[i] = nprec + ilogb(cabsmax) + 1 + 1;  

    if(MDG_DEBUG_FLAG)fprintf(stderr, "Min: %le Max: %le Absmax : %le Nbit : %d sign: %d\n",
	    cmin, cmax, cabsmax, t->cbit[i], t->sign[i]);
  }

  double *cp = t->cp;
  int64_t *cip = t->ci;
  for(int i=0;i<t->nentries;++i) {
    for(int j=0;j<=t->norder;++j) {
      *cip = lround(scalbn(*cp,nprec));
      ++cip;
      ++cp;
    }
  }
}

void
mdg_table_integer_conversion_coulomb
(
 int one,     // value for "1"
 mdg_func_table *t
 ) {
  t->nprec=one;
  t->floating=0;
  //  int cbit_total=0;
  check_sign(t);
  for(int i=0;i<=t->norder;++i) {
    double cmin=t->cp[i];
    double cmax=t->cp[i];
    double cabsmax;
    for(int j=1;j<t->nentries;++j) {
      int k=j*(t->norder+1);
      if (t->cp[k+i]>cmax) cmax=t->cp[k+i];
      if (t->cp[k+i]<cmin) cmin=t->cp[k+i];
    }
    cabsmax = (fabs(cmax)>fabs(cmin)) ? fabs(cmax) : fabs(cmin);

    // ilogb returns 0 when [1,2)
    // c in [1/2,1)   -> Nprecision - Ntable bit
    // c in [1,2)     -> Nprecision - Ntable + 1bit
    // c in [1/4,1/2) -> Nprecision - Ntable - 1bit
    // Here, we always include sign bit, so plus 1.
    t->cbit[i] = ilogb(cabsmax*one) + 1 + 1;  

    if(MDG_DEBUG_FLAG)fprintf(stderr, "Min: %le Max: %le Absmax : %le Nbit : %d sign: %d\n",
	    cmin, cmax, cabsmax, t->cbit[i], t->sign[i]);
  }

  double *cp = t->cp;
  int64_t *cip = t->ci;
  for(int i=0;i<t->nentries;++i) {
    for(int j=0;j<=t->norder;++j) {
      *cip = lround(*cp * one );
      ++cip;
      ++cp;
    }
  }
}

void
mdg_table_integer_conversion_float
(
 int nprec,    // precision for fractions
 mdg_func_table *t
 ) {
  t->nprec=nprec;
  t->floating=1;
  //  int cbit_total=0;
  int n=t->norder+1;
  // check whether sign of each coefficient changes 
  int *sign = t->sign;
  for(int i=0;i<n;++i) sign[i]=2;
  for(int j=0;j<t->nentries;++j) {
    for(int i=0;i<n;++i) {
      if (sign[i]==2) {
	sign[i] =
	  ((t->cp[j*n+i]) > 0) ? 1 :
	  ( ((t->cp[j*n+i]) < 0) ? -1 : 0);
      } if ( (sign[i]==1)||(sign[i]==-1) ) {
	if (sign[i]*t->cp[j*n+i] < 0.0) {
	  sign[i] = 0;
	}
      }
    }
  }
  for(int i=0; i<n; ++i) t->cbit[i] = 0;

  double *cp = t->cp;
  int64_t * cip = t->ci;
  t->cbit[0]=nprec+2; // with sign
  for(int j=0;j<t->nentries;++j) {
    // In some case c0<cn, so search the max coefficient
    double cmax=0.0;
    int imax=-1;
    for(int i=0;i<=t->norder;++i) {
      if (fabs(cp[i])>cmax) {
	cmax = fabs(cp[i]);
	imax = i;
      }
    }
    if(MDG_DEBUG_FLAG)fprintf(stderr, "cmax %f at %d\n",cmax,imax);
    
    // note : if c0==-2^k, it will fit in k+1 bits, but here k+2 bits will be used.
    int c0exp = ilogb(cmax); // 2^e <= c0 < 2^(e+1)
    // scale by 2^(nprec-1-e) (without sign case)
    int shift = nprec-c0exp;
    t->exp[j] = -shift;
    
    for(int i=0;i<=t->norder;++i) {
      double cn = *cp;
      cn = scalbn(cn, shift);
      int cnbits = ilogb(fabs(cn));
      if (cn==-(1<<cnbits)) cnbits +=1;
      else cnbits +=2;
      if (cnbits > t->cbit[i]) t->cbit[i]=cnbits;
      *cip = llround(cn);
      ++cp;
      ++cip;
    }
  }

  for(int i=0;i<n;++i) {
    if(MDG_DEBUG_FLAG)fprintf(stderr, "order %d: Nbit %d sign: %d\n",
	    i, t->cbit[i], t->sign[i]);
  }

  t->hsign[0] = t->sign[0];
  for(int i=1;i<=t->norder;++i) {
    t->hsign[i] = t->sign[i] * t->sign[i-1];
  }

  // Process exponents
  int expmax = t->exp[0];
  int expmin = expmax;
  for(int j=1;j<t->nentries;++j) {
    if (t->exp[j] > expmax) expmax = t->exp[j];
    if (t->exp[j] < expmin) expmin = t->exp[j];
  }
  int dexpn = expmax-expmin+1; // include exponent for 0
  int expbits = 0;
  while (dexpn) {
    ++expbits;
    dexpn>>=1;
  }
  t->expbits = expbits;
  for(int j=0;j<t->nentries;++j)
    t->exp[j] += -expmin + 1; // expmin = 1
  t->exporg = expmin-1;

  if(MDG_DEBUG_FLAG)fprintf(stderr, "exponent Nbit %d: Max %d\n", expbits, expmax-expmin);
}

static void
get_minmax(mdg_func_table *t, double *xmin, double *xmax)
{
  *xmin = *xmax = t->xmin[0];
  for(int i=0;i<t->nentries;++i) {
    if (t->xmin[i]+t->w[i]>*xmax) *xmax = t->xmin[i]+t->w[i];
    else if (t->xmin[i]<*xmin) *xmin = t->xmin[i];
  }
}

// Evaluation in double
double
mdg_table_eval(double x, int k, mdg_func_table *t) {
  double *c = t->cp+k*(t->norder+1);
  double d = 2.0*(x - t->xmin[k])/t->w[k] - 1.0; // d = -1~1

  double z=c[t->norder];
  for(int j=t->norder-1;j>=0;--j)
    z = z*d + c[j];

  return z;
}

// Evaluation in integer. input and output are also integer.
int64_t
mdg_table_eval_ii(int64_t dx, int dxprec, int k, mdg_func_table *t) {
  // Note: in C, a result of division is rounded toward 0.
  // That's why I use double for this calculation..

  double dxf = dx;
  int64_t *c = t->ci+k*(t->norder+1);
  int *cbit = t->cbit;
  double dxl;
  double dxw;

  double z = c[t->norder];
  for (int i=t->norder; i>0; --i) {
    if (dxprec > cbit[i]) {
      dxw = 1<<cbit[i];
      dxl = floor(dxf / (1<<(dxprec-cbit[i])));
    } else {
      dxw = 1<<dxprec;
      dxl = dxf;
    }
    double z0 = z;
    z = c[i-1] + floor(dxl * z / dxw);
    if(MDG_DEBUG_FLAG)printf("%d: %"PRIx64" %"PRIx64" %"PRIx64" %f %f\n", i, (int64_t) z, c[i-1], (int64_t)z0, dxl, dxw);
    int cw = cbit[i-1];
    if (t->allow_carryup) ++cw; // allow carryup after addition
    int64_t cmax, cmin;
    cmax = (1<<(cw-1))-1;
    cmin = -cmax-1;
    if ( (z>cmax) || (z<cmin) ) {
      double x = t->xmin[k] + t->w[k]*( scalbn(dxf,-dxprec)+1.0);
      fprintf(stdout, "INFO (%s) : Range error: x=%f, n=%d, dx=%"PRIx64", z=%"PRIx64"(%le), c=%"PRIx64", cmin=%"PRIx64", cmax=%"PRIx64"\n",
	      __func__, x, i-1, dx & MASK(dxprec+1), (int64_t) z,z,
	      c[i-1] & MASK(cbit[i-1]),
	      cmin & MASK(cbit[i-1]),
	      cmax & MASK(cbit[i-1]));
    }
  }
  return (int64_t) z;
}

// Evaluation in integer. input and output are double.
double
mdg_table_eval_i(double x, int dxprec, mdg_func_table *t) {
  int k=-1;
  for(int i=0; i<t->nentries; ++i) {
    if ((x>=t->xmin[i])&&(x<t->xmin[i]+t->w[i])) {
      k=i;
      break;
    }
  }
  if (k<0) {
    fprintf(stdout, "ERROR (%s): input range error : x=%f, min=%f, max=%f\n",
	    __func__, x, t->xmin[0], t->xmin[t->nentries-1]+ t->w[t->nentries-1]);
    return 0;
  }
  double d = 2.0*(x - t->xmin[k])/t->w[k] - 1.0; // d = -1~1
  int64_t dx = scalbn(d, dxprec);

  int64_t z = mdg_table_eval_ii(dx,dxprec, k, t);

  double zd;
  if (t->floating) {
    zd = scalbn((double)z, t->exp[k]+t->exporg);
  } else {
    zd = scalbn((double)z, -t->nprec);
  }
  return zd;
}

#ifndef GNUPLOT_EXE
#define GNUPLOT_EXE "/usr/local/bin/gnuplot"
//#define GNUPLOT_EXE "/usr/bin/gnuplot"
#endif

void
mdg_table_plot(int n, int ylog, double (*f)(double), mdg_func_table *t)
{
  char *datafilename = (char *) malloc(sizeof(char)*(strlen(t->name)+100));
  char *plotfilename = (char *) malloc(sizeof(char)*(strlen(t->name)+100));
  strcpy(datafilename, t->name);
  strcpy(plotfilename, t->name);
  strcat(datafilename, "_err.dat");
  strcat(plotfilename, "_err.gp");
  FILE *fp = fopen(datafilename, "w");
  double xmin, xmax;
  get_minmax(t, &xmin, &xmax);

  // Index 0: Theoretical Error, absolute and relative
  for(int i=0; i<t->nentries; ++i) {
    double err = fabs(t->c1[(t->norder+2)*(i+1)-1]);
    double errlsb = scalbn(err, t->nprec);
    if (t->floating) {
      errlsb = scalbn(err, -(t->exp[i]+t->exporg));
    }
    double relerr = 0.0;
    double x = t->xmin[i]+t->w[i]*0.5;
    double exact = (*f)(x);
    if (exact!=0.0) relerr = fabs(err/exact);
    fprintf(fp,"%f %le %le %le\n", x, err, errlsb, relerr);
  }
  fprintf(fp,"\n\n");

  // Index 1: Error based on calculations in double and int
  fprintf(fp,"#1:x 2:err(dbl) 3:err(int) 4:err(dbl/lsb) 5:err(int/lsb) 6:y(exact) 7:y(d) 8:y(i) 9:relerr(d) 10:relerr(i)\n");
  for(int i=0; i<t->nentries; ++i) {
    double x = t->xmin[i];
    double w = t->w[i]/n;
    double lsb = scalbn(1.0, t->exp[i]+t->exporg);
    for(int j=0; j<n; ++j) {
      double y_exact  = (*f)(x);
      double y_poly_d = mdg_table_eval(x,i,t);
      double y_poly_i = mdg_table_eval_i(x,18,t);
      double err_d = fabs(y_poly_d-y_exact);
      double err_i = fabs(y_poly_i-y_exact);
      double relerr_d = 0.0;
      double relerr_i = 0.0;
      if (y_exact!=0.0) {
	relerr_d = fabs(err_d/y_exact);
	relerr_i = fabs(err_i/y_exact);
      }
      fprintf(fp,"%f %le %le %le %le %le %le %le %le %le\n",x,
	      err_d, err_i, err_d/lsb, err_i/lsb, y_exact, y_poly_d, y_poly_i, relerr_d, relerr_i);
      x += w;
    }
  }
  fprintf(fp,"e\n");
  fclose(fp);

  //fp = popen(GNUPLOT_EXE, "w");
  //pclose(fp);
  fp = fopen(plotfilename, "w");

  fprintf(fp,"set xrange [%f:%f]\n", xmin, xmax);
  if (ylog) fprintf(fp,"set logscale y\n");
  fprintf(fp,"file='%s'\n", datafilename);
  fprintf(fp,"plot file index 0 using 1:2 with steps lt rgb 'black' t 'Theoretical',");
  //  fprintf(fp,"     file index 1 using 1:2 with lines lt rgb 'blue' ,");
  fprintf(fp,"     file index 1 using 1:3 with lines lt rgb 'red' t 'Actual'\n");
  fprintf(fp,"pause -1\n");
  fprintf(fp,"plot file index 0 using 1:3 with steps lt rgb 'black' t 'Theoretical/lsb',");
  //  fprintf(fp,"     file index 1 using 1:4 with lines lt rgb 'blue' ,");
  fprintf(fp,"     file index 1 using 1:5 with lines lt rgb 'red' t 'Actual/lsb'\n");
  fprintf(fp,"pause -1\n");
  fprintf(fp,"plot file index 0 using 1:4 with steps lt rgb 'black' t 'Theoretical relative',");
  fprintf(fp,"     file index 1 using 1:10 with lines lt rgb 'red' t 'Actual relative'\n");
  fprintf(fp,"pause -1\n");
  if (ylog) fprintf(fp,"set nologscale y\n");
  fprintf(fp,"plot file index 0 using 1:4 with steps lt rgb 'black' t 'Theoretical relative'\n");
  fprintf(fp,"pause -1\n");
  
  fclose(fp);
  
}

// VERILOG Generation
//   if sign of each coefficient changes, sign bit is automatically added
//     and output signed value; in the case, with_sign has no effect.
//   sign_mode=0 : always include sign bit
//   sign_mode=1 : 2's complement and no sign bit (if possible)
//   sign_mode=2 : absolute and no sign bit (if possible)
void 
mdg_table_output_verilog
(char *name, int sign_mode, mdg_func_table *t, char *command) {
  char *s = (char *)malloc(sizeof(char)*(strlen(name)+3));
  strcpy(s, name);
  strcat(s, ".v");
  FILE *fp = fopen(s, "w");
  if (!fp) {
    fprintf(stderr, "Output file %s cannot open.\n", name);
    exit(1);
  }
  int nbit_adr = ilogb(t->nentries);
  while (t->nentries>(1<<nbit_adr)) ++nbit_adr;
  for(int i=0;i <= t->norder; ++i) {
    int width = t->cbit[i];
    int sign = 1;
    if ( sign_mode && t->sign[i]) {
      if (sign_mode==2) sign = t->sign[i];
      --width;
    }

    fprintf(fp,"//  %s\n",command);
    fprintf(fp,"//  sign : %d\n\n",sign);
    fprintf(fp,"module %s_%d(input [%d:0] adr, output reg [%d:0] c);\n",
	    name, i, nbit_adr-1, width-1);
    fprintf(fp,"  always_comb begin\n");
    fprintf(fp,"    case (adr)\n");
    for(int j=0;j<t->nentries; ++j) {
      fprintf(fp,"      %d'h%x: c = %d'h%llx;\n",
	      nbit_adr, j, width, (sign*t->ci[j*(t->norder+1)+i]) & MASK64(width));
    }
    fprintf(fp,"      default: c = %d'bx;\n",width);
    fprintf(fp,"    endcase;\n");
    fprintf(fp,"  end\n");
    fprintf(fp,"endmodule // %s_%d\n\n", name, i);
  }
  fclose(fp);
}

// C Generation
void
mdg_table_output_c(char *name, mdg_func_table *t, char *command) {
  char *s = (char *)malloc(sizeof(char)*(strlen(name)+3));
  strcpy(s, name);
  strcat(s, ".h");
  FILE *fp = fopen(s, "w");
  if (!fp) {
    fprintf(stderr, "Output file %s cannot open.\n", name);
    exit(1);
  }

  fprintf(fp,"//  %s\n\n",command);
  fprintf(fp,"#define NORDER_%s %d\n", name, t->norder);
  fprintf(fp,"#define ADRMAX_%s %d\n", name, t->nentries);
  fprintf(fp,"#define NAME_%s \"%s\"\n", name, name);
  fprintf(fp,"int %s_sign[%d] = {", name, t->norder+1);
  for(int i=0;i <= t->norder; ++i) {
    fprintf(fp," %d", t->sign[i]);
    if (i==t->norder) fprintf(fp,"};\n"); else fprintf(fp,",");
  }
  fprintf(fp,"int %s_cbit[%d] = {", name, t->norder+1);
  for(int i=0;i <= t->norder; ++i) {
    fprintf(fp," %d", t->cbit[i]);
    if (i==t->norder) fprintf(fp,"};\n"); else fprintf(fp,",");
  }

  int int64flag = 0;
  for(int i=0;i <= t->norder; ++i) 
    if (t->cbit[i]>31) int64flag = 1;
  if (int64flag) {
    fprintf(fp,"int64_t ");
  } else {
    fprintf(fp,"int ");
  }

  fprintf(fp, "%s_table[%d][%d] = {\n",name, t->nentries, t->norder+1);
  
  int64_t *ci = t->ci;
  for(int j=0;j<t->nentries; ++j) {
    fprintf(fp,"  { ");
    for(int i=0;i <= t->norder; ++i) {
      fprintf(fp, "%"PRId64, *ci);
      ++ci;
      if (i==t->norder) fprintf(fp,"}");
      else fprintf(fp, ", ");
    }
    if (j<t->nentries-1) fprintf(fp,",");
    fprintf(fp,"// [%f, %f) d=%f\n",
	    t->xmin[j], t->xmin[j]+t->w[j], t->w[j]);
  }
  fprintf(fp,"};\n");

  if (t->floating) {
    fprintf(fp, "%s_exp_table[%d] = {\n",name, t->norder+1);
  
    for(int j=0;j<t->nentries; ++j) {
      fprintf(fp, "%d", t->exp[j]);
      if (j==t->nentries-1) fprintf(fp,"};\n");
      else fprintf(fp,",\n");
    }
  }
  fclose(fp);
}

