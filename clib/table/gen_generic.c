//
// Table Generation for general function
//
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <dlfcn.h>

#include "gentable.h"
int MDG_DEBUG_FLAG;	//Added (Komatsu) I hope no side effect happen due to this addition.
#include "fpmd.h"	//Added (Komatsu) for debug

//double log2(double x) { return log(x)*M_LOG2E; }

char *c_filename = "/tmp/tbl_func.c";
char *so_filename = "/tmp/tbl_func.so";
char *Funcname = "tbl_func";

int compile(char *func)
{
  int i;
  char cmd[1024];
  FILE *fp;

  //snprintf(s, 256, "/tmp/mktbl_func_%d.c",getpid());
  fprintf(stderr, "Temporal file name = %s\n", c_filename);
  fp = fopen(c_filename, "w");
  if (!fp) {
    fprintf(stderr, "Error: cannot open temporal file %s\n", c_filename);
    return 1;
  }
  fprintf(fp, "#include <stdio.h>\n");
  fprintf(fp, "#include <math.h>\n");
  //fprintf(fp, "void _init(void) {}\n");
  //fprintf(fp, "void _fini(void) {}\n");
  fprintf(fp, "double %s (double x) {\n", Funcname);
  fprintf(fp, "  return %s;\n", func);
  fprintf(fp, "}\n");
  fclose(fp);

  //  snprintf(cmd, 1024, "osarch=`uname -s`-`uname -m`; cc -g -O2 -fPIC -shared -o %s -lc -lm -L ../lib/$osarch -rpath ../lib/$osarch -lmd %s",
  snprintf(cmd, 1024, "osarch=`uname -s`-`uname -m`; cc -g -O2 -fPIC -shared -o %s -lc -lm -L `hg root`/lib/$osarch -Wl,-rpath `hg root`/lib/$osarch -lmd %s",
	  so_filename, c_filename);

  if (system(cmd)) {
    fprintf(stderr, "Error: in compile\n");
    return 2;
  }
  return 0;
}

void *Handle;
double (* Func)(double);

int dynamic_load_func()
{
  char *error;

  Handle=dlopen(so_filename, RTLD_NOW);
  if (!Handle) {
    fprintf(stderr, "Error: opening shared object file %s\n", so_filename);
    return 1;
  }
  dlerror();
  *(void **) (&Func) = dlsym(Handle, Funcname);
  if ((error=dlerror())) {
    fputs(error, stderr);
    return 2;
  }
  return 0;
}

void usage(char *progname)
{
  printf("%s: Generate Function Table\n",progname);
  printf("Options\n");
  printf("  -f func : function definition in C\n");
  printf("  -b #    : bit width for fraction part of resuls\n");
  printf("  -o #    : polynomial order\n");
  printf("  -w #    : number of segment = 2^#\n");
  printf("  -e #    : input lower bound (log) = 2^-#\n");
  printf("  -l #    : input lower bound\n");
  printf("  -r #    : input upper bound\n");
  printf("  -s #    : input floating mode\n");
  printf("            0: fixed, 1: floating, 2: floating with an additional segment [0,xmin)\n");
  printf("  -t #    : output floating mode\n");
  printf("            0: fixed, 1: floating\n");
  printf("  -z #    : sign mode for verilog output\n");
  printf("            0: always include sign bit\n");
  printf("            1: no sign bit added but number is in 2's compliment\n");
  printf("            2: absolute value\n");
  printf("            if coefficient has both signs, a sign bit will be added automatically.\n");
  printf("  -n name : table name\n");
  printf("  -h      : print this help\n");
}

#define MAXLEN_NAME 1024
#define MAXLEN_FUNC 2048
#define NDIVMAX_FLOAT 4096

int
main(int argc, char *argv[]) {
  
  int nbit    = 24;
  int norder  = 2;
  double xmin = 1.0;
  double xmin_e = 10.0;
  double xmax = 2.0;
  int ndivlog = 6;
  char name[MAXLEN_NAME];
  char func[MAXLEN_FUNC];
  int floating = 0;
  int output_floating = 0;

  int sign_mode = 0;
  int allow_carryup = 0;

  mdg_func_table t;

  strcpy(name, "log_table");
  strcpy(func, "log2(x)");

  int optlen = 0;
  for(int i=0;i<argc;++i) optlen += strlen(argv[i]);
  optlen += argc+1;
  char *cmdline = (char *)malloc(sizeof(char)*optlen);
  strcpy(cmdline, argv[0]);
  for(int i=1;i<argc;++i) {
    strcat(cmdline, " ");
    strcat(cmdline, argv[i]);
  }

  //  printf("%s\n", cmdline);

  char ch;
  while ((ch=getopt(argc, argv, "chb:n:o:w:e:l:r:s:t:f:z:"))!=-1) {
    switch (ch) {
    case 'c' : allow_carryup=1; break;
    case 'h' : usage(argv[0]); exit(0);
    case 'b' : sscanf(optarg, "%d", &nbit);     break;
    case 'n' : strncpy(name, optarg, MAXLEN_NAME); break;
    case 'o' : sscanf(optarg, "%d", &norder);   break;
    case 'w' : sscanf(optarg, "%d", &ndivlog);  break;
    case 'e' : 
      sscanf(optarg, "%lf", &xmin_e);
      xmin = pow(2.0, -xmin_e);
      break;
    case 'l' : sscanf(optarg, "%lf", &xmin);    break;
    case 'r' : sscanf(optarg, "%lf", &xmax);    break;
    case 's' : sscanf(optarg, "%d", &floating); break;
    case 't' : sscanf(optarg, "%d", &output_floating); break;
    case 'f' : 
      if (strlen(optarg)>MAXLEN_FUNC) {
	fprintf(stderr,"Error: function definition too long.\n");
	exit(1);
      } else {
	strncpy(func, optarg, MAXLEN_FUNC);
      } break;
    case 'z' :
      sscanf(optarg, "%d", &sign_mode);
      if ((sign_mode>2)||(sign_mode<0)) {
	fprintf(stderr,"Error: illegal sign mode %d\n", sign_mode);
	exit(1);
      }
      break;
    }
  }

  int ndiv    = 1<<ndivlog;
  name[MAXLEN_NAME-1]=0;

  compile(func);
  if (dynamic_load_func()) {
    fprintf(stderr,"Error: dynamic load failed.\n");
    exit(1);
  }

 if (floating) {
    printf("Input floating mode\n");
    printf("x range = [%f,%f)\n", xmin, xmax);
    // for example:
    //  (xmin, xmax) = (0.5, 2.0)
    //   -> (e_min, e_max) = (-1, 1)
    //  (xmin, xmax) = (0.5, 2.1)
    //   -> (e_min, e_max) = (-1, 2)
    // e_max is actually e_max+1 
    int e_min = ilogb(xmin);
    int e_max = ceil(log2(xmax));
    int de = e_max-e_min;
    if (floating==2) {
      printf("mode 2: [0, xmin) is included as a fixed-point range.\n");
      de++;
    }
    int wexp = ceil(log2(de));
    int wman = ndivlog - wexp;
    printf("exp_min=%d exp_max=%d exp_bits=%d man_bits=%d\n", e_min, e_max, wexp, wman);
    // if 2^wexp > e_max-e_min, we can extend range for large or small values.
    // Here, we always extend for smaller values than xmin.
    if ((1<<wexp)>de) {
      e_min = e_max - (1<<wexp);
      if (floating==2) ++e_min;
      printf("new exp_min=%d\n", e_min);
    }
    double xmin_new = scalbn(1.0, e_min);
    double xmax_new = scalbn(1.0, e_max);
    printf("new range = [%f,%f)\n", xmin_new, xmax_new);
    printf("new range sqrt = [%f,%f)\n", sqrt(xmin_new), sqrt(xmax_new));

    //init_table(&t, norder, ndiv);
    mdg_table_init(&t, norder, NDIVMAX_FLOAT, name, allow_carryup);

    int nblock = 1<<wman;
    printf("man_bits=%d nblock=%d\n", wman, nblock);
    int n=0;
    if (floating==2) {
      // Add interval from 0~xmin
      double w = scalbn(1.0, e_min-wman);
      mdg_table_add_intervals(Func, nblock, 0.0, w, &t);
      printf("Block %d : xmin=0.0 width=%e xmax=%f\n",n, w, nblock*w);
      ++n;
    }
    for(int ex=e_min; ex<e_max;++ex,++n) {
      double x = scalbn(1.0, ex);
      double w = scalbn(1.0, ex-wman);
      printf("Block %d : xmin=%f width=%e xmax=%f\n", n, x, w, x+nblock*w);
      mdg_table_add_intervals(Func, nblock, x, w, &t);
    }
  } else {
   mdg_table_init(&t, norder, ndiv, name, allow_carryup);

   double w = (xmax-xmin)/ndiv;
   mdg_table_add_intervals(Func, ndiv, xmin, w, &t);
 }

 if (output_floating) {
   printf("Floating output\n");
   mdg_table_integer_conversion_float(nbit, &t);
 } else {
   printf("Fixed output\n");
   mdg_table_integer_conversion(nbit, &t);
 }

 mdg_table_plot(10, 1, Func, &t);

 mdg_table_output_c(name, &t, cmdline);
 mdg_table_output_verilog(name, sign_mode, &t, cmdline);

}
