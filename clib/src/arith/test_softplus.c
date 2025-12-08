#include <stdio.h>

#include "softplus.h"

int
main(int argc, char *argv[]) {
  for(int i=-1000;i<=1000;++i) {
    double x=i/100.0;
    printf("%le", x);
    printf(" %le", rial_softplus(x));
    printf(" %le", rial_softplus_deriv(x));
    for(int a=1;a<=8;++a) {
      printf(" %le", rial_softplus_parabolic(x,(double)a));
      printf(" %le", rial_softplus_parabolic_deriv(x,(double)a));
      printf(" %le", rial_softplus_quadratic(x,(double)a));
      printf(" %le", rial_softplus_quadratic_deriv(x,(double)a));
    }
    printf("\n");
  }
}
