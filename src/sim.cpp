
#include <RcppArmadillo.h>

using namespace Rcpp;


//' Filler function to allow package to build until real cpp functions are added.
//'
//' @noRd
//'
//[[Rcpp::export]]
int foo(int a, int b) {
    return a + b;
}
