#include <cstdlib>
#include <cmath>
#include <iostream>

#define N 840

int main(void) {
  double pi, exactpi;
  std::cout << "Computing approximation to pi using N = " << N;

  pi = 0.0;

  for (int i = 1; i <= N; i++) {
    pi = pi + 1.0/( 1.0 + pow( (((double) i) - 0.5) / ((double) N), 2.0));
  }

  pi = pi * 4.0 / ((double) N);

  exactpi = 4.0 * atan(1.0);

  double error = fabs(100.0 * (pi- exactpi) / exactpi);
  
  std::cout << "pi = " << pi << ", error = " << error << "%" << std::endl;

  return 0;
}
