#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>

#define N 84000000

int main(void) {
  std::chrono::time_point<std::chrono::system_clock> startTime;
  std::chrono::time_point<std::chrono::system_clock> endTime;
  startTime = std::chrono::system_clock::now();

  constexpr auto max_precision{std::numeric_limits<long double>::digits10 + 1}; 

  long double pi, exactpi;
  std::cout << "Computing approximation to pi using N = " << N << std::endl;

  pi = 0.0;

  for (int i = 1; i <= N; i++) {
    pi = pi + 1.0/( 1.0 + pow((((long double) i) - 0.5) / ((long double) N), 2.0));
  }

  pi = pi * 4.0 / ((long double) N);

  exactpi = 4.0 * atan(1.0);

  long double error = fabs(100.0 * (pi- exactpi) / exactpi);
  
  std::cout << "pi = " << std::setprecision(max_precision) << pi << ", error = " << error << "%" << std::endl;

  std::cout << "Time taken = " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - startTime).count() << "ms" << std::endl;
  return 0;
}
