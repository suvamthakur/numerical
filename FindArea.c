// simpson,  trapezoidal
// F(x) = 3x^2 + 2x + 5;

#include <stdio.h>

float simpson(int, int, int);
float trapezoidal(int, int, int);
float function(float);

int main() {
  int ll, ul, n;
  float areaS, areaT;

  printf("Enter lower limit and upper limit: ");
  scanf("%d %d", &ll, &ul);
  printf("Enter no of division: ");
  scanf("%d", &n);

  areaS = simpson(ll, ul, n);
  printf("Simpson's 1/3 method: %g", areaS);

  areaT = trapezoidal(ll, ul, n);
  printf("\nTrapezoidal method: %g", areaT);

  return 0;
}

float function(float x) {
  return 3 * x * x + 2 * x + 5;
}

float simpson(int ll, int ul, int n) {
  float y[15], h, result, x, even = 0, odd = 0;
  int i = 0;

  h = (float)(ul - ll) / n;
  x = ll;
  while (x <= ul) {
    y[i] = function(x);
    x = x + h;
    i++;
  }

  for (int i = 1; i < ul; i++) {
    if (i % 2 == 0) {
      even = even + y[i];
    } else {
      odd = odd + y[i];
    }
  }

  result = (h / 3) * ((y[0] + y[ul]) + 2 * odd + 4 * even);
  return result;
}

float trapezoidal(int ll, int ul, int n) {
  float y[15], x, h, sum = 0, result;
  int i = 0;

  h = (float)(ul - ll) / n;
  x = ll;
  while (x <= ul) {
    y[i] = function(x);
    x = x + h;
    i++;
  }

  for (int i = 1; i < ul; i++) {
    sum = sum + y[i];
  }

  result = (h / 2) * ((y[0] + y[ul]) + 2 * sum);
  return result;
}