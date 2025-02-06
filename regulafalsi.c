#include <stdio.h>
#include<math.h>

#define E 2.718281828459045

float regulafalsi(float, float);
float f(float);

int main() {
    float xn, xp;
    
    printf("Enter negative and positive no: ");
    scanf("%f %f", &xn, &xp);
    
    float ans = regulafalsi(xn, xp);
    printf("\nRoot is: %f", ans);
    return 0;
}

float f(float x) {
    return x*x*x - 2*x - 5;
}

float regulafalsi(float xn, float xp) {
    float xpre, xnext;
    float c = (xn*f(xp) - xp*f(xn))/(f(xp) - f(xn));
    do {
        float res = f(c);
        if(res < 0) {
            xn = c;
        }
        else {
            xp = c;
        }
        printf("\n%f is c", c);
        xpre = round(c*1000);
        c = (xn*f(xp) - xp*f(xn))/(f(xp) - f(xn));
        xnext = round(c*1000);
    } while(xpre != xnext);
    return c;
}
