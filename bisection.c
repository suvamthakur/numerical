// Online C compiler to run C program online
#include <stdio.h>
#include<math.h>

float bisection(float, float);
float function(float);

int main() {
    float xn, xp;
    
    printf("Enter negative and positive no: ");
    scanf("%f %f", &xn, &xp);
    
    float ans = bisection(xn, xp);
    printf("Root is: %f", ans);
    return 0;
}

float function(float x) {
    return (x*x*x - x - 7);
}

float bisection(float xn, float xp) {
    float xpre, xnext;
    float mid = (xn+xp)/2;
    do {
        float res = function(mid);
        if(res < 0) {
            xn = mid;
        }
        else {
            xp = mid;
        }
        xpre = round(mid*1000);
        mid = (xn+xp)/2;
        xnext = round(mid*1000);
    } while(xpre != xnext);
    return mid;
}
