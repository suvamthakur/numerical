// Online C compiler to run C program online
#include <stdio.h>
#include<math.h>

float newton(float);
float dfunction(float);
float function(float); 

int main() {
    float xn, xp;
    
    printf("Enter negative and positive no: ");
    scanf("%f %f", &xn, &xp);
    
    float ans = newton((xn+xp)/2);
    printf("Root is: %f", ans);
    return 0;
}

float dfunction(float x) {
    return 3*x*x - 1;
}

float function(float x) {
    return (x*x*x - x - 7);
}

float newton(float x) {
    float xnext = x, xpre;
    
    do {
        xpre = xnext;
        float fx = function(xpre);
        float dfx = dfunction(xpre);
        
        xnext = xpre - (fx/dfx);
    } while(round(xpre*1000) != round(xnext*1000));
    return xpre;
}
