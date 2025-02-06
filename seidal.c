#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include <stdbool.h>

bool check(float a[][10], int n) {
    for(int i=0; i<n; i++) {
        float sum = 0;
        for(int j=0; j<n; j++) {
            if(i != j) {
                sum += fabs(a[i][j]);
            }
        }
        if(fabs(a[i][i]) < sum) return false;
    }
    return true;
}

void seidal(float a[][10],float b[10], int n) {
    float x[10], xn, sum, flag, error;
    float err = 0.0001;
    
    // x[] -> 0
    for (int i = 0; i < n; i++) {
        x[i] = 0;
    }
    printf("x[0] | x[1] | x[2] \n");
    do {
        error = 0;
        // Calculate
        for(int i=0; i<n; i++) {
            sum = b[i];
            for(int j=0; j<n; j++) {
                if(i != j) {
                    sum = sum - (a[i][j]*x[j]);
                }
            }
            xn = sum / a[i][i];
            error += fabs(xn - x[i]);
            x[i] = xn;
        }
        
        for(int i=0; i<n; i++) {
            printf("%10f", x[i]);
        }
        printf("\n");
    }while(error >= err);
    
    printf("Solution is: \n");
    for (int i = 0; i < n; i++) {
        printf("%8.5f ", x[i]);
    }
}


int main() {
    float a[10][10], b[10];
    int n;
    
    printf("Enter no of var: ");
    scanf("%d", &n);
    
    printf("Enter values row-wise: ");
    for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) {
            scanf("%f", &a[i][j]);
        }
    }
    
    printf("\nEnter right hand vector: ");
    for (int i = 0; i < n; i++) {
        scanf("%f", &b[i]);
    }
    
    if(check(a, n)) {
        seidal(a, b, n);
    }
    else {
        printf("The coefficient matrix is not diagonally dominant");
    }

    return 0;
}
