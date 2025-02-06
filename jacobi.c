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

void jacobi(float a[][10],float b[10], int n) {
    float x[10], xn[10], sum, flag;
    float err = 0.0001;
    
    // x[] -> 0
    for (int i = 0; i < n; i++) {
        x[i] = 0;
    }
    printf("x[0] | x[1] | x[2] \n");
    do {
        // Calculate
        for(int i=0; i<n; i++) {
            sum = b[i];
            for(int j=0; j<n; j++) {
                if(i != j) {
                    sum = sum - (a[i][j]*x[j]);
                }
            }
            xn[i] = sum / a[i][i];
        }
        
        // check
        flag = 0;
        for(int i=0; i<n; i++) {
            if(fabs(xn[i] - x[i]) >= err) {
                flag = 1;
                break;
            }
        }
        
        if(flag == 1) {
            for(int i=0; i<n; i++) {
                printf("%10f", xn[i]);
                x[i] = xn[i];
            }
        }
        printf("\n");
    }while(flag == 1);
    
    printf("Solution is: \n");
    for (int i = 0; i < n; i++) {
        printf("%8.5f ", xn[i]);
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
        jacobi(a, b, n);
    }
    else {
        printf("The coefficient matrix is not diagonally dominant");
    }

    return 0;
}
