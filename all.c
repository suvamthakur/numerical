// newtons_forward.c
#include <stdio.h>
#include <math.h>

void main() {
    int n, i, j;
    float x[20], y[20][20], xp, yp = 0, p;
    
    printf("\nNewton's Forward Interpolation\n");
    printf("Enter number of data points: ");
    scanf("%d", &n);
    
    printf("Enter data points:\n");
    for(i = 0; i < n; i++) {
        printf("x[%d] = ", i);
        scanf("%f", &x[i]);
        printf("y[%d] = ", i);
        scanf("%f", &y[i][0]);
    }
    
    // Calculate forward differences
    for(j = 1; j < n; j++) {
        for(i = 0; i < n-j; i++) {
            y[i][j] = y[i+1][j-1] - y[i][j-1];
        }
    }
    
    printf("\nEnter interpolation point: ");
    scanf("%f", &xp);
    
    // Calculating interpolated value
    p = (xp - x[0]) / (x[1] - x[0]);
    
    yp = y[0][0];
    float term = p;
    int fact = 1;
    
    for(j = 1; j < n; j++) {
        yp += (term * y[0][j]) / fact;
        term *= (p - j);
        fact *= (j + 1);
    }
    
    printf("\nInterpolated value at x = %.4f is %.4f\n", xp, yp);
}



// newtons_backward.c
#include <stdio.h>
#include <math.h>

void main() {
    int n, i, j;
    float x[20], y[20][20], xp, yp = 0, p;
    
    printf("\nNewton's Backward Interpolation\n");
    printf("Enter number of data points: ");
    scanf("%d", &n);
    
    printf("Enter data points:\n");
    for(i = 0; i < n; i++) {
        printf("x[%d] = ", i);
        scanf("%f", &x[i]);
        printf("y[%d] = ", i);
        scanf("%f", &y[i][0]);
    }
    
    // Calculate backward differences
    for(j = 1; j < n; j++) {
        for(i = n-1; i >= j; i--) {
            y[i][j] = y[i][j-1] - y[i-1][j-1];
        }
    }
    
    printf("\nEnter interpolation point: ");
    scanf("%f", &xp);
    
    p = (xp - x[n-1]) / (x[1] - x[0]);
    
    yp = y[n-1][0];
    float term = p;
    int fact = 1;
    
    for(j = 1; j < n; j++) {
        yp += (term * y[n-1][j]) / fact;
        term *= (p + j);
        fact *= (j + 1);
    }
    
    printf("\nInterpolated value at x = %.4f is %.4f\n", xp, yp);
}


// divided_difference.c
#include <stdio.h>
#include <math.h>

void main() {
    int n, i, j;
    float x[20], y[20][20], xp, yp, p;
    
    printf("\nNewton's Divided Difference Interpolation\n");
    printf("Enter number of data points: ");
    scanf("%d", &n);
    
    printf("Enter data points:\n");
    for(i = 0; i < n; i++) {
        printf("x[%d] = ", i);
        scanf("%f", &x[i]);
        printf("y[%d] = ", i);
        scanf("%f", &y[i][0]);
    }
    
    // Calculate divided differences
    for(j = 1; j < n; j++) {
        for(i = 0; i < n-j; i++) {
            y[i][j] = (y[i+1][j-1] - y[i][j-1]) / (x[i+j] - x[i]);
        }
    }
    
    printf("\nEnter interpolation point: ");
    scanf("%f", &xp);
    
    // Initial term
    yp = y[0][0];
    
    // Apply formula
    for(j = 1; j < n; j++) {
        p = 1;
        for(i = 0; i < j; i++) {
            p = p * (xp - x[i]);
        }
        yp = yp + p * y[0][j];
    }
    
    printf("\nInterpolated value at x = %.4f is %.4f\n", xp, yp);
}



// lagrange_interpolation.c
#include <stdio.h>

void main() {
    float x[20], y[20], xp, yp = 0, p;
    int i, j, n;
    
    printf("\nLagrange's Interpolation\n");
    printf("Enter number of data points: ");
    scanf("%d", &n);
    
    printf("Enter data points:\n");
    for(i = 0; i < n; i++) {
        printf("x[%d] = ", i);
        scanf("%f", &x[i]);
        printf("y[%d] = ", i);
        scanf("%f", &y[i]);
    }
    
    printf("\nEnter interpolation point: ");
    scanf("%f", &xp);
    
    for(i = 0; i < n; i++) {
        p = 1;
        for(j = 0; j < n; j++) {
            if(i != j) {
                p = p * (xp - x[j])/(x[i] - x[j]);
            }
        }
        yp = yp + p * y[i];
    }
    
    printf("\nInterpolated value at x = %.4f is %.4f\n", xp, yp);
}


// numerical_differentiation.c
#include <stdio.h>

void forward_diff();
void backward_diff();

void main() {
    int choice;
    
    printf("\nNumerical Differentiation\n");
    printf("1. Forward Differentiation\n");
    printf("2. Backward Differentiation\n");
    printf("Enter choice (1/2): ");
    scanf("%d", &choice);
    
    if(choice == 1)
        forward_diff();
    else if(choice == 2)
        backward_diff();
    else
        printf("Invalid choice!\n");
}

void forward_diff() {
    float x[20], y[20], h, x_val, derivative;
    int n, i;
    
    printf("\nForward Differentiation\n");
    printf("Enter number of points: ");
    scanf("%d", &n);
    
    printf("Enter data points:\n");
    for(i = 0; i < n; i++) {
        printf("x[%d] = ", i);
        scanf("%f", &x[i]);
        printf("y[%d] = ", i);
        scanf("%f", &y[i]);
    }
    
    printf("Enter the point at which derivative is required: ");
    scanf("%f", &x_val);
    
    h = x[1] - x[0];
    
    // Find the position of x_val in the array
    for(i = 0; i < n; i++) {
        if(x[i] == x_val) {
            if(i == n-1) {
                printf("Cannot calculate forward derivative at last point!\n");
                return;
            }
            derivative = (y[i+1] - y[i])/h;
            printf("Forward derivative at x = %.4f is %.4f\n", x_val, derivative);
            return;
        }
    }
    printf("Point not found in data!\n");
}

void backward_diff() {
    float x[20], y[20], h, x_val, derivative;
    int n, i;
    
    printf("\nBackward Differentiation\n");
    printf("Enter number of points: ");
    scanf("%d", &n);
    
    printf("Enter data points:\n");
    for(i = 0; i < n; i++) {
        printf("x[%d] = ", i);
        scanf("%f", &x[i]);
        printf("y[%d] = ", i);
        scanf("%f", &y[i]);
    }
    
    printf("Enter the point at which derivative is required: ");
    scanf("%f", &x_val);
    
    h = x[1] - x[0];
    
    // Find the position of x_val in the array
    for(i = 0; i < n; i++) {
        if(x[i] == x_val) {
            if(i == 0) {
                printf("Cannot calculate backward derivative at first point!\n");
                return;
            }
            derivative = (y[i] - y[i-1])/h;
            printf("Backward derivative at x = %.4f is %.4f\n", x_val, derivative);
            return;
        }
    }
    printf("Point not found in data!\n");
}


// euler_method.c
#include <stdio.h>
#include <math.h>

// Define the differential equation dy/dx = f(x,y)
float f(float x, float y) {
    return x + y; // Example: dy/dx = x + y
}

void main() {
    float x0, y0, h, xn, x, y;
    int i, n;
    
    printf("\nEuler's Method\n");
    printf("Enter initial values:\n");
    printf("x0 = ");
    scanf("%f", &x0);
    printf("y0 = ");
    scanf("%f", &y0);
    
    printf("Enter calculation point xn = ");
    scanf("%f", &xn);
    printf("Enter step size h = ");
    scanf("%f", &h);
    
    n = (int)((xn - x0)/h);
    
    x = x0;
    y = y0;
    
    printf("\nx\t\ty\n");
    printf("------------------------\n");
    printf("%.4f\t%.4f\n", x, y);
    
    for(i = 1; i <= n; i++) {
        y = y + h * f(x, y);
        x = x + h;
        printf("%.4f\t%.4f\n", x, y);
    }
}



// taylor_series.c
#include <stdio.h>
#include <math.h>

// Define the function and its derivatives
float f(float x) { return sin(x); }
float f1(float x) { return cos(x); }
float f2(float x) { return -sin(x); }
float f3(float x) { return -cos(x); }

void main() {
    float x0, x, h;
    int n, i;
    
    printf("\nTaylor Series Method\n");
    printf("Enter x0 around which to expand: ");
    scanf("%f", &x0);
    printf("Enter point x at which to evaluate: ");
    scanf("%f", &x);
    printf("Enter number of terms (1-4): ");
    scanf("%d", &n);
    
    h = x - x0;
    float result = f(x0);
    
    if(n >= 2) result += h * f1(x0);
    if(n >= 3) result += (h*h/2) * f2(x0);
    if(n >= 4) result += (h*h*h/6) * f3(x0);
    
    printf("\nTaylor series approximation = %.6f\n", result);
    printf("Actual value = %.6f\n", f(x));
    printf("Error = %.6f\n", fabs(f(x) - result));
}



// runge_kutta.c
#include <stdio.h>
#include <math.h>

float f(float x, float y) {
    return x + y; // Example: dy/dx = x + y
}

void main() {
    float x0, y0, h, xn, k1, k2, k3, k4, x, y;
    int i, n;
    
    printf("\nRunge-Kutta 4th Order Method\n");
    printf("Enter initial values:\n");
    printf("x0 = ");
    scanf("%f", &x0);
    printf("y0 = ");
    scanf("%f", &y0);
    
    printf("Enter calculation point xn = ");
    scanf("%f", &xn);
    printf("Enter step size h = ");
    scanf("%f", &h);
    
    n = (int)((xn - x0)/h);
    
    x = x0;
    y = y0;
    
    printf("\nx\t\ty\n");
    printf("------------------------\n");
    printf("%.4f\t%.4f\n", x, y);
    
    for(i = 1; i <= n; i++) {
        k1 = h * f(x, y);
        k2 = h * f(x + h/2, y + k1/2);
        k3 = h * f(x + h/2, y + k2/2);
        k4 = h * f(x + h, y + k3);
        
        y = y + (k1 + 2*k2 + 2*k3 + k4)/6;
        x = x + h;
        
        printf("%.4f\t%.4f\n", x, y);
    }
}



// newton_nonlinear.c
#include <stdio.h>
#include <math.h>

#define EPSILON 0.0001
#define MAX_ITER 100

// Define the system of equations
float f1(float x, float y) {
    return x*x + y*y - 25;    // Example: x² + y² = 25
}

float f2(float x, float y) {
    return x*y - 9;           // Example: xy = 9
}

// Define the Jacobian elements
float df1_dx(float x, float y) { return 2*x; }
float df1_dy(float x, float y) { return 2*y; }
float df2_dx(float x, float y) { return y; }
float df2_dy(float x, float y) { return x; }

void main() {
    float x0, y0, x1, y1, det, dx, dy;
    int iter = 0;
    
    printf("\nNewton's Method for Nonlinear Equations\n");
    printf("Enter initial guess:\n");
    printf("x0 = ");
    scanf("%f", &x0);
    printf("y0 = ");
    scanf("%f", &y0);
    
    do {
        // Calculate function values
        float f1_val = f1(x0, y0);
        float f2_val = f2(x0, y0);
        
        // Calculate Jacobian
        float j11 = df1_dx(x0, y0);
        float j12 = df1_dy(x0, y0);
        float j21 = df2_dx(x0, y0);
        float j22 = df2_dy(x0, y0);
        
        // Calculate determinant
        det = j11*j22 - j12*j21;
        
        // Calculate changes
        dx = (-f1_val*j22 + f2_val*j12)/det;
        dy = (-f2_val*j11 + f1_val*j21)/det;
        
        // Update values
        x1 = x0 + dx;
        y1 = y0 + dy;
        
        printf("Iteration %d: x = %.6f, y = %.6f\n", iter+1, x1, y1);
        
        if(fabs(dx) < EPSILON && fabs(dy) < EPSILON) {
            printf("\nSolution found:\nx = %.6f\ny = %.6f\n", x1, y1);
            return;
        }
        
        x0 = x1;
        y0 = y1;
        iter++;
        
    } while(iter < MAX_ITER);
    
    printf("Solution not found within %d iterations\n", MAX_ITER);
}


// gauss_thomas.c
#include <stdio.h>
#include <math.h>

void main() {
    int n, i;
    float a[20], b[20], c[20], d[20], x[20], p[20], q[20];
    
    printf("\nGauss Thomas Method for Tri-diagonal System\n");
    printf("Enter size of matrix: ");
    scanf("%d", &n);
    
    printf("\nEnter elements of tri-diagonal matrix:\n");
    printf("Enter lower diagonal elements (a):\n");
    for(i = 1; i < n; i++) {
        printf("a[%d] = ", i);
        scanf("%f", &a[i]);
    }
    
    printf("\nEnter main diagonal elements (b):\n");
    for(i = 0; i < n; i++) {
        printf("b[%d] = ", i);
        scanf("%f", &b[i]);
    }
    
    printf("\nEnter upper diagonal elements (c):\n");
    for(i = 0; i < n-1; i++) {
        printf("c[%d] = ", i);
        scanf("%f", &c[i]);
    }
    
    printf("\nEnter right hand side elements (d):\n");
    for(i = 0; i < n; i++) {
        printf("d[%d] = ", i);
        scanf("%f", &d[i]);
    }
    
    // Forward elimination
    p[0] = b[0];
    q[0] = d[0];
    
    for(i = 1; i < n; i++) {
        p[i] = b[i] - (a[i]*c[i-1])/p[i-1];
        q[i] = d[i] - (a[i]*q[i-1])/p[i-1];
    }
    
    // Back substitution
    x[n-1] = q[n-1]/p[n-1];
    
    for(i = n-2; i >= 0; i--) {
        x[i] = (q[i] - c[i]*x[i+1])/p[i];
    }
    
    printf("\nSolution:\n");
    for(i = 0; i < n; i++) {
        printf("x[%d] = %.4f\n", i, x[i]);
    }
}



// cubic_spline.c
#include <stdio.h>
#include <stdlib.h>

void main() {
    int n, i;
    float *x, *y, *h, *a, *b, *c, *d, *l, *u, *z, *alpha;
    float xp, yp;
    
    printf("\nCubic Spline Interpolation\n");
    printf("Enter number of data points: ");
    scanf("%d", &n);
    
    // Allocate memory
    x = (float*)malloc(n * sizeof(float));
    y = (float*)malloc(n * sizeof(float));
    h = (float*)malloc((n-1) * sizeof(float));
    a = (float*)malloc(n * sizeof(float));
    b = (float*)malloc(n * sizeof(float));
    c = (float*)malloc(n * sizeof(float));
    d = (float*)malloc(n * sizeof(float));
    l = (float*)malloc(n * sizeof(float));
    u = (float*)malloc(n * sizeof(float));
    z = (float*)malloc(n * sizeof(float));
    alpha = (float*)malloc(n * sizeof(float));
    
    printf("Enter data points:\n");
    for(i = 0; i < n; i++) {
        printf("x[%d] = ", i);
        scanf("%f", &x[i]);
        printf("y[%d] = ", i);
        scanf("%f", &y[i]);
    }
    
    // Step 1: Calculate h
    for(i = 0; i < n-1; i++) {
        h[i] = x[i+1] - x[i];
    }
    
    // Step 2: Calculate alpha
    for(i = 1; i < n-1; i++) {
        alpha[i] = 3/h[i]*(y[i+1]-y[i]) - 3/h[i-1]*(y[i]-y[i-1]);
    }
    
    // Step 3: Calculate l, u, z
    l[0] = 1;
    u[0] = 0;
    z[0] = 0;
    
    for(i = 1; i < n-1; i++) {
        l[i] = 2*(x[i+1]-x[i-1]) - h[i-1]*u[i-1];
        u[i] = h[i]/l[i];
        z[i] = (alpha[i] - h[i-1]*z[i-1])/l[i];
    }
    
    l[n-1] = 1;
    z[n-1] = 0;
    c[n-1] = 0;
    
    // Step 4: Calculate coefficients
    for(i = n-2; i >= 0; i--) {
        c[i] = z[i] - u[i]*c[i+1];
        b[i] = (y[i+1]-y[i])/h[i] - h[i]*(c[i+1]+2*c[i])/3;
        d[i] = (c[i+1]-c[i])/(3*h[i]);
        a[i] = y[i];
    }
    
    printf("\nEnter point for interpolation: ");
    scanf("%f", &xp);
    
    // Find appropriate interval
    for(i = 0; i < n-1; i++) {
        if(xp >= x[i] && xp <= x[i+1]) {
            float dx = xp - x[i];
            yp = a[i] + b[i]*dx + c[i]*dx*dx + d[i]*dx*dx*dx;
            printf("Interpolated value at x = %.4f is %.4f\n", xp, yp);
            break;
        }
    }
    
    // Free allocated memory
    free(x); free(y); free(h); free(a); free(b);
    free(c); free(d); free(l); free(u); free(z); free(alpha);
}



// romberg_integration.c
#include <stdio.h>
#include <math.h>

// Define the function to integrate
float f(float x) {
    return x*x;  // Example: f(x) = x²
}

void main() {
    int i, j, k, n;
    float a, b, h, R[20][20], sum;
    
    printf("\nRomberg Integration\n");
    printf("Enter lower limit a: ");
    scanf("%f", &a);
    printf("Enter upper limit b: ");
    scanf("%f", &b);
    printf("Enter number of iterations: ");
    scanf("%d", &n);
    
    // First approximation R(0,0)
    h = b - a;
    R[0][0] = h * (f(a) + f(b)) / 2;
    
    printf("\nRomberg Integration Table:\n");
    printf("R[0][0] = %.6f\n", R[0][0]);
    
    // Calculate other values
    for(i = 1; i < n; i++) {
        h = h/2;
        sum = 0;
        
        // Calculate R(i,0)
        for(k = 1; k <= pow(2,i)-1; k += 2) {
            sum += f(a + k*h);
        }
        R[i][0] = R[i-1][0]/2 + h*sum;
        
        // Calculate R(i,j)
        for(j = 1; j <= i; j++) {
            R[i][j] = R[i][j-1] + (R[i][j-1] - R[i-1][j-1])/(pow(4,j)-1);
        }
        
        // Print row
        printf("Row %d: ", i);
        for(j = 0; j <= i; j++) {
            printf("%.6f ", R[i][j]);
        }
        printf("\n");
    }
    
    printf("\nFinal Result = %.6f\n", R[n-1][n-1]);
}



// finite_difference.c
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void main() {
    int n, i, j, max_iter = 100;
    float h, k, error = 0.0001;
    float *x, *y, *temp;
    
    printf("\nFinite Difference Method for Boundary Value Problems\n");
    printf("Solving: y'' = f(x,y) with boundary conditions\n");
    printf("Enter number of internal points: ");
    scanf("%d", &n);
    
    // Allocate memory
    x = (float*)malloc((n+2) * sizeof(float));
    y = (float*)malloc((n+2) * sizeof(float));
    temp = (float*)malloc((n+2) * sizeof(float));
    
    printf("Enter boundary conditions:\n");
    printf("y(0) = ");
    scanf("%f", &y[0]);
    printf("y(1) = ");
    scanf("%f", &y[n+1]);
    
    // Calculate step size
    h = 1.0/(n+1);
    
    // Initialize x values
    for(i = 0; i <= n+1; i++) {
        x[i] = i*h;
    }
    
    // Initialize y values
    for(i = 1; i <= n; i++) {
        y[i] = 0;  // Initial guess
    }
    
    // Iterative solution
    int iter = 0;
    float max_diff;
    
    do {
        // Save current values
        for(i = 0; i <= n+1; i++) {
            temp[i] = y[i];
        }
        
        // Update y values
        for(i = 1; i <= n; i++) {
            // Using central difference formula
            y[i] = ((y[i-1] + y[i+1])/2) + (h*h/2)*x[i];  // Example: y'' = x
        }
        
        // Check convergence
        max_diff = 0;
        for(i = 1; i <= n; i++) {
            if(fabs(y[i] - temp[i]) > max_diff) {
                max_diff = fabs(y[i] - temp[i]);
            }
        }
        
        iter++;
        
    } while(max_diff > error && iter < max_iter);
    
    printf("\nSolution:\n");
    printf("x\t\ty\n");
    printf("-------------------\n");
    for(i = 0; i <= n+1; i++) {
        printf("%.4f\t%.4f\n", x[i], y[i]);
    }
    
    printf("\nNumber of iterations: %d\n", iter);
    
    // Free allocated memory
    free(x);
    free(y);
    free(temp);
}
