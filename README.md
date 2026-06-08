# Numerical Analysis Methods
Numerical Analysis methods are included:
1. Bisection
2. Regula-Falsi
3. Newton-Rapshon
4. NxN matrix inverse 
5. Gauss elimination
6. Gauss-Seidel
7. Numerical differentiation
8. Simpson
9. Trapez
10. Gregory-Newton interpolation 

# Differential Equations Runge Kutta4 Method
Numerical solution of linear ordinary differential equations with constant coefficients using the Runge-Kutta 4 method. 

## Compile & Run

```bash
# Numerical analysis
gcc Numerical-Analysis/numerical_analysis.c -o numerical -lm
./numerical

# Runge-Kutta
gcc Runge-Kutta4-Method/runge_kutta4.c -o rk4 -lm
./rk4
