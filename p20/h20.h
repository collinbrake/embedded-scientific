struct x
{
    double x1;
    double x2;
    double x3;
    double x4;
};

static struct x xkp1;

double h20(double uk)
{
    // x[k] is x[k+1] from last iteration
    struct x xk;
    xk.x1 = xkp1.x1;
    xk.x2 = xkp1.x2;
    xk.x3 = xkp1.x3;
    xk.x4 = xkp1.x4;

    // state update for Subsystem 1 (Second-order block 1: complex conjugate pair)
    xkp1.x1 = xk.x2;
    xkp1.x2 = -0.36717889 * xk.x1 + 1.13041431 * xk.x2 + uk;

    // state update for Subsystem 2 (Second-order block 2: real pole pair)
    xkp1.x3 = xk.x4;
    xkp1.x4 =  0.26898333 * xk.x3 + 5.36458569 * xk.x4 + uk;

    // output calculation (D = 0)
    return 0.04186632 * xk.x1 + 0.28720744 * xk.x2 + 0.83763351 * xk.x3 + 0.71279256 * xk.x4;
}
