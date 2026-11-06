struct x
{
    double x1;
    double x2;
    double x3;
    double x4;
};

static struct x xkp1;

double hd(double uk)
{
    // x[k] is x[k+1] from last iteration
    struct x xk;
    xk.x1 = xkp1.x1;
    xk.x2 = xkp1.x2;
    xk.x3 = xkp1.x3;
    xk.x4 = xkp1.x4;

    // state update
    xkp1.x1 =   0.803278688525*xk.x1 + 0.163934426229*xk.x2 + 153.6818140292*uk; // block 1
    xkp1.x2 =  -0.163934426229*xk.x1 + 0.803278688525*xk.x2 + 222.6818121642*uk;
    xkp1.x3 =   0.818181818182*xk.x3 - 378.7173393922*uk;                        // block 2
    xkp1.x4 =   0.666666666667*xk.x4 + 190.6206957916*uk;                        // block 2

    // output
    return -0.001602281012*xk.x1 - 0.000607936458*xk.x2 - 0.001745779652*xk.x3 - 0.001457227803*xk.x4 + 0.0002483854942871*uk;
}
