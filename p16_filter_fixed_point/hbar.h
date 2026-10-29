struct w
{
    float zn1;
    float zn2;
};

static struct w w;

// This function expects u in the [16, -6] fixed point representation, so multipled by 64
//  it returns y in the same representation
float hbar(float u)
{
    float wCurrent = -103*w.zn1 - 42*w.zn2 + u;
    float y = wCurrent + 2*w.zn1 + w.zn2; // we are already working with 64*w[k], don't need 64*coefficients
    w.zn2 = w.zn1;
    w.zn1 = wCurrent;
    return y;
}
