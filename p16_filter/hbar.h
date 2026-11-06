struct w
{
    float zn1;
    float zn2;
};

static struct w w;

float hbar(float u)
{
    float wCurrent = -1.609375*w.zn1 - 0.658441*w.zn2 + u;
    float y = wCurrent + 2*w.zn1 + w.zn2;
    w.zn2 = w.zn1;
    w.zn1 = wCurrent;
    return y;
}
