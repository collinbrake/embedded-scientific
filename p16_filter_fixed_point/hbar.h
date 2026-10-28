#define PATTERN 0xB97AE000;

struct w
{
    signed short zn1;
    signed short zn2;
};

static struct w w;

signed short hbar(signed short u)
{
    signed short wCurrent = -103*w.zn1 - 42*w.zn2 + u;
    return 64*wCurrent + 128*w.zn1 + 64*w.zn2;
}
