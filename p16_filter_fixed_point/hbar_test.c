#define PATTERN 0xB97AE000;

unsigned short d;
unsigned long a;
signed long b;
float c;

int main(void) {
    a = PATTERN;
    b = PATTERN;
    c = PATTERN;
    d = PATTERN;
    printf("32-bit short: %u\n", d);
    printf("32-bit unsigned: %u\n", a);
    printf("32-bit signed: %d\n", b);
    printf("32-bit float: %f\n", c);
    return 0;
}
