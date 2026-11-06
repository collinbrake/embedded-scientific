#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "hd.h"

static void test_step_response(void)
{
    printf("\n");
    printf("TEST: Step Response & DC Gain Evaluation (u[k] = 1.0)\n");

    for (int k = 0; k < 115; k++) {
        float y = hd(1);

        printf("%-4d | %-6d | %-14f |\n",
               k, 1, y);
    }
}

int main(void)
{
    test_step_response();

    return 0;
}
