#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "h20.h"

static void test_step_response(void)
{
    printf("\n");
    printf("TEST: Step Response (u[k] = 1.0)\n");

    for (int k = 0; k < 25; k++) {
        double y = h20(1.0);

        printf("%-4d | %-6.1f | %-14f |\n",
               k, 1.0, y);
    }
}

int main(void)
{
    test_step_response();

    return 0;
}
