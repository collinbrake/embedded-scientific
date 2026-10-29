#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

/* Include the user's implementation exactly as is */
#include "hbar.h"

/*
 * ==============================================================================
 * Problem 16 Filter Test Harness
 *
 * This test suite evaluates the implementation in hbar.h against the
 * requirements specified in Problem 16:
 *   - Transfer function: H(z) with 2 zeros at z = -1 and poles at
 *     z = 0.8046875 +/- 0.1044925 i
 *   - Direct Form II difference equation structure
 *   - Fixed-point coefficient approximation [16, -6] (scale factor = 64)
 * ==============================================================================
 */

/* Theoretical Direct Form II reference with exact [16, -6] coefficients:
 * b = [1.0, 2.0, 1.0]
 * a1 = -103/64 = -1.609375
 * a2 = 42/64 = 0.65625
 * w[k] = u[k] + (103/64)*w[k-1] - (42/64)*w[k-2]
 * y[k] = w[k] + 2*w[k-1] + w[k-2]
 */
static void theoretical_hbar_step(double u_real, double *w1, double *w2, double *y_out, double *w_curr_out)
{
    double w_curr = u_real + (103.0 / 64.0) * (*w1) - (42.0 / 64.0) * (*w2);
    *y_out = w_curr + 2.0 * (*w1) + (*w2);
    *w_curr_out = w_curr;
    *w2 = *w1;
    *w1 = w_curr;
}

/* Theoretical infinite-precision reference H(z) */
static void theoretical_h_exact_step(double u_real, double *w1, double *w2, double *y_out)
{
    /* H(z) denominator: 1 - 1.609375*z^-1 + 0.6584406552125*z^-2 */
    double w_curr = u_real + 1.609375 * (*w1) - 0.6584406552125 * (*w2);
    *y_out = w_curr + 2.0 * (*w1) + (*w2);
    *w2 = *w1;
    *w1 = w_curr;
}

static void print_separator(char c, int len)
{
    for (int i = 0; i < len; i++) putchar(c);
    putchar('\n');
}

/* -------------------------------------------------------------------------- */
/* Test 1: State Retention and Update Check                                   */
/* -------------------------------------------------------------------------- */
static void test_state_updates(void)
{
    printf("\n");
    print_separator('=', 78);
    printf("TEST 1: State Retention & Delay Line Update Check\n");
    print_separator('=', 78);
    printf("Checking if internal states w.zn1 (w[k-1]) and w.zn2 (w[k-2]) are updated\n");
    printf("when hbar(u) is executed:\n\n");

    /* Reset states manually through access to static struct w */
    w.zn1 = 0;
    w.zn2 = 0;

    printf("Initial states: w.zn1 = %d, w.zn2 = %d\n", w.zn1, w.zn2);

    /* Step 1: Feed input u = 64 (representing 1.0 in [16, -6]) */
    signed short u0 = 64;
    signed short y0 = hbar(u0);
    printf("After hbar(%d): returned y = %d | States: w.zn1 = %d, w.zn2 = %d\n",
           u0, y0, w.zn1, w.zn2);

    /* Step 2: Feed input u = 0 */
    signed short u1 = 0;
    signed short y1 = hbar(u1);
    printf("After hbar(%d):  returned y = %d | States: w.zn1 = %d, w.zn2 = %d\n",
           u1, y1, w.zn1, w.zn2);

    /* Step 3: Feed input u = 0 again */
    signed short u2 = 0;
    signed short y2 = hbar(u2);
    printf("After hbar(%d):  returned y = %d | States: w.zn1 = %d, w.zn2 = %d\n\n",
           u2, y2, w.zn1, w.zn2);

    if (w.zn1 == 0 && w.zn2 == 0) {
        printf(">> DIAGNOSTIC RESULT: [BUG DETECTED] Delay states are NEVER updated!\n");
        printf("   In hbar.h, w.zn1 and w.zn2 are read but never written to.\n");
        printf("   Without updating (w.zn2 = w.zn1; w.zn1 = wCurrent;), the filter\n");
        printf("   acts as a static memoryless gain: y[k] = 64 * u[k].\n");
    } else {
        printf(">> DIAGNOSTIC RESULT: Delay states are being updated.\n");
    }
}

/* -------------------------------------------------------------------------- */
/* Test 2: Impulse Response Evaluation                                        */
/* -------------------------------------------------------------------------- */
static void test_impulse_response(void)
{
    printf("\n");
    print_separator('=', 78);
    printf("TEST 2: Impulse Response Comparison (u[0] = 64 = 1.0 in Q6, u[k>=1] = 0)\n");
    print_separator('=', 78);
    printf("%-4s | %-6s | %-12s | %-12s | %-14s | %-14s\n",
           "k", "u[k]", "w.zn1", "w.zn2", "hbar() Actual", "Hbar Expected");
    print_separator('-', 78);

    w.zn1 = 0;
    w.zn2 = 0;

    double theo_w1 = 0.0, theo_w2 = 0.0;

    for (int k = 0; k < 12; k++) {
        signed short u_k = (k == 0) ? 64 : 0;
        double u_real = (k == 0) ? 1.0 : 0.0;

        signed short zn1_before = w.zn1;
        signed short zn2_before = w.zn2;
        signed short y_actual = hbar(u_k);

        double y_theo, w_theo;
        theoretical_hbar_step(u_real, &theo_w1, &theo_w2, &y_theo, &w_theo);

        printf("%-4d | %-6d | %-12d | %-12d | %-14d | %-14.4f\n",
               k, u_k, zn1_before, zn2_before, y_actual, y_theo * 64.0);
    }

    printf("\n>> OBSERVATION:\n");
    printf("   Actual hbar() outputs [4096, 0, 0, 0, ...] (which is 64*64).\n");
    printf("   Expected Hbar(z) response peaks around k=4 at ~518 (real ~8.09) and\n");
    printf("   decays back toward 0 (stable pole at |z| ~ 0.81).\n");
}

/* -------------------------------------------------------------------------- */
/* Test 3: Unit Step Response and DC Gain Evaluation                          */
/* -------------------------------------------------------------------------- */
static void test_step_response(void)
{
    printf("\n");
    print_separator('=', 78);
    printf("TEST 3: Step Response & DC Gain Evaluation (u[k] = 64 = 1.0 in Q6)\n");
    print_separator('=', 78);

    w.zn1 = 0;
    w.zn2 = 0;

    double theo_w1 = 0.0, theo_w2 = 0.0;
    double y_theo = 0.0, w_theo = 0.0;

    printf("Theoretical DC Gain Calculation:\n");
    printf("  H(1) (infinite precision) = 4 / (1 - 1.609375 + 0.658440655) = 81.5234\n");
    printf("  Hbar(1) (fixed-point [16,-6]) = 4 / (1 - 103/64 + 42/64) = 256/3 = 85.3333\n");
    printf("  In Q6 representation (scaled by 64): 85.3333 * 64 = 5461.33\n\n");

    printf("%-4s | %-6s | %-14s | %-16s | %-16s\n",
           "k", "u[k]", "hbar() Actual", "Hbar Expected Q6", "Hbar Expected Real");
    print_separator('-', 78);

    for (int k = 0; k < 15; k++) {
        signed short y_actual = hbar(64);
        theoretical_hbar_step(1.0, &theo_w1, &theo_w2, &y_theo, &w_theo);

        printf("%-4d | %-6d | %-14d | %-16.1f | %-16.4f\n",
               k, 64, y_actual, y_theo * 64.0, y_theo);
    }

    printf("\n>> OBSERVATION:\n");
    printf("   Actual hbar() consistently outputs 4096 (static 64 * 64).\n");
    printf("   A properly working Direct Form II filter converges to ~5461 (~85.33).\n");
}

/* -------------------------------------------------------------------------- */
/* Test 4: Simulation of Updating States With Current Coefficients (Overflow) */
/* -------------------------------------------------------------------------- */
static void test_overflow_and_scaling_analysis(void)
{
    printf("\n");
    print_separator('=', 78);
    printf("TEST 4: Root Cause Analysis - Effect of Scaling & Sign on Feedback\n");
    print_separator('=', 78);
    printf("This test demonstrates what happens if states ARE updated using the exact\n");
    printf("formulas in hbar.h without the /64 division: wCurrent = -103*w1 - 42*w2 + u\n\n");

    signed short mock_zn1 = 0;
    signed short mock_zn2 = 0;

    printf("%-4s | %-6s | %-12s | %-12s | %-14s | %-14s\n",
           "k", "u[k]", "zn1", "zn2", "wCurrent", "y");
    print_separator('-', 78);

    for (int k = 0; k < 5; k++) {
        signed short u = (k == 0) ? 1 : 0;
        signed short zn1 = mock_zn1;
        signed short zn2 = mock_zn2;

        /* Current logic in hbar.h */
        signed short wCurrent = -103 * zn1 - 42 * zn2 + u;
        signed short y = 64 * wCurrent + 128 * zn1 + 64 * zn2;

        printf("%-4d | %-6d | %-12d | %-12d | %-14d | %-14d\n",
               k, u, zn1, zn2, wCurrent, y);

        /* If delay states were shifted: */
        mock_zn2 = mock_zn1;
        mock_zn1 = wCurrent;
    }

    printf("\n>> CRITICAL FINDING:\n");
    printf("   Notice how at k=2, wCurrent reaches 10567, and at k=3 it overflows 16-bit\n");
    printf("   short (-32768..32767)! Why?\n");
    printf("   1) Feedback multipliers (-103 and -42) represent -103/64 and -42/64.\n");
    printf("      Without dividing by 64 (or >> 6), the feedback loop gain is 64x too large,\n");
    printf("      making the filter wildly unstable!\n");
    printf("   2) The sign of -a1 should be +103, not -103.\n");
    printf("   3) The output coefficients (64, 128, 64) also need division by 64.\n");
}

/* -------------------------------------------------------------------------- */
/* Test 5: Corrected Direct Form II Demonstration                              */
/* -------------------------------------------------------------------------- */
static void test_corrected_demonstration(void)
{
    printf("\n");
    print_separator('=', 78);
    printf("TEST 5: Verification of Corrected Direct Form II Fixed-Point Logic\n");
    print_separator('=', 78);
    printf("Implementing corrected logic:\n");
    printf("  wCurrent = u + (103 * zn1 - 42 * zn2) / 64;\n");
    printf("  y = (64 * wCurrent + 128 * zn1 + 64 * zn2) / 64; (= wCurrent + 2*zn1 + zn2)\n");
    printf("  zn2 = zn1;  zn1 = wCurrent;\n\n");

    int16_t c_zn1 = 0, c_zn2 = 0;
    double ref_w1 = 0.0, ref_w2 = 0.0;
    double ref_y = 0.0, ref_w = 0.0;

    printf("%-4s | %-6s | %-10s | %-10s | %-12s | %-12s | %-10s\n",
           "k", "u[k]", "Fixed Q6", "Fixed Real", "Float Ref", "Abs Error", "Status");
    print_separator('-', 78);

    double max_err = 0.0;
    for (int k = 0; k < 15; k++) {
        int16_t u = (k == 0) ? 64 : 0; /* Unit impulse = 1.0 * 64 */
        double u_float = (k == 0) ? 1.0 : 0.0;

        /* Corrected fixed-point step */
        int32_t feedback = (103L * (int32_t)c_zn1 - 42L * (int32_t)c_zn2) / 64L;
        int16_t wCurrent = (int16_t)((int32_t)u + feedback);
        int32_t y = (64L * (int32_t)wCurrent + 128L * (int32_t)c_zn1 + 64L * (int32_t)c_zn2) / 64L;
        c_zn2 = c_zn1;
        c_zn1 = wCurrent;

        /* Float reference */
        theoretical_hbar_step(u_float, &ref_w1, &ref_w2, &ref_y, &ref_w);

        double y_real = (double)y / 64.0;
        double err = fabs(y_real - ref_y);
        if (err > max_err) max_err = err;

        printf("%-4d | %-6d | %-10d | %-10.4f | %-12.4f | %-12.4f | %s\n",
               k, u, y, y_real, ref_y, err, (err < 0.25) ? "MATCH" : "DIFF");
    }

    printf("\n>> Max error between Q6 fixed-point and float reference: %.4f (within ~1 LSB)\n",
           max_err);
}

/* -------------------------------------------------------------------------- */
/* Main entry point                                                           */
/* -------------------------------------------------------------------------- */
int main(void)
{
    print_separator('#', 78);
    printf("         PROBLEM 16 DIGITAL FILTER TEST SUITE (hbar_test)\n");
    print_separator('#', 78);

    test_state_updates();
    test_impulse_response();
    test_step_response();
    test_overflow_and_scaling_analysis();
    test_corrected_demonstration();

    printf("\n");
    print_separator('#', 78);
    printf("                         SUMMARY & NEXT STEPS\n");
    print_separator('#', 78);
    printf("1. hbar.h currently lacks delay state updates: w.zn1 and w.zn2 are unchanged.\n");
    printf("2. The feedback sign on w.zn1 is -103 instead of +103 (matches sign error\n");
    printf("   in p16_filter_fixed_point.m where [1, +1.609, ...] was used instead of -1.609).\n");
    printf("3. Fixed-point multipliers (103, 42, 64, 128, 64) require division by 64 (>> 6)\n");
    printf("   to maintain proper signal scaling.\n");
    print_separator('#', 78);
    printf("\n");

    return 0;
}
