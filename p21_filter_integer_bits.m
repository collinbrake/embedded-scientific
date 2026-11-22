% state space representation
A = [37 5. 0;...
        -5 37 0;...
        0 0 16] / 64;
B = [52; 78; 64] / 64;
C = [96 32 64] / 64;

% [10, -6] fixed point representation -- already there bc / 64 = / 2^6
% f = 6;
% A = round(A*2^f) / (2^f);
% B = round(B*2^f) / (2^f);
% C = round(C*2^f) / (2^f);

disp("============================= (A) ===============================")
% transfer function and eigenvalues
[num, den] = ss2tf(A, B, C, 0);
G = tf(num, den, -1)
eigenVals = eig(A)
disp("")

disp("============================= (B) ===============================")
disp("A more parallel implementation with real-valued multipliers")
disp("cannot be found since the third order system has a pair of")
disp("complex poles and is already block diagional")
disp("")

disp("============================= (D) ===============================")
[num, den] = ss2tf(A, B, [1 0 0], 0);
Gux1 = tf(num, den, -1)
[num, den] = ss2tf(A, B, [0 1 0], 0);
Gux2 = tf(num, den, -1)
[num, den] = ss2tf(A, B, [0 0 1], 0);
Gux3 = tf(num, den, -1)
disp("")

disp("============================= (E) ===============================")
figure(1)
hold on

[gux1, k] = impulse(Gux1);
plot(k, gux1)

[gux2, k] = impulse(Gux2);
plot(k, gux2)

[gux3, k] = impulse(Gux3);
plot(k, gux3)

disp("============================= (F) ===============================")
gux1_l1 = max(abs(gux1))
gux2_l1 = max(abs(gux2))
gux3_l1 = max(abs(gux3))

disp("============================= (G) ===============================")

u_linf = 5;

x1_linf = gux1_l1 * u_linf
x2_linf = gux2_l1 * u_linf
x3_linf = gux3_l1 * u_linf

% add 1 for sign bit
x1_nInt = ceil(log2(ceil(x1_linf))) + 1
x2_nInt = ceil(log2(ceil(x2_linf))) + 1
x3_nInt = ceil(log2(ceil(x3_linf))) + 1

disp("============================= (H) ===============================")

% special ones we need a transfer function for (others are linear
% combinations of other l-inf norms)
[num, den] = ss2tf(A, B, [0 5/64 0], 52/64);
Gur4 = tf(num, den, -1)
[num, den] = ss2tf(A, B, [-5/64 0 0], 78/64);
Gur5 = tf(num, den, -1)
[gur4, k] = impulse(Gur4);
plot(k, gur4)
[gur5, k] = impulse(Gur5);
plot(k, gur5)
legend('g_{u->x1}[k]', 'g_{u->x2}[k]', 'g_{u->x3}[k]', 'g_{u->r4}[k]', 'g_{u->r5}[k]')
gur4_l1 = max(abs(gur4))
gur5_l1 = max(abs(gur5))
r4_linf = gur4_l1 * u_linf
r5_linf = gur5_l1 * u_linf

% number of integer bits
r1_nInt = ceil(log2(ceil( u_linf * 52/64 ))) + 1
r2_nInt = ceil(log2(ceil( u_linf * 78/64 ))) + 1
r3_nInt = ceil(log2(ceil( u_linf ))) + 1
r4_nInt = ceil(log2(ceil( r4_linf ))) + 1
r5_nInt = ceil(log2(ceil( r5_linf ))) + 1
disp("r6 redundant with x_3")
r7_nInt = ceil(log2(ceil( x1_linf * 37/64 ))) + 1
r8_nInt = ceil(log2(ceil( x1_linf * 5/64 ))) + 1
r9_nInt = ceil(log2(ceil( x2_linf * 37/64 ))) + 1
r10_nInt = ceil(log2(ceil( x2_linf * 5/64 ))) + 1
r11_nInt = ceil(log2(ceil( x3_linf * 16/64 ))) + 1
r12_linf = x1_linf * 96/64;
r12_nInt = ceil(log2(ceil( r12_linf ))) + 1
r13_linf = x2_linf * 32/64;
r13_nInt = ceil(log2(ceil( r13_linf ))) + 1
r14_nInt = ceil(log2(ceil( x3_linf ))) + 1
r15_nInt = ceil(log2(ceil( r12_linf + r13_linf ))) + 1