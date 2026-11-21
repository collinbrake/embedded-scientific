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

legend('g_{u->x1}[k]', 'g_{u->x2}[k]', 'g_{u->x3}[k]')

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
x1_NumIntBits = ceil(log2(ceil(x1_linf))) + 1
x2_NumIntBits = ceil(log2(ceil(x2_linf))) + 1
x3_NumIntBits = ceil(log2(ceil(x3_linf))) + 1

disp("============================= (H) ===============================")

