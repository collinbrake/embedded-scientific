% ====================================================================== %
% Problem 20: A Fourth-Order Filter
% Embedded Scientific Computing 2026
% Collin Brake
% ====================================================================== %

clear all
close all
clc

% Transfer function definition
% H(z) = (z^3 - 1.467z^2 - 0.987z + 0.2963) / (z^4 - 6.495z^3 + 6.1624z^2 - 1.6657z - 0.098765)

num = [1.0, -1.467, -0.987, 0.2963];
den = [1.0, -6.495, 6.1624, -1.6657, -0.098765];

% ============================== (A) ================================== %
% Direct-form II state-space representation
% (Numerator in C matrix, negative denominator in bottom row of A matrix)
% Note D = 0 since numerator degree (3) < denominator degree (4)

a3 = -6.495;
a2 = 6.1624;
a1 = -1.6657;
a0 = -0.098765;

b3 = 1.0;
b2 = -1.467;
b1 = -0.987;
b0 = 0.2963;

A = [0, 1, 0, 0; ...
     0, 0, 1, 0; ...
     0, 0, 0, 1; ...
    -a0, -a1, -a2, -a3];

B = [0; 0; 0; 1];

C = [b0, b1, b2, b3];

D = 0;

disp('(A) Direct-Form II State-Space:');
disp('A ='); disp(A);
disp('B ='); disp(B);
disp('C ='); disp(C);
disp('D ='); disp(D);

% Verify transfer function matches original
H_df2 = tf(ss(A, B, C, D, -1));
disp('Verified Transfer Function:');
disp(H_df2);

% ============================== (B) ================================== %
% Eigenvalues (poles) and eigenvectors of A
[V, D_eig] = eig(A);
lambda = diag(D_eig);

disp('(B) Eigenvalues and Eigenvectors:');
disp('Eigenvalues lambda ='); disp(lambda);
disp('Eigenvector matrix V ='); disp(V);

% Display in polar form
rho = abs(lambda);
theta = angle(lambda);
polar_form = compose('%.6f * exp(i * %.5f)', rho, theta);
T_B = table(lambda, rho, theta, polar_form, 'VariableNames', {'Eigenvalue', 'Magnitude_rho', 'Angle_theta_rad', 'Polar_Form'});
disp(T_B);

% ============================== (C) ================================== %
% Diagonal form
% A_diag = inv(V)*A*V, B_diag = inv(V)*B, C_diag = C*V
P_diag = V;
Adiag = inv(P_diag) * A * P_diag;
Bdiag = inv(P_diag) * B;
Cdiag = C * P_diag;

disp('(C) Diagonal Form (Complex Multipliers):');
disp('Adiag ='); disp(Adiag);
disp('Bdiag ='); disp(Bdiag);
disp('Cdiag ='); disp(Cdiag);
disp('Why not implementable: Complex conjugate eigenvalues yield complex multipliers in Adiag, Bdiag, and Cdiag.');

% ============================== (D) ================================== %
% Coupled-form state-space representation ("as diagonal as possible")
% Couple complex-conjugate eigenvalues into 2x2 block, keep real eigenvalues separate.
% All entries must be real-valued.

% Identify complex and real eigenvalues
idx_c = find(abs(imag(lambda)) > 1e-8);
idx_r = find(abs(imag(lambda)) <= 1e-8);

% Complex pair eigenvector real and imaginary parts
v_c = V(:, idx_c(1));
v_R = real(v_c);
v_I = imag(v_c);

% Real eigenvectors
v_r1 = real(V(:, idx_r(1)));
v_r2 = real(V(:, idx_r(2)));

% Similarity transformation matrix
P_coup = [v_R, v_I, v_r1, v_r2];

A_coup = inv(P_coup) * A * P_coup;
B_coup = inv(P_coup) * B;
C_coup = C * P_coup;

disp('(D) Coupled-Form State-Space (Real-Valued):');
disp('A_coup ='); disp(A_coup);
disp('B_coup ='); disp(B_coup);
disp('C_coup ='); disp(C_coup);

% ============================== (F) ================================== %
% Parallel-form state-space representation with two second-order blocks
% Block 1: Complex conjugate pair
% Block 2: Two real poles

% Partial fraction expansion using residue
[r, p, k_pfe] = residue(num, den);

% Find indices in residue output matching complex and real poles
idx_pfe_c = find(abs(imag(p)) > 1e-8);
idx_pfe_r = find(abs(imag(p)) <= 1e-8);

% Subsystem 1 (complex pair): poles p(idx_pfe_c)
r_c1 = r(idx_pfe_c(1));
p_c1 = p(idx_pfe_c(1));

% H1(z) = r_c1/(z - p_c1) + conj(r_c1)/(z - conj(p_c1))
%       = (2*real(r_c1)*z - 2*real(r_c1*conj(p_c1))) / (z^2 - 2*real(p_c1)*z + abs(p_c1)^2)
num1 = [2*real(r_c1), -2*real(r_c1 * conj(p_c1))];
den1 = [1, -2*real(p_c1), abs(p_c1)^2];

% Subsystem 2 (real poles): poles p(idx_pfe_r)
r_r1 = r(idx_pfe_r(1)); p_r1 = p(idx_pfe_r(1));
r_r2 = r(idx_pfe_r(2)); p_r2 = p(idx_pfe_r(2));

% H2(z) = r_r1/(z - p_r1) + r_r2/(z - p_r2)
num2 = [real(r_r1 + r_r2), -real(r_r1*p_r2 + r_r2*p_r1)];
den2 = conv([1, -real(p_r1)], [1, -real(p_r2)]);

% State space matrices for each second-order subsystem (Direct Form II canonical)
A1 = [0, 1; -den1(3), -den1(2)];
B1 = [0; 1];
C1 = [num1(2), num1(1)];

A2 = [0, 1; -den2(3), -den2(2)];
B2 = [0; 1];
C2 = [num2(2), num2(1)];

% Combined parallel form
A_par = [A1, zeros(2, 2); zeros(2, 2), A2];
B_par = [B1; B2];
C_par = [C1, C2];

disp('(F) Parallel-Form State-Space (Two 2nd-Order Blocks):');
disp('A_par ='); disp(A_par);
disp('B_par ='); disp(B_par);
disp('C_par ='); disp(C_par);

% ============================== (H) ================================== %
% Transfer functions for each of the second-order subsystems
H1 = tf(num1, den1, -1);
H2 = tf(num2, den2, -1);

disp('(H) Transfer functions of second-order subsystems:');
disp('H1(z) ='); disp(H1);
disp('H2(z) ='); disp(H2);
