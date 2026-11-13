% ============================== (A) ================================== %

clear all

dt = 0.2;

den_common = dt^4 + 5*dt^3 + 10*dt^2 + 10*dt + 4;

b4 = (dt^4) / den_common;
a3 = (4*dt^4 + 10*dt^3 - 20*dt - 16) / den_common;
a2 = (6*dt^4 - 20*dt^2 + 24)          / den_common;
a1 = (4*dt^4 - 10*dt^3 + 20*dt - 16) / den_common;
a0 = (dt^4 - 5*dt^3 + 10*dt^2 - 10*dt + 4) / den_common;

num_d = b4 * [1, 4, 6, 4, 1];
den_d = [1, a3, a2, a1, a0];

num_c = [0, 0,  0,  0, 4];
den_c = [1, 5, 10, 10, 4];

OMEGA = linspace(0, pi, 10001); % normalized digital frequency [0, pi]
omega = OMEGA/dt;                % continuous frequency [0, pi/T]
sFreq = 1i.*omega;
zFreq = exp(1i.*OMEGA);
sMat = [sFreq.^4; sFreq.^3; sFreq.^2; sFreq.^1; sFreq.^0];
zMat = [zFreq.^4; zFreq.^3; zFreq.^2; zFreq.^1; zFreq.^0];

Yexact = (num_c*sMat) ./ (den_c*sMat);

Yd = (num_d*zMat) ./ (den_d*zMat);

Yerror = abs(Yexact - Yd);

figure(1)
subplot(3, 1, 1)
plot(OMEGA, 20*log10(abs(Yexact)), ...
    'Color', [0.15, 0.25, 0.55], 'LineWidth', 2.2, 'LineStyle', '-');
hold on;
plot(OMEGA, 20*log10(abs(Yd)), ...
    'Color', [0.95, 0.60, 0.20], 'LineWidth', 1.8, 'LineStyle', '--');
hold off;
xlabel('rad'); ylabel('Magnitude (dB)'); grid on
legend('20*log10|H_{d,exact}(e^{i\Omega})| = 20*log10|H_c(i\Omega/dt)|', '20*log10|H_d(e^i^\Omega)|', 'Location', 'southwest')
title('Magnitude Response (dB)')

subplot(3, 1, 2)
plot(OMEGA, unwrap(angle(Yexact)), ...
    'Color', [0.15, 0.25, 0.55], 'LineWidth', 2.2, 'LineStyle', '-');
hold on;
plot(OMEGA, unwrap(angle(Yd)), ...
    'Color', [0.95, 0.60, 0.20], 'LineWidth', 1.8, 'LineStyle', '--');
hold off;
xlabel('rad'); ylabel('Phase (rad)'); grid on
legend('\angle H_{d,exact}(e^{i\Omega}) = \angle H_c(i\Omega/dt)', '\angle H_d(e^{i\Omega})', 'Location', 'southwest');
title('Phase Response (rad)')

subplot(3, 1, 3)
plot(OMEGA, Yerror, 'Color', [0.95, 0.60, 0.20], 'LineWidth', 1.8)
xlabel('rad'); ylabel('Magnitude Error (Unitless)'); grid on
title('Magnitude Response Error (Unitless)')
sgtitle('Homework 19 - Frequency Response (Discrete Time)')

% ============================== (B) ================================== %

% Direct Form II implementation
A = [0  1   0   0; ...
    0   0   1   0; ...
    0   0   0   1; ...
    -a0 -a1 -a2 -a3]
B = [0; 0; 0; 1]

% find the strictly proper tf part numerator
beta3 = b4 * (4 - a3);
beta2 = b4 * (6 - a2);
beta1 = b4 * (4 - a1);
beta0 = b4 * (1 - a0);
C = [beta0, beta1, beta2, beta3]

D = b4

% ============================== (C) ================================== %

% Eigenvalues and distance to unit circle
lambda = eig(A);
rho = abs(lambda);
theta = angle(lambda);
d = abs(1 - rho);

polar_form = compose('%.6f * exp(i * %.5f)', rho, theta);

T = table(polar_form, d, 'VariableNames', {'Polar_Form', 'Distance_d'});
disp(T)

% ============================== (E) ================================== %

% Diagonal form

[P, Adiag] = eig(A);

Adiag
Bdiag = inv(P)*B
Cdiag = C*P
Ddiag = D

% ============================== (F) ================================== %

% Parallel as possible

Ppar = [real(P(:,1)), imag(P(:,1)), P(:,3), P(:,4)];

Apar = inv(Ppar) * A * Ppar
Bpar = inv(Ppar) * B
Cpar = C*Ppar
Dpar = D

% ============================== (I) ================================== %

% Parallel as possible - implemented in fixed point [16, -10]
% representation

K = 1000;
Pbar = Ppar*K;

A_inf = inv(Pbar) * A * Pbar;
B_inf = inv(Pbar) * B;
C_inf = C*Pbar;
D_inf = D;

% Scale factor for [16, -10] fixed-point representation
scale = 2^10; % 1024

% Direct round-to-nearest fixed-point conversion
Abar = round(A_inf * scale) / scale;
Bbar = round(B_inf * scale) / scale;
Cbar = round(C_inf * scale) / scale;
Dbar = 1 / scale; % truncates to zero

% Display quantized matrices
disp('Abar:'); disp(Abar);
disp('Bbar:'); disp(Bbar);
disp('Cbar:'); disp(Cbar);
disp('Dbar:'); disp(Dbar);

% ============================== (J) ================================== %
% Eigenvalues of approximated Abar in polar form
lambda_bar = eig(Abar);
rho_bar = abs(lambda_bar);
theta_bar = angle(lambda_bar);

polar_form_bar = compose('%.6f * exp(i * %.5f)', rho_bar, theta_bar);
T_J = table(polar_form_bar, 'VariableNames', {'Approximated_Polar_Form'});
disp(T_J)

% ============================== (K) ================================== %
% Distance d_bar between each approximated eigenvalue and the unit circle
d_bar = abs(1 - rho_bar);

T_K = table(polar_form_bar, d_bar, 'VariableNames', {'Polar_Form_Bar', 'Distance_d_bar'});
disp(T_K)

% ============================== (L) ================================== %
% Relative movement of each eigenvalue: X = |d_l - d_bar_l| / d_l
% (Sort both sets by real/imag parts to ensure matching pairs)
[~, idx_orig] = sort(real(lambda) + imag(lambda));
[~, idx_bar]  = sort(real(lambda_bar) + imag(lambda_bar));

d_sorted     = d(idx_orig);
d_bar_sorted = d_bar(idx_bar);

X = abs(d_sorted - d_bar_sorted) ./ d_sorted;
X_pct = X * 100;

for l = 1:length(X)
    fprintf('eigenvalue lambda_%d moved a total of %.4f%% of the distance between lambda_%d and the unit circle.\n', ...
        l, X_pct(l), l);
end