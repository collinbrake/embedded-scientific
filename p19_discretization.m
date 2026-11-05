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