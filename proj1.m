% ============================== (A) ================================== %

fc = 698.46;
wc = 2*pi*fc % corner frequency, rad/s

s = tf('s');

p1 = wc*exp(1i*5*pi/8);
p2 = conj(p1);
p3 = wc*exp(1i*7*pi/8);
p4 = conj(p3);

% verify conjugate poles:
distp2 = abs(p2 - wc*exp(1i*11*pi/8))
distp4 = abs(p4 - wc*exp(1i*9*pi/8))

Hc1 = 1 / ((s - p1)*...
            (s - p2)) % part of Hc containing one pair of cmplx conj poles

Hc2 = 1 / ((s - p3)*...
            (s - p4)) % part of Hc containing one pair of cmplx conj poles

Hc_unscaled = Hc1*Hc2;

D = Hc_unscaled.Denominator{1};
invDC = D(end); % DC gain, s = 0 eliminates all but last term

KHc = invDC*4 % select K_Hc so that overall DC gain is 4

Hc = KHc*Hc_unscaled;

% ============================== (B) ================================== %

T = 1 / (90*2*10^3) % sampling frequency has to be at least double 2kHz signal

Hd = c2d(Hc, T, 'Tustin')

Hd_exact = c2d(Hc, T);

% ========================= Responses ================================= %

OMEGA = linspace(0, pi, 1000001); % normalized digital frequency [0, pi]
omega = OMEGA/T;                % continuous frequency [0, pi/T]
sFreq = 1i.*omega;
zFreq = exp(1i.*OMEGA);
sMat = [sFreq.^4; sFreq.^3; sFreq.^2; sFreq.^1; sFreq.^0];
zMat = [zFreq.^4; zFreq.^3; zFreq.^2; zFreq.^1; zFreq.^0];

nc = Hc.Numerator{1};
dc = Hc.Denominator{1};
Yc = (nc*sMat) ./ (dc*sMat);

nd = Hd.Numerator{1};
dd = Hd.Denominator{1};
Yd = (nd*zMat) ./ (dd*zMat);

nd = Hd_exact.Numerator{1};
dd = Hd_exact.Denominator{1};
Yd_exact = (nd*zMat) ./ (dd*zMat);

Yerror_norm = abs(Yd_exact - Yd) / max(abs(Yd_exact));

% (B) - not for report, just to iterate on sample time T
figure(1)
subplot(3, 1, 1)
semilogx(omega, 20*log10(abs(Yd)), omega, 20*log10(abs(Yd_exact)))
xlabel('rad/s'); ylabel('Magnitude (dB)'); grid on
legend('H_d(s)', 'H_d_e_x_a_c_t(z)', 'Location', 'southwest')
title('Magnitude Response (dB)')
subplot(3, 1, 2)
semilogx(omega, unwrap(angle(Yd)), omega, unwrap(angle(Yd_exact)))
xlabel('rad/s'); ylabel('Phase (rad)'); grid on
legend('H_d(s)', 'H_d_e_x_a_c_t(z)', 'Location', 'southwest')
title('Phase Response (rad)')
subplot(3, 1, 3)
semilogx(omega, abs(Yerror_norm*100))
xlabel('rad/s'); ylabel('Magnitude Error (%)'); grid on
title('Magnitude Response Error (%)')
sgtitle('Project (1) - Frequency Response of 4th-order Butterworth Filter (Discrete Time)')

% (E) - unitless, with phase
figure(2)
subplot(2, 1, 1)
semilogx(omega, abs(Yc), ...
    'Color', [0.15, 0.25, 0.55], 'LineWidth', 2.2, 'LineStyle', '-');
hold on;
semilogx(omega, abs(Yd), ...
    'Color', [0.95, 0.60, 0.20], 'LineWidth', 1.8, 'LineStyle', '--');
hold off;
xlabel('rad/s'); ylabel('Magnitude (Unitless)'); grid on
legend('|H_c(i\omega)|', '|H_d(e^i^\Omega)|', 'Location', 'southwest')

subplot(2, 1, 2)
semilogx(omega, unwrap(angle(Yc)), ...
    'Color', [0.15, 0.25, 0.55], 'LineWidth', 2.2, 'LineStyle', '-');
hold on;
semilogx(omega(1:end-100), unwrap(angle(Yd(1:end-100))), ...
    'Color', [0.95, 0.60, 0.20], 'LineWidth', 1.8, 'LineStyle', '--');
hold off;
xlabel('rad/s'); ylabel('Phase (rad)'); grid on
legend([char(8736) 'H_c(i\omega)'], [char(8736) 'H_d(e^i^\Omega)'], 'Location', 'southwest')
sgtitle('Project (1) - Frequency Response of 4th-order Butterworth Filter (Continuous Time versus Discrete Time)')

% (E) - dB
figure(3)
semilogx(omega, 20*log10(abs(Yc)), ...
    'Color', [0.15, 0.25, 0.55], 'LineWidth', 2.2, 'LineStyle', '-');
hold on;
semilogx(omega, 20*log10(abs(Yd)), ...
    'Color', [0.95, 0.60, 0.20], 'LineWidth', 1.8, 'LineStyle', '--');
hold off;
xlabel('rad/s'); ylabel('Magnitude (dB)'); grid on
legend('20*log10|H_c(i\omega)|', '20*log10|H_d(e^i^\Omega)|', 'Location', 'southwest')
sgtitle('Project (1) - Frequency Response of 4th-order Butterworth Filter (Continuous Time versus Discrete Time)')