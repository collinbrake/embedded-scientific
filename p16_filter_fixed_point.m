% Parameters
T = 0.1

% Filter - infinite precision
H = tf([1, 2, 1], [1, 1.609375, 0.658440655212500], T)

% Filter - fixed point [16, -6] representation
Hbar = tf([1, 2, 1], [1, 1.609375, 0.65625], T)

% Responses
omega = linspace(0, pi, 10001);
z = exp(1i.*omega);
zMat = [z.^2; z.^1; z.^0];

n = H.Numerator{1};
d = H.Denominator{1};
Y = (n*zMat) ./ (d*zMat);

n = Hbar.Numerator{1};
d = Hbar.Denominator{1};
Ybar = (n*zMat) ./ (d*zMat);

Yerror = Ybar - Y;

figure(1)
subplot(2, 1, 1)
plot(omega, abs(Y))
xlabel('rad'); ylabel('Magnitude (unitless)'); grid on
legend('H(z) - infinite precision', 'Location', 'southwest')
subplot(2, 1, 2)
plot(omega, unwrap(angle(Y)))
xlabel('rad'); ylabel('Phase (rad)'); grid on
sgtitle('Homework 16 (E) - Frequency Response of Infinite Precision H(z)')

figure(2)
subplot(2, 1, 1)
plot(omega, abs(Y), omega, abs(Ybar))
xlabel('rad'); ylabel('Magnitude (unitless)'); grid on
legend('H(z) - infinite precision', 'Hbar(z) - fixed-point', 'Location', 'southwest')
subplot(2, 1, 2)
plot(omega, unwrap(angle(Y)), omega, unwrap(angle(Ybar)))
xlabel('rad'); ylabel('Phase (rad)'); grid on
legend('H(z) - infinite precision', 'Hbar(z) - fixed-point', 'Location', 'southwest')
sgtitle('Homework 16 (G) - Frequency Response of Fixed-Point [-16, -6] Hbar(z)')

figure(3)
subplot(2, 1, 1)
semilogx(omega, 20*log10(abs(Y)), omega, 20*log10(abs(Ybar)))
xlabel('rad'); ylabel('Magnitude (dB)'); grid on
legend('H(z) - infinite precision', 'Hbar(z) - fixed-point', 'Location', 'southwest')
subplot(2, 1, 2)
semilogx(omega, unwrap(angle(Y)), omega, unwrap(angle(Ybar)))
xlabel('rad'); ylabel('Phase (rad)'); grid on
legend('H(z) - infinite precision', 'Hbar(z) - fixed-point', 'Location', 'southwest')
sgtitle('Homework 16 (G) - Frequency Response of Fixed-Point [-16, -6] Hbar(z)')

figure(4)
subplot(2, 1, 1)
plot(omega, abs(Yerror))
xlabel('rad'); ylabel('Magnitude (unitless)'); grid on
subplot(2, 1, 2)
plot(omega, unwrap(angle(Yerror)))
xlabel('rad'); ylabel('Phase (rad)'); grid on
sgtitle('Homework 16 (H) - Magnitude/Phase Response Error for Fixed Point Approximation of H(z)')