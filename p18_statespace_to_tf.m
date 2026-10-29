A = [0.4103125 0.3578125 -0.3996875; 0 -1 1; 0.5896875 -0.3578125 0.3996875]
B = [0.5; 0; -0.5]
C = [0.125 0.75 1]
D = 1.35

z = tf('z');
I = eye(size(A));

H = C * ((z*I - A) \ B) + D;
H = minreal(H, 1e-5)

Hexact = tf(ss(A, B, C, D, -1))

poles = roots(H.Denominator{1})
poles_exact = roots(Hexact.Denominator{1})
eigenVals = eig(A)

% Responses
omega = linspace(0, pi, 10001);
zr = exp(1i.*omega);
zMat = [zr.^3; zr.^2; zr.^1; zr.^0];

n = H.Numerator{1};
d = H.Denominator{1};
Y = (n*zMat) ./ (d*zMat);

n = Hexact.Numerator{1};
d = Hexact.Denominator{1};
Yexact = (n*zMat) ./ (d*zMat);

% Magnitude and Phase Plots
figure(1)
subplot(2, 1, 1)
plot(omega, abs(Y), omega, abs(Yexact), '--')
xlabel('rad'); ylabel('Magnitude (unitless)'); grid on
legend('H(z) - C * ((z*I - A)^-^1 * B) + D', 'H(z) - exact', 'Location', 'southwest')
subplot(2, 1, 2)
plot(omega, unwrap(angle(Y)), omega, unwrap(angle(Yexact)), '--')
xlabel('rad'); ylabel('Phase (rad)'); grid on
legend('H(z) - C * ((z*I - A)^-^1 * B) + D', 'H(z) - exact', 'Location', 'southwest')
sgtitle('Homework 18 (C) - State Space Rep to Transfer Function')
