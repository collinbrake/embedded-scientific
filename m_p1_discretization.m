% Parameters
T = 500*10^-6

% Euler - fixed point
Ke = 2^21
He = tf([0, 0, 0.125]/Ke, [1, -2, 1], T)

% Tustin - fixed point
Kt = 2^23
Ht = tf([0.125, 0.25, 0.125]/Kt, [1, -2, 1], T)

% Exact discretization with ZOH
Hcontinuous = tf([2], [1, 2, 2]);
Hexact = c2d(Hcontinuous, T)

omega = linspace(0, pi, 1001);
z = exp(1i.*omega);

Yexact = response(omega, Hexact);
Yeuler = response(omega, He);
Ytustin = response(omega, Ht);

figure(1)
subplot(2, 1, 1)
plot(omega, 20*log10(abs(Yeuler)), omega, 20*log10(abs(Ytustin)), omega, 20*log10(abs(Yexact)))
ylabel('Magnitude (dB)'); grid on
legend('Euler, fixed-point', 'Tustin, fixed-point', 'Exact (ZOH)', 'Location', 'southwest')
subplot(2, 1, 2)
plot(omega, unwrap(angle(Yeuler))-2*pi, omega, unwrap(angle(Ytustin)), omega,...
        unwrap(angle(Yexact)))
xlabel('\omega (rad)'); ylabel('Phase (rad)'); grid on
sgtitle('Problem 01 (G) - Frequency Responses')

function Y = response(omega, H)
    n = H.Numerator{1};
    d = H.Denominator{1};
    z = exp(1i.*omega);
    zMat = [z.^2; z.^1; z.^0];
    Y = (n*zMat) ./ (d*zMat);
end
