% Parameters
T = 500*10^-6

% Euler
He = tf([0, 0, 2*T^2],...
            [1, -2+2*T, 1-2*T+2*T^2], T)

% Tustin
Htustin = tf([w0^2*T^2, 2*w0^2*T^2, w0^2*T^2],...
                [4+w0^2*T^2+2*sqrt(2)*w0*T,...
                    -8+2*w0^2*T^2, 4-2*sqrt(2)*w0*T+w0^2*T^2], T)

% Exact discretization with ZOH
Hcontinuous = tf([2], [1, 2, 2]);
Hexact = c2d(Hcontinuous, T)
Ht = c2d(Hcontinuous, T, 'tustin')

omega = linspace(0, pi, 1001);
z = exp(1i.*omega);

Yexact = response(omega, Hexact);
Yeuler = response(omega, He);
Ytustin = response(omega, Ht);

figure(1)
subplot(2, 1, 1)
plot(omega, abs(Yeuler), omega, abs(Ytustin), omega, abs(Yexact))
subplot(2, 1, 2)
plot(omega, unwrap(angle(Yeuler)), omega, unwrap(angle(Ytustin)), omega,...
        unwrap(angle(Yexact)))

function Y = response(omega, H)
    n = H.Numerator{1};
    d = H.Denominator{1};
    z = exp(1i.*omega);
    zMat = [z.^2; z.^1; z.^0];
    Y = (n*zMat) ./ (d*zMat);
end
