
w0 = 5*pi*10^3;
Hcontinuous = tf([w0^2], [1 sqrt(2)*w0 w0^2]);
Hexact = c2d(Hcontinuous, 167*10^-6)


omega = linspace(0, pi, 1001);
z = exp(1i.*omega);
% Yexact = (n(1).*z.^2 + n(2).*z + n(3)) ./ (d(1).*z.^2 + d(2).*z + d(3));
Yexact = response(omega, Hexact);
Heuler = 6.88133 ./ (z.^2 + 1.70981.*z + 4.17152);
Htustin = 0.37601.*(z.^2 + 2.*z + 1) ./ (z.^2 + 0.314883.*z + 0.18915);

figure(1)
subplot(2, 1, 1)
plot(omega, abs(Heuler), omega, abs(Htustin), omega, abs(Yexact))
subplot(2, 1, 2)
plot(omega, angle(Heuler), omega, angle(Htustin), omega, angle(Yexact))

function Y = response(omega, H)
    n = H.Numerator{1}; % transpose to multiply by z matrix
    d = H.Denominator{1};
    z = exp(1i.*omega);
    zMat = [z.^2; z.^1; z.^0];
    Y = (n*zMat) ./ (d*zMat);
end
