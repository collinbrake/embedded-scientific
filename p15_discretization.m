omega = linspace(0, 2*pi, 1001);
z = exp(i.*omega);
Heuler = 6.88133 ./ (z.^2 + 1.70981.*z + 4.17152);
figure(1)
subplot(2, 1, 1)
plot(omega, abs(Heuler))
plot(omega, angle(Heuler))
