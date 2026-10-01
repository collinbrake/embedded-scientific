x = linspace(0, 1, 101);

% the natural mollifier function to interpolate on range [0, 1]
f = exp(1./(x.^2 - 1));

% numerator of the lagrange barycentric polynomial interpolation
n = exp(-1)./x - 3*exp(-9/8)./(x - 1/3) + 3*exp(-9/5)./(x - 2/3);

% denominator
d = 1./x - 3./(x - 1/3) + 3./(x - 2/3) - 1./(x - 1);

% interpolating polynomial
p = n./d;

% plot
figure(1)
subplot(2, 1, 1)
plot(x, f, x, p)
xlabel("x")
legend("f(x)", "p(x)")

subplot(2, 1, 2)
plot(x, p-f)
xlabel("x")
legend("Error p(x) - f(x)")