% define the highest harmonic the approximations should use
harmonic = 3
harmonic2 = 10

% define one period of the function to approximate
t = linspace(-1, 1, 20001); t(end) = [];
t = t';
f = t; % sawtooth
T = 1;

X = ones(size(t));
X2 = X;

for l = 1:harmonic
    X = [X, sin(l*pi/T*t), cos(l*pi/T*t)];
end

for l = 1:harmonic2
    X2 = [X2, sin(l*pi/T*t), cos(l*pi/T*t)];
end

A = inv(X'*X)*X'*f
error = f - X*A;

A2 = inv(X2'*X2)*X2'*f
error2 = f - X2*A2;

subplot(2, 1, 1)
plot(t, X*A, t, X2*A2, t, f)
legend("LS approx, harmonic " + harmonic, "LS approx, harmonic " + harmonic2, 'Original function')
title("Fourier series approximations (least squares coefficient fit)")
grid on

subplot(2, 1, 2)
plot(t, error, t, error2)
legend("Error for harmonic " + harmonic, "Error for harmonic " + harmonic2)
title("Approximation errors")
grid on