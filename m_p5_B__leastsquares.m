% define the highest harmonic the approximations should use
harmonic = 3
harmonic2 = 10
harmonic3 = 50 % the assignment says 101 basis functions and 100th harmonic -- going with 101 basis functions

% define one period of the function to approximate
t = linspace(-1, 1, 20001); t(end) = [];
t = t';
f1 = (t+1).^(1/2); % sharkfin
f2 = (t+1).^(2);
f = [f1(1:10000); f2(10001:end)];
T = 1;

X = ones(size(t));
X2 = X;
X3 = X;

for l = 1:harmonic
    X = [X, sin(l*pi/T*t), cos(l*pi/T*t)];
end

for l = 1:harmonic2
    X2 = [X2, sin(l*pi/T*t), cos(l*pi/T*t)];
end

for l = 1:harmonic3
    X3 = [X3, sin(l*pi/T*t), cos(l*pi/T*t)];
end

A = inv(X'*X)*X'*f
error = f - X*A;

A2 = inv(X2'*X2)*X2'*f
error2 = f - X2*A2;

A3 = inv(X3'*X3)*X3'*f
error3 = f - X3*A3;

subplot(2, 1, 1)
plot(t, X*A, t, X2*A2, t, X3*A3, t, f)
legend("LS approx, harmonic " + harmonic, "LS approx, harmonic " + harmonic2, "LS approx, harmonic " + harmonic3, 'Original function')
title("Fourier series approximations (least squares coefficient fit)")
grid on

subplot(2, 1, 2)
plot(t, error, t, error2, t, error3)
legend("Error for harmonic " + harmonic, "Error for harmonic " + harmonic2, "Error for harmonic " + harmonic3)
title("Approximation errors")
grid on