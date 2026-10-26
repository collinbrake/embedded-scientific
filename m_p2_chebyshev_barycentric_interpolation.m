% range [-2, 5.5]
a = -2;
b = 5.5;
x = linspace(a, b, 1001);

% function to interpolate on range [-2, 5.5]
f = 2*sin(x) - cos(exp(3*x.^2 ./ 4));

% 22 fitting points -> 21st-order polynomial p21(x)
nNodes = 22;
nOrder = nNodes - 1;
j = 0:nOrder;

% Chebyshev nodes on [-1, 1], mapped onto the domain [a, b]
xj = (a + (b-a)/2*(cos(j*pi/nOrder) + 1));
fj = 2*sin(xj) - cos(exp(3*xj.^2 ./ 4));

% barycentric weights for chebyshev
w = (-1).^j;
w(1) = 1/2*w(1);
w(nNodes) = 1/2*w(nNodes);

% barycentric form of the interpolating polynomial
p21 = zeros(size(x));
for k = 1:numel(x)
    xk = x(k);
    match = find(xk == xj, 1);
    if ~isempty(match)
        p21(k) = fj(match);
    else
        terms = w ./ (xk - xj);
        p21(k) = sum(terms .* fj) / sum(terms);
    end
end

% Error/RMS calcs
error_p21 = p21-f;
rms_p21 = sqrt(sum((error_p21).^2)/length(x))

% plot
figure(1)
subplot(2, 1, 1)
plot(x, f, x, p21, xj, fj, 'o')
xlabel("x")
legend("f(x)", "p_{21}(x)", "(xj, fj)")

subplot(2, 1, 2)
plot(x, error_p21)
xlabel("x")
legend("Error p_{21}(x)")

sgtitle("Problem 02 (A)")

% 102 fitting points -> 101st-order polynomial p101(x)
nNodes = 102;
nOrder = nNodes - 1;
j = 0:nOrder;

% Chebyshev nodes on [-1, 1], mapped onto the domain [a, b]
xj = (a + (b-a)/2*(cos(j*pi/nOrder) + 1));
fj = 2*sin(xj) - cos(exp(3*xj.^2 ./ 4));

% barycentric weights for chebyshev
w = (-1).^j;
w(1) = 1/2*w(1);
w(nNodes) = 1/2*w(nNodes);

% barycentric form of the interpolating polynomial
p101 = zeros(size(x));
for k = 1:numel(x)
    xk = x(k);
    match = find(xk == xj, 1);
    if ~isempty(match)
        p101(k) = fj(match);
    else
        terms = w ./ (xk - xj);
        p101(k) = sum(terms .* fj) / sum(terms);
    end
end

% Error/RMS calcs
error_p101 = p101-f;
rms_p101 = sqrt(sum((error_p101).^2)/length(x))

% plot
figure(2)
subplot(2, 1, 1)
plot(x, f, x, p101)
xlabel("x")
legend("f(x)", "p_{101}(x)")

subplot(2, 1, 2)
plot(x, error_p101)
xlabel("x")
legend("Error p_{101}(x)")

sgtitle("Problem 02 (B)")