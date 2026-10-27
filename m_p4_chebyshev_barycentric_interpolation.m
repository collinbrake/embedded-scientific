% range [-2, 5.5]
a = -3;
b = 1;
x = a:0.0001:b;
size(x)

% function to interpolate on range [-3, 1]
f = cos(x.^2 - 1./(x-1.17).^2);

% fitting points
nNodes = 97;
nOrder = nNodes - 1;
j = 0:nOrder;

% Chebyshev nodes on [-1, 1], mapped onto the domain [a, b]
xj = (a + (b-a)/2*(cos(j*pi/nOrder) + 1));
fj = cos(xj.^2 - 1./(xj-1.17).^2);

% barycentric weights for chebyshev
w = (-1).^j;
w(1) = 1/2*w(1);
w(nNodes) = 1/2*w(nNodes);

% barycentric form of the interpolating polynomial
p = zeros(size(x));
for k = 1:numel(x)
    xk = x(k);
    match = find(xk == xj, 1);
    if ~isempty(match)
        p(k) = fj(match);
    else
        terms = w ./ (xk - xj);
        p(k) = sum(terms .* fj) / sum(terms);
    end
end

% Error/RMS calcs
error_p = p-f;
rms_p = sqrt(sum((error_p).^2)/length(x))

% plot
figure(1)
subplot(2, 1, 1)
plot(x, f, x, p)
xlabel("x")
legend("f(x)", "p(x)")
title("Problem 04 (B)")

subplot(2, 1, 2)
plot(x, error_p)
xlabel("x")
legend("Error p(x)")
title("Problem 04 (C)")