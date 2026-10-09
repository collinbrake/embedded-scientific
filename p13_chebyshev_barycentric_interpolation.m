x = linspace(0, 1, 401);

% the natural mollifier function to interpolate on range [0, 1]
% (defined to be exactly 0 at x = 1, where the formula is singular)
f = exp(1./(x.^2 - 1));
f(x >= 1) = 0;

% numerator of the chebyshev barycentric polynomial interpolation
n = -1*exp(-8/5)./(x-3/4) + 1*exp(-8/7)./(x-1/4) - 1/2*exp(-1)./x;

% denominator
d = 1/2./(x-1) - 1./(x-3/4) + 1./(x-1/4) - 1/2./x;

% interpolating polynomial (n./d is 0/0-removable at node x = 0)
p = n./d;
p(x == 1) = 0;
p(x == 3/4) = exp(-8/5);
p(x == 1/4) = exp(-8/7);
p(x == 0) = exp(-1);

% 13 fitting points -> 12-order polynomial p12(x)
nNodes = 13;
nOrder = nNodes - 1;
j = 0:nOrder;
% Chebyshev nodes on [-1, 1], mapped onto the domain [0, 1]
xj = (1 - cos(j*pi/nOrder))/2;
fj = exp(-1./(1 - xj.^2));
fj(xj >= 1) = 0;

% barycentric weights for chebyshev
w = (-1).^j;
w(1) = 1/2*w(1);
w(nNodes) = 1/2*w(nNodes);

% barycentric form of the interpolating polynomial
p12 = zeros(size(x));
for k = 1:numel(x)
    xk = x(k);
    match = find(xk == xj, 1);
    if ~isempty(match)
        p12(k) = fj(match);
    else
        terms = w ./ (xk - xj);
        p12(k) = sum(terms .* fj) / sum(terms);
    end
end

% Error/RMS calcs
error_p = p-f;
error_p12 = p12-f;
rms_p = sqrt(sum(error_p.^2)/length(x))
rms_p12 = sqrt(sum((error_p12).^2)/length(x))

% plot
figure(1)
subplot(2, 1, 1)
plot(x, f, x, p, x, p12, '--')
xlabel("x")
legend("f(x)", "p_{3}(x)", "p_{12}(x)")

subplot(2, 1, 2)
plot(x, error_p, x, error_p12)
xlabel("x")
legend("Error p_{3}(x)", "Error p_{12}(x)")