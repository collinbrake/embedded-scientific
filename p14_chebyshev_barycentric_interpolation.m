% Problem 14: (Chebyshev) Barycentric Polynomial Approximation
%
% Runge's example f(x) = 1 / (1 + 25x^2) on [-1, 1], approximated with
% Chebyshev-Lobatto fitting points at increasing polynomial order.

a = -1;
b = 1;
x = linspace(a, b, 1001);
f = @(x) 1 ./ (1 + 25*x.^2);

% (A) sixth-order approximation, seven Chebyshev fitting points
figure(1)
subplot(2, 3, 1)
chebyshev_interp(f, a, b, 7, x, "fit");
title("(A) p_6(x) vs f(x)")

subplot(2, 3, 4)
chebyshev_interp(f, a, b, 7, x, "error");
title("(A) Error f(x) - p_6(x)")

% (B) twelfth-order approximation, thirteen Chebyshev fitting points
subplot(2, 3, 2)
chebyshev_interp(f, a, b, 13, x, "fit");
title("(B) p_{12}(x) vs f(x)")

subplot(2, 3, 5)
chebyshev_interp(f, a, b, 13, x, "error");
title("(B) Error f(x) - p_{12}(x)")

% (C) 100th-order approximation, one-hundred-and-one Chebyshev fitting points
subplot(2, 3, 3)
chebyshev_interp(f, a, b, 101, x, "fit");
title("(C) p_{100}(x) vs f(x)")

subplot(2, 3, 6)
chebyshev_interp(f, a, b, 101, x, "error");
title("(C) Error f(x) - p_{100}(x)")

% (D) compare the three approximations by RMS error
fprintf("%-18s %s\n", "Fitting points", "RMS error");
for nNodes = [7, 13, 101]
    p = chebyshev_interp(f, a, b, nNodes, x);
    rmsError = sqrt(mean((f(x) - p).^2));
    fprintf("%-18d %.6g\n", nNodes, rmsError);
end

function p = chebyshev_interp(f, a, b, nNodes, x, plotMode)
if nargin < 6
    plotMode = "none";
end

n = nNodes - 1;
j = 0:n;
% Chebyshev nodes on [-1, 1], mapped onto the domain [a, b]
xj = a + (b - a)*(1 - cos(j*pi/n))/2;
fj = f(xj);

% barycentric weights for chebyshev
w = (-1).^j;
w(1) = 1/2*w(1);
w(nNodes) = 1/2*w(nNodes);

% barycentric form, falling back to fj at the nodes themselves (0/0)
p = zeros(size(x));
for i = 1:numel(x)
    diffs = x(i) - xj;
    exact = find(diffs == 0, 1);
    if isempty(exact)
        terms = w ./ diffs;
        p(i) = sum(terms .* fj) / sum(terms);
    else
        p(i) = fj(exact);
    end
end

switch plotMode
    case "fit"
        plot(x, f(x), x, p, '--')
        xlabel("x")
        legend("f(x)", sprintf("p_{%d}(x)", nNodes - 1))
    case "error"
        plot(x, f(x) - p)
        xlabel("x")
        legend(sprintf("f(x) - p_{%d}(x)", nNodes - 1))
    case "none"
        % return p without plotting
    otherwise
        error("chebyshev_interp:badMode", "Unknown plotMode '%s'", plotMode)
end

end