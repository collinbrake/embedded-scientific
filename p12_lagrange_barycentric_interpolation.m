% Problem 12: (Lagrange) Barycentric Polynomial Approximation
%
% Runge's example f(x) = 1 / (1 + 25x^2) on [-1, 1], approximated with
% evenly spaced fitting points at increasing polynomial order.

a = -1;
b = 1;
x = linspace(a, b, 1001);
f = @(x) 1 ./ (1 + 25*x.^2);

% (A) sixth-order approximation, seven evenly spaced fitting points
figure(1)
subplot(2, 3, 1)
barycentric_interp(f, a, b, 7, x, "fit");
title("(A) p_6(x) vs f(x)")

subplot(2, 3, 4)
barycentric_interp(f, a, b, 7, x, "error");
title("(A) Error f(x) - p_6(x)")

% (B) twelfth-order approximation, thirteen evenly spaced fitting points
subplot(2, 3, 2)
barycentric_interp(f, a, b, 13, x, "fit");
title("(B) p_{12}(x) vs f(x)")

subplot(2, 3, 5)
barycentric_interp(f, a, b, 13, x, "error");
title("(B) Error f(x) - p_{12}(x)")

% (C) 100th-order approximation, one-hundred-and-one evenly spaced fitting points
subplot(2, 3, 3)
barycentric_interp(f, a, b, 101, x, "fit");
title("(C) p_{100}(x) vs f(x)")

subplot(2, 3, 6)
barycentric_interp(f, a, b, 101, x, "error");
title("(C) Error f(x) - p_{100}(x)")

% (D) compare the three approximations by RMS error
fprintf("%-18s %s\n", "Fitting points", "RMS error");
for nNodes = [7, 13, 101]
    p = barycentric_interp(f, a, b, nNodes, x);
    rmsError = sqrt(mean((f(x) - p).^2));
    fprintf("%-18d %.6g\n", nNodes, rmsError);
end

function p = barycentric_interp(f, a, b, nNodes, x, plotMode)
%BARYCENTRIC_INTERP Evenly spaced barycentric Lagrange interpolation.
%   p = BARYCENTRIC_INTERP(f, a, b, nNodes, x, plotMode) fits nNodes
%   evenly spaced points on [a, b] to f, evaluates the barycentric
%   interpolant at x, and returns p. plotMode is "fit" (plot f and p),
%   "error" (plot f - p), or "none" (default, no plot) on the current axes.

if nargin < 6
    plotMode = "none";
end

xj = linspace(a, b, nNodes);
fj = f(xj);

% barycentric weights for evenly spaced nodes (second form): lambda_j = (-1)^j * nchoosek(n, j)
n = nNodes - 1;
j = 0:n;
w = (-1).^j .* arrayfun(@(k) nchoosek(n, k), j);

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
        error("barycentric_interp:badMode", "Unknown plotMode '%s'", plotMode)
end

end