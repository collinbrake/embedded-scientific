% Selected clean segment
x_fit = [x_j(1:6:39), x_j(39)]; % make sure to get the last point
f_fit = [f_j(1:6:39), f_j(39)]; % make sure to get the last point
M = length(x_fit);

w = zeros(1, M);
for j = 1:M
    diffs = x_fit(j) - x_fit([1:j-1, j+1:M]);
    w(j) = 1 / prod(diffs);
end

% display fitting points, function values, and Barycentric weights:
fprintf('Fitting points:\n');
x_fit
fprintf('Function values:\n');
f_fit
fprintf('Barycentric weights w:\n');
fprintf('%1.4e\n', w);
% x_fit
% f_fit
