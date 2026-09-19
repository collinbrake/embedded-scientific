x = [0:pi/99:pi]; 
f = cos(2*cos(2*cos(x)));

[m, n] = findpeaks(f, x) % print maxima
[m, n] = findpeaks(-f, x) % minima

% fitting points
xj = [0, 0.6664, 1.5549, 2.4752, pi]
fj = cos(2*cos(2*cos(xj)))

% cubic spline
cs = spline(xj, [0, fj, 0]); % 0's before and after y specify natural spline - second deriv 0

%% Part A - Cubic spline table

coef = cs.coefs;

intervals = strings(length(xj)-1,1);
cubics = strings(length(xj)-1,1);

for i = 1:length(xj)-1

    % Coefficients in shifted form
    a = coef(i,1);
    b = coef(i,2);
    c = coef(i,3);
    d = coef(i,4);

    % Expand:
    % a(x-xi)^3 + b(x-xi)^2 + c(x-xi) + d
    xi = xj(i);

    A = a;
    B = b - 3*a*xi;
    C = 3*a*xi^2 - 2*b*xi + c;
    D = -a*xi^3 + b*xi^2 - c*xi + d;

    intervals(i) = sprintf('[%.4f, %.4f]', xj(i), xj(i+1));

    cubics(i) = sprintf( ...
        '%.6fx^3 %+.6fx^2 %+.6fx %+.6f', ...
        A, B, C, D);
end

SplineTable = table(intervals, cubics, ...
    'VariableNames', {'Interval','Cubic'});

disp(SplineTable)

% plotting for part (B)
plot(xj, fj, 'ko', ...
     x, f, 'b-', ...
     x, ppval(cs,x), 'r--', ...
     'LineWidth', 1.5)

legend('Data points', 'Original function', 'Spline')
xlabel('x')
ylabel('f(x)')
grid on