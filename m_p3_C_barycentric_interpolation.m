% Fine evaluation grid across the segment interval
x_eval = linspace(x_fit(1), x_fit(end), 1000);
p_eval = zeros(size(x_eval));

tol = 1e-14;

for i = 1:length(x_eval)
    x = x_eval(i);
    d = x - x_fit;
    
    % Removable singularity check
    exact_match = find(abs(d) < tol, 1);
    if ~isempty(exact_match)
        p_eval(i) = f_fit(exact_match);
    else
        terms = w ./ d;
        p_eval(i) = sum(terms .* f_fit) / sum(terms);
    end
end

% Plotting results
figure;
plot(x_eval, p_eval, 'b-', 'LineWidth', 1.5); hold on;
plot(x_fit, f_fit, 'ro', 'MarkerSize', 5, 'LineWidth', 1.2);
xlabel('x');
ylabel('p(x)');
title('Second Barycentric Polynomial Interpolation (Verified Segment)');
legend('Interpolating Polynomial', 'Data Points', 'Location', 'best');
grid on;