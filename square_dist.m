function dist = square_dist(h_all,h_hat)

sum = 0;

% disp([h_all,h_hat]);
% sorted_h_all = sort(h_all);
% sorted_h_hat = sort(h_hat);

sorted_h_all = sort(h_all,'ComparisonMethod','real');
sorted_h_hat = sort(h_hat,'ComparisonMethod','real');

% disp('square_dist');
% 
% disp([sorted_h_all,sorted_h_hat]);

for i=1:length(h_all)
        sum = sum + (real(sorted_h_all(i)) - real(sorted_h_hat(i)))^2 + (imag(sorted_h_all(i)) - imag(sorted_h_hat(i)))^2;
end

dist = sum/length(h_all);

