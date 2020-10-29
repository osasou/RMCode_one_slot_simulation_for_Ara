function [h_value] = rayleigh_rand(num)
%RAYLY ���̊֐��̊T�v�������ɋL�q
%   �ڍא����������ɋL�q
h_value = [];
x = [];
for i = 1:num
    x = [x, i];
    r1 = normrnd(0,0.5);
    r2 = normrnd(0,0.5);
    h = sqrt(r1^2 + r2^2);
    h_value = [h_value,h];
end

%scatter(x, h_value)
%histogram(h_value,100);

%{
�������ł�����
h = [];
h = raylrnd(0.5,1000000,1);
histogram(h,100)
%}

end

