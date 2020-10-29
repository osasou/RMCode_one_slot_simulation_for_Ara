% upper_indices    Outputs a list of entries in the upper half of a matrix.
%                  Complex includes the diagonal, real does not.
%
% m          size of m x m matrix    
% re         logical: false=complex chirps, true=real chirps
%
%[i j]       [row col] coordinates
%
% AJT (12/9/18)

function [i j] = upper_indices(M,re)

inds = 0:M^2-1;
j = mod(inds,M) + 1;
i = (inds - j + 1)/M + 1;
%{
disp(inds)
disp(j)
disp(i)
%}
if (re==0)
    upper = i<=j;
else
    upper = i<j;
end
%{
disp(upper)
disp(i(upper))
disp(j(upper))
%}
%i(upper)�͔z��upper�ł�1�̕����̗v�f�̏ꏊ�ɓ�����z��i�����o�����z��B
i = i(upper);
j = j(upper);