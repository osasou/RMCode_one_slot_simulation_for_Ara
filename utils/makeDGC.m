function K = makeDGC(re,m)

[i j] = upper_indices(m,re);
%disp([i,j])
%makeDGC(1,3)�̏ꍇ�A
%i=[1 1 2]
%j=[2 3 3]�Ɠ���

for p = 1:length(i)
    
    K(:,:,p) = zeros(m,m);
    % i�sj���j�si���1��������
    % K_1�ɂ́A1�s2���2�s1��
    % K_2�ɂ́A1�s3���3�s1��
    % K_3�ɂ́A2�s3���3�s2��
    K(i(p),j(p),p) = 1;
    K(j(p),i(p),p) = 1;
    
end

%disp(K)