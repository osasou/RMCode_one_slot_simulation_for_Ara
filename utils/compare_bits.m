function propfound = compare_bits(input,output)

%propfound = input�̗�Ō��āAoutput�ɂ�������count+1���Ă����Acount/size(input,2)  (input�̗�)
k1 = size(input,2);
k2 = size(output,2);

found = zeros(k1,1);
%{
disp("k1")
disp(k1)

disp("k2")
disp(k2)

disp("input")
disp(input)
disp("output")
disp(output)
%}
for i = 1:k1
    %{
    disp("input")
    disp(i)
    disp(input(:,i))
    %}
    for j = 1:k2
        %{
        disp("output")
        disp(j)
        disp(output(:,j))
        %}
        try
            if any(output(:,j)~=input(:,i)) % ~= means != in C lang
                continue
            end
            found(i,1) = j;

            break
        catch
            disp("input and output size errorrrr in compare_bits.m")
            disp("input size  output size")
            disp([length(input'), length(output')])
        end
    end
end
if(k1 < sum(found>0))
    disp("warninggggg")
end
propfound = sum(found>0)/k1; 