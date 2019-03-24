function [A] = generate_A(multiplier)
A=zeros(20,20);
for i=1:1:20
    for j=1:1:20
        if i==j
            A(i,j)=3*multiplier;
            continue;
        end
        if abs(i-j)==1
            A(i,j)=-1/2;
            continue;
        end
        if abs(i-j)==2
            A(i,j)=-1/4;
            continue;
        end
    end
end
end

