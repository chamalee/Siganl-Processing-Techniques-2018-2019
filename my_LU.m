function [L,U] = my_LU(A)

[m,n] = size(A);
L = eye(size(A));
U = zeros(size(A));

for j=1:n
    for i=1:j
        alpha = A(i,j);
        if(i==1)
           alpha = alpha;
        else
            for k=1:i-1
                alpha = alpha - L(i,k)*U(k,j);
            end
        end
            U(i,j) =alpha;       
    end
    
    for i=j+1:n
        alpha = A(i,j);
        for k=1:j-1
            alpha = alpha - L(i,k)*U(k,j);
        end
        L(i,j) = alpha/U(j,j);
    end
end


end

