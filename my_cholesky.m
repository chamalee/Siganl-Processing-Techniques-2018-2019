function [L] =  my_cholesky(A)

[m,n]=size(A); %get the input matrix dimentions
L = eye(m); 
for i=1:n
  
         for j=1:i
             if(i == j) % diagonal
                 sum = 0;
                 
                 for k=1:j-1
                     sum = sum + L(j,k)*conj(L(j,k));
                 end
                 L(j,j) = sqrt(A(j,j)-sum);
                 
     
             else % non diagonal
                  sum = 0;
                 
                 for k=1:j-1
                     sum = sum + L(i,k)*conj(L(j,k));
                 end
                  L(i,j) = (A(i,j)-sum)/L(j,j);
         end
             end
        end
end
  
  
    
