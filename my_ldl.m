function [L,D] =  my_ldl(A)

[m,n]=size(A); %get the input matrix dimentions
L = eye(m); 
D = eye(m); 
for i=1:n
  
         for j=1:i
             if(i == j) % diagonal
                 sum = 0;
                 
                 for k=1:j-1
                     sum = sum +( L(j,k)*conj(L(j,k))*D(k,k));
                 end
                 D(j,j) = A(j,j)-sum;
                 
     
             else % non diagonal
                  sum = 0;
                 
                 for k=1:j-1
                     sum = sum + (L(i,k)*conj(L(j,k))*D(k,k));
                 end
                  L(i,j) = (A(i,j)-sum)/D(j,j);
         end
             end
        end
end
  
  
    
