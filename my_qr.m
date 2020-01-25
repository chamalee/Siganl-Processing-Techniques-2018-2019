function [Q,R] = my_qr(A,economy)

[m,n]=size(A); %get the input matrix dimentions
Q = eye(m); %let the initial Q matrix to be identity matrix
R=A;%let the initial R matrix to be A then A = QR

stnd_basis = eye(m,m);%standard basis of the given space of dimention m

for j=1:n %goes through coloumns 
    x = R(1:end,j);%take x vector as the jth coloumn of R
    basis = stnd_basis(:,j);%standard basis correspond to the jth coloumn
    
    %except for the first vector make elements below jth row zero
    if(j~=1) 
        x(1:j-1) = 0;      
    end
    %factor 
    alpha = -exp(1i*angle(x(j,1)))*norm(x);%vector norm
    
    u = x - alpha*basis;%householder vector
    
    qq = eye(m,m) - 2*(u*u')/norm(u)^2;%house holder matrix
    
    Q = qq*Q; % new Q matrix by pre multiplying by the householder matrix
    R = qq*R; % new R matrix transformed by the householder matrix
    
    %makes very small values in the matrix zero
    Q(abs(Q)<1e-10) = 0;
    R(abs(R)<1e-10) = 0;
    
    if(istriu(R)) %if already a upper triangular matrix exit the loop
        break;
    end
end
if(economy)
    [m,n] =size(A);
    if(m>n)
        Q = Q';
        R(n+1:m,:) = [];
        Q(:,n+1:m)=[];
    end
else
    Q = Q';
    R = R;
end


end

