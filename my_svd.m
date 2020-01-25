function [U_O,D_O,V_O] = my_svd(A)

[m,n] = size(A);
%when A is a single value
if(m==1 && n==1)
    D_O = abs(A);
    U_O = A/abs(A);
    V_O = 1;
    return;
end

if(m > n) 
    [Q,B]=my_qr(A); %take the QR decomposition for tall matrices
elseif(m < n)
    B = [A;zeros(n-m,n)];
else
    B =A;
end
U = eye(m,m);
V = eye(n,n);

N_2 =0;
for i=1:m
    for j=1:n
        N_2 = N_2+abs(B(i,j))^2;
    end
end
e= 1e-10;
s_er = N_2;
first = true;
while((s_er > e*N_2)|| first == true)
    
    s_er=0;
    first = false;
    for i = 1: m-1
        for j = i+1 : n 
            b_sub =  [B(i,i) B(i,j) ; B(j,i) B(j,j)];
            [J11,D,J22] = my_jacobi(b_sub);%jacobian algorithm for 2x2 matrix
            
            %cyclic jacobian matrix left multiplication
            J1 = eye(m,m);
            
            J1(i,i) = J11(1,1);
            J1(i,j) = J11(1,2);
            J1(j,i) = J11(2,1);
            J1(j,j) = J11(2,2);
            
            %cyclic jacobian matrix right multiplication
            J2 = eye(n,n);
            
            J2(i,i) = J22(1,1);
            J2(i,j) = J22(1,2);
            J2(j,i) = J22(2,1);
            J2(j,j) = J22(2,2);
            
            B = J1'*B*J2;%resultant matrix
            U = U*J1;%left transformation
            V = V*J2;%right transformation
        end
    end
    for i=1:m-1
        for j=i+1:n
            s_er = s_er + abs(B(i,j))^2+abs(B(j,i))^2;%error calculation
        end
    end
end
B(abs(B)<1e-4) = 0;
D=real(B);
if(m > n) %tall matrix case
    U_O = Q*U(1:m,1:m);
    V_O = V(1:n,1:n);
    D_O = D(1:m,1:n);
elseif(m < n)%Fat matrix case
    U_O = U(1:m,1:m);
    V_O = V(1:n,1:n);
    D_O = D(1:m,1:n);
else % for a square matrix
    U_O = U;
    V_O = V;
    D_O = D;
end


end

