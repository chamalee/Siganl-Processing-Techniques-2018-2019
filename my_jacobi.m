function [J1,D,J2] = my_jacobi(A)
W = A;
[m,n] = size(W);
%% %Transform to a real matrix
if(~isreal(A))
    [theta] = angle(W);%angles of all the complex values of the matrix
    [rho] = abs(W);
    theta_1 = theta(1,1); %angle of the complex value at (1,1)
    theta = theta - theta_1;%make a11 a real value
    
    B_ = [rho(1,1) rho(1,2)*exp(1i*theta(1,2));rho(2,1)*exp(1i*theta(2,1)) rho(2,2)*exp(1i*theta(2,2))];%new complex matrix
    
    %triangularization
    c0 = abs(B_(1,1))/(sqrt(abs(B_(1,1))^2+abs(B_(2,1))^2));
    s0 = (B_(2,1)/B_(1,1))*c0;
    
    Q0 = [c0 conj(s0);-s0 c0];
    
    X_= Q0*B_;
    
    %phase cancellation
    phase_l = [exp(-1i*angle(X_(1,2))) 0;0 exp(-1i*angle(X_(2,2)))];%left side phase cancellation
    phase_r = [exp(1i*angle(X_(1,2))) 0;0 1];%right side phase cancellation
    
    X = phase_l*X_*phase_r;
else
    X = A;
    phase_l = eye(2,2);
    phase_r=eye(2,2);
    Q0 = eye(2,2);
    theta_1  =0;
end
%% %Transform to a symmetric matrix
X(abs(X)<1e-10) = 0; %make very small values zero
WW = real(X);%X is a real matrix now
W = X;
if(W(1,2)~=W(2,1) && WW(2,1)== 0 && WW(2,2)==0)%special case 1
    rho = (W(1,1)+abs(W(2,2)))/(W(1,2)-W(2,1));
    
    s1 = sign(rho)/sqrt(1+rho^2);
    c1 = rho*s1;
    
    Q1 = [c1 s1;-s1 c1];
    Y = W*Q1';
elseif(W(1,2)~=W(2,1))%non symmetric case
    rho = (W(1,1)+abs(W(2,2)))/(W(1,2)-W(2,1));
    
    s1 = sign(rho)/sqrt(1+rho^2);
    c1 = rho*s1;
    
    Q1 = [c1 s1;-s1 c1];
    Y = Q1'*W;
else
    Y = W; %if already symmetric
    Q1 = eye(m,m);
end

%% % SVD decmposition
yy=Y;
if((WW(1,2)== 0 && WW(2,2)==0))%[a1 0;a2 0] case 
    D =real(yy);%for the better reading of the matrix
    J1 = exp(1i*theta_1)*(Q0'*phase_l'*Q1');%pre multiplication matrix
    J2 = (phase_r);%pre multiplication matrix
    
elseif((WW(2,1)== 0 && WW(2,2)==0))%special case 2: A = [a1 a2 ; 0 0]
    D =real(yy);%for the better reading of the matrix
    J1 = exp(1i*theta_1)*(Q0'*phase_l');
    J2 = (phase_r*Q1');
    
elseif((yy(2,1)== 0 && yy(1,2)==0))%resultant matrix already in diagonal form
    D =real(yy);%for the better reading of the matrix
    J1 = exp(1i*theta_1)*(Q0'*phase_l');%pre multiplication matrix
    J2 = (phase_r);  %post multiplication matrix
    
else %SVD deomposition of the symmetric matrix
    Q2 =eye(2,2);
    
    e= (yy(2,2) -yy(1,1))/(2*yy(2,1));
    if(e~=0)
        t = -sign(e)*((abs(e)+sqrt(1+e^2)));
    else
        t = 1;
    end
    c2 = 1/sqrt(1+t^2);
    s2 = t*c2;
    
    qq2 = [c2 s2;-s2 c2];
    yy = qq2'*yy*qq2;
    
    Q2 = Q2*qq2;
    D =real(yy);%for the better reading of the matrix
    J1 = exp(1i*theta_1)*(Q0'*phase_l'*Q1*Q2);%overall left transformation matrix
    J2 = (phase_r*Q2);%overall right transformation matrix
end
end