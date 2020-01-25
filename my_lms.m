function [c,e,X] = my_lms(mu,y,v,M,L,T,EqD)

c = zeros(L+1,1);
X=zeros(L+1,M); 
y = [y zeros(1,L+1)];

e=zeros(1,T-EqD); 

%windowing the signal to be multiplied by estimated wiener coefficients
    for i=1:M
        ytemp = y(L+i:-1:i).';
        X(:,i)=ytemp;
    end

    D = zeros(size(v));%new desired sequence initialization
    for i=EqD:M-EqD
        D(1,i) = v(1,i-EqD+1);%delaying the input signal by EqD
    end

    for i=1:T-EqD
       
        e(i)=D(i+EqD-1)-c'*X(:,i); % instant error with the delayed desired sequence
        c=c+mu*conj(e(i))*X(:,i); % update filter or equalizer coefficient
    end
end




