function [xhat] = my_MMSE (H,y,N0)

[m, n] = size(H);
 % Gramian 
HH = H'*H;
% MMSE detector
xhat = inv(HH +N0*eye(n))*H'*y;

end

