function [f,Sd] = my_8qam_decisor_circular(Sr)
Sd = [];
f =[];
%symbol vector
symbols = [(1+sqrt(3))+0j,1-1j,-1-1j,0-(1+sqrt(3))*1i,1+1j,0+(1+sqrt(3))*1i,-(1+sqrt(3))+0j,-1+1j];
Q = log2(length(symbols)); % number of bits per symbol
bits = de2bi(0:length(symbols)-1,Q,'left-msb');% code book of bits
for xhat=Sr
    [~,idx] = min(abs(xhat*ones(1,length(symbols))-symbols).^2,[],2);%find the minimum distance point and its index
    Sd  = [Sd symbols(idx)];%store to the putput symbols
    f = [f bits(idx,:)];%store to the output bi stream
end
end