function [ Sd ] = my_8qam_circular( St )
%symbol vector
symbols = [(1+sqrt(3))+0j,1-1j,-1-1j,0-(1+sqrt(3))*1i,1+1j,0+(1+sqrt(3))*1i,-(1+sqrt(3))+0j,-1+1j];
Q = log2(length(symbols)); % number of bits per symbol
k=1;
%Goes through the input bitstream
for m=1:Q:length(St)
    idx = bi2de(St(1,m:m+Q-1),'left-msb')+1; %decimal value of teh three bits selected
    Sd(k) = symbols(idx).'; %Select symbol from the symbol matrix
    k=k+1;
end
end
