function dft_out = my_fourier(x_n)
%% %Definition of variables
N = length(x_n);%length of the input signal
dft_out = zeros(N,1);%fourier transform output
%% % DFT calculation
for n=1:N
    for k=1:N
    dft_out(n) = dft_out(n)+ x_n(k)*exp(-1i*2*pi*(n-1)*(k-1)/N);% from the DFT equation
    end
end 

end

