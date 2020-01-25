clc;
clear all;
figure;
channel = [1 2 4];
tic;
for n_Trans = channel
trials = 100;  % Number of Monte-Carlo Simulations
nR = n_Trans;        % Number of Received Antenna
nT = n_Trans;        % Number of Transmit Antenna
SNR_dB = [-10:5:40];  % SNR values in dB

sumRate = zeros(1,length(SNR_dB));

snrIndex = 0;
for iSNR = SNR_dB
    
    xRate = 0;
    snrIndex = snrIndex + 1;
    txPower = 10^(iSNR / 10);  %% linear SNR
    
    for i = 1:trials
        
        H = complex(randn(nR,nT),randn(nR,nT)) * sqrt(0.5);  % Rayleigh Fading Channel  
        
        %%%%%%%%WRITE YOUR OWN SVD FUNCTION HERE%%%%%%%%%%
        [my_U,my_D,my_V] = my_svd(H); % my SVD function
        [U,D,V] = svd(H); % Matlab built-in SVD function
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        txPrecoder = V * sqrt(txPower / nT);
        rxBeamformer = U';
        %% caculate the data rate
        ant =min(nR,nT);
        for j = 1: ant
        xRate = xRate + log2(1+(D(j,j)^2)*((10^(iSNR / 10))/nT)); %Summing up rates
        end       
    end
    
    hold all;
    sumRate(1,snrIndex) = xRate / trials;
    
end


 plot(SNR_dB,sumRate,'-+');

hold on
legend('1x1 SISO','2x2 MIMO','4x4 MIMO')
xlabel('SNR [dB]')
ylabel('Channel Capacity (bits/Hz)')
title('Channel Capacity Vs SNR')

end
toc