%% Matlab implementation of OFDM for one transmission block.
% Transmitter and receiver have one antenna.
% - 159744 bits of input data.
% - 8-QAM symbols.
% - 128 subcarriers.
% - Use of Cyclic Prefix.
% - AWGN
% SNR = 10 dBm
clear; clc;close all;
rng(12345);
%% initialization  
N=159744;% #bits
SC=128; % #subcarriers
CP=16; % #cycle prefix
SNR=10 ;%-15:2.5:15; %SNR
z=zeros(1,length(SNR));
L = 5; % number taps of the channels
k=1;
M=3;
%% L-tap channel coefficients
h = randn(1,L);% + 1i*randn(1,L); %cannel 3
smat=zeros(SC,(N/M)/SC);
smatifft=zeros(SC,(N/M)/SC);
xt=zeros(1,SC);
xr=zeros(1,SC);
t=zeros(SC,(N/M)/SC);
xscat=zeros(SC,(N/M)/SC);
xrfftscat=zeros(SC,(N/M)/SC);
rmat=zeros(SC,(N/M)/SC);
xtx=zeros(1,N/M);
ts=zeros(1,N/M);
vectcp=zeros(1,SC+CP);
smatcp=zeros(SC+CP,(N/M)/SC);
vectrcp=zeros(SC,1);
rmatcp=zeros(SC+CP,(N/M)/SC);
I=eye(SC);
W=zeros(1,SC);

%create input chain of bits
x=randi(2,1,N);
x=x-1;

%convert bits into symbols
s = my_8qam_circular(x);
for d=1:N/M
    smat(d)=s(d);
end

%ifft
for e=1:((N/M)/SC)
    xt=smat(:,e);
    xtifft=ifft(xt);
    smatifft(:,e)=xtifft;
end
%cyclic prefix
for e=1:((N/M)/SC)
    smatcp(1:CP,e)=smatifft(SC-CP+1:SC,e);
    smatcp(1+CP:SC+CP,e)=smatifft(1:SC,e);
end
%paralel to serial
for d=1:((N/M)+CP*((N/M)/SC))
    xtx(d)=smatcp(d);
end
%variance tx
potx=0;
for m=1:((N/M)+CP*((N/M)/SC))
    potx=potx + (abs(xtx(1,m)))^2;
end
px=potx/((N/M)+CP*((N/M)/SC));
%channel
xtxch=conv(xtx,h);
trch=zeros(1,(N/M)+CP*((N/M)/SC));
for d=1:((N/M)+CP*((N/M)/SC))
    trch(d)=xtxch(d);
end
%variance after the channel to perform the AWGN
pot=0;
for m=1:((N/M)+CP*((N/M)/SC))
    pot=pot + (abs(trch(1,m)))^2;
end
for isnr=SNR
ps=pot/((N/M)+CP*((N/M)/SC));
var=ps/(8*(10^(isnr/10)));

%noise
nr=randn(1,((N/M)+CP*((N/M)/SC)));
nr=sqrt(var)*nr;
ni=randn(1,((N/M)+CP*((N/M)/SC)));
ni=sqrt(var)*ni;
ni= 1i*ni;
n=nr+ni;
%variance of the noise
potn=0;
for m=1:((N/M)+CP*((N/M)/SC))
    potn=potn + (abs(n(1,m))^2);
end
pn=potn/((N/M)+CP*((N/M)/SC));
r = trch + n ;
%serial to paralel
for d=1:((N/M)+CP*((N/M)/SC))
rmatcp(d)=r(d);
end

%remove cyclic prefix
for e=1:((N/M)/SC)
rmat(:,e)=rmatcp(1+CP:SC+CP,e);
end

%coeficients for  equalizer
H=fft(h,128);

Ckmmse=zeros(1,SC);

for d=1:SC
    Ckmmse(d)=   (conj(H(d)))/((abs(H(d))^2)+pn/px);
    Ckzf(d) = 1/H(d);
end

%fft
MMSE =1;
for e=1:((N/M)/SC)
    xr=rmat(:,e);
    xrfft=fft(xr);
    xrfftscat(:,e)=xrfft;
    for d=1:SC
        if(MMSE == 1)
            xrfft(d)= xrfft(d)*(Ckmmse(d)) ;% or Ckzf(d)); %equalization
        else
            xrfft(d)= xrfft(d)*Ckzf(d); %equalization
        end
    end
    xscat(:,e)=xrfft;
    t(:,e)=xrfft;
end

for d=1:N/M
    ts(d)=t(d);
end

f=zeros(1,N);
fr=real(ts);
fi=imag(ts);

%decisor
%convert symbols into bits
[f,Sd] = my_8qam_decisor_circular(ts);

%BER
d=1;
count=0;
for d=1:N
    if x(d)==f(d)
    else count=count+1;
    end
end
BER=count/N;
z(k)=BER;
k=k+1;
end
%theoretical BER
yid=SNR;
z_id = (4/3)*qfunc(sqrt(9*10.^(yid/10)/7));

%BER plots

y=1:1:length(SNR);
semilogy(yid,z_id, 'r-')
hold on
semilogy(yid,z, 'b*')
title('SNR vs BER')
grid on;
hold off;

figure;
semilogy(10*log(abs(H)));
title('Channel Gains (dB)')

if(MMSE ==1)
%plot inputs of different sub carriers
figure;
subplot(421);plot(xrfftscat(1,:),'.');title('MMSE Input, Sub Carrier  = 1');grid on;
subplot(422);plot(xrfftscat(8,:),'.');title('MMSE Input, Sub Carrier  = 8');grid on;
subplot(423);plot(xrfftscat(24,:),'.');title('MMSE Input, Sub Carrier  = 24');grid on;
subplot(424);plot(xrfftscat(46,:),'.');title('MMSE Input, Sub Carrier  = 46');grid on;
subplot(425);plot(xrfftscat(66,:),'.');title('MMSE Input, Sub Carrier  = 66');grid on;
subplot(426);plot(xrfftscat(86,:),'.');title('MMSE Input, Sub Carrier  = 86');grid on;
subplot(427);plot(xrfftscat(106,:),'.');title('MMSE Input, Sub Carrier  = 106');grid on;
subplot(428);plot(xrfftscat(118,:),'.');title('MMSE Input, Sub Carrier  = 118');grid on;
grid on;

figure;
subplot(421);plot(xscat(1,:),'.');title('MMSE Output, Sub Carrier  = 1');grid on;
subplot(422);plot(xscat(8,:),'.');title('MMSE Output, Sub Carrier  = 8');grid on;
subplot(423);plot(xscat(24,:),'.');title('MMSE Output, Sub Carrier  = 24');grid on;
subplot(424);plot(xscat(46,:),'.');title('MMSE Output, Sub Carrier  = 46');grid on;
subplot(425);plot(xscat(66,:),'.');title('MMSE Output, Sub Carrier  = 66');grid on;
subplot(426);plot(xscat(86,:),'.');title('MMSE Output, Sub Carrier  = 86');grid on;
subplot(427);plot(xscat(106,:),'.');title('MMSE Output, Sub Carrier  = 106');grid on;
subplot(428);plot(xscat(118,:),'.');title('MMSE Output, Sub Carrier  = 118');grid on;
grid on;
else
%plot inputs of different sub carriers
figure;
subplot(421);plot(xrfftscat(1,:),'.');title('ZF Input, Sub Carrier  = 1');grid on;
subplot(422);plot(xrfftscat(8,:),'.');title('ZF Input, Sub Carrier  = 8');grid on;
subplot(423);plot(xrfftscat(24,:),'.');title('ZF Input, Sub Carrier  = 24');grid on;
subplot(424);plot(xrfftscat(46,:),'.');title('ZF Input, Sub Carrier  = 46');grid on;
subplot(425);plot(xrfftscat(66,:),'.');title('ZF Input, Sub Carrier  = 66');grid on;
subplot(426);plot(xrfftscat(86,:),'.');title('ZF Input, Sub Carrier  = 86');grid on;
subplot(427);plot(xrfftscat(106,:),'.');title('ZF Input, Sub Carrier  = 106');grid on;
subplot(428);plot(xrfftscat(118,:),'.');title('ZF Input, Sub Carrier  = 118');grid on;

%plot outputs of different sub carriers
figure;
subplot(421);plot(xscat(1,:),'.');title('ZF Output, Sub Carrier  = 1');grid on;
subplot(422);plot(xscat(8,:),'.');title('ZF Output, Sub Carrier  = 8');grid on;
subplot(423);plot(xscat(24,:),'.');title('ZF Output, Sub Carrier  = 24');grid on;
subplot(424);plot(xscat(46,:),'.');title('ZF Output, Sub Carrier  = 46');grid on;
subplot(425);plot(xscat(66,:),'.');title('ZF Output, Sub Carrier  = 66');grid on;
subplot(426);plot(xscat(86,:),'.');title('ZF Output, Sub Carrier  = 86');grid on;
subplot(427);plot(xscat(106,:),'.');title('ZF Output, Sub Carrier  = 106');grid on;
subplot(428);plot(xscat(118,:),'.');title('ZF Output, Sub Carrier  = 118');grid on;
end