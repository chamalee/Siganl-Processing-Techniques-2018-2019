%% Matlab implementation of OFDM for one transmission block.
% Transmitter and receiver have one antenna.
% - 159744 bits of input data.
% - 4-QAM symbols.
% - 64 subcarriers.
% - Use of Cyclic Prefix.
% - AWGN
% SNR = 10 dBm
clear; clc;
close all;
%% initialization  
kk=1;
M=3;
Nr = 2;
Nt = 2;

N=159744*Nt;% #bits
SC=128; % #subcarriers
CP=16; % #cycle prefix
SNR=-15:5:15; %SNR
z=zeros(1,length(SNR));
L = 5; % number taps of the channels

%% L-tap channel coefficients
H = randn(Nr,Nt) + 1i*randn(Nr,Nt); %cannel 3
h = rand(1,L);
smat_=zeros(SC,(N/M)/SC);

smat=zeros(SC,(N/M)/SC/Nt,Nt);
smatifft=zeros(SC,(N/M)/SC/Nt,Nt);
xt=zeros(1,SC);
xr=zeros(1,SC);
t=zeros(SC,(N/M)/SC);
xscat=zeros(SC,(N/M)/SC/Nt);
xrfftscat=zeros(SC,(N/M)/SC/Nt);
rmat=zeros(SC,(N/M)/SC/Nt,Nr);
rmat_=zeros(SC,(N/M)/SC);
ts=zeros(1,N/M);
vectcp=zeros(1,SC+CP);
smatcp=zeros(SC+CP,(N/M)/SC/Nt,Nt);
vectrcp=zeros(SC,1);
rmatcp=zeros(SC+CP,(N/M)/SC/Nt,Nr);
I=eye(SC);
W=zeros(1,SC);

%create input chain of bits
x=randi(2,1,N);
x=x-1;

%convert bits into symbols
%s = my_8qam(x);
s = my_8qam_circular(x);

for d=1:N/M
   smat_(d)=s(d);
end
for m=1:Nt
    smat(:,:,m) = smat_(:,(m-1)*(N/M)/SC/Nt+1:m*(N/M)/SC/Nt);
end

%ifft
for m=1:Nt
    for e=1:((N/M)/SC/Nt)
        xt=smat(:,e);
        xtifft=ifft(xt);
        smatifft(:,e,m)=xtifft;
    end
end
%cyclic prefix
for m=1:Nt
    for e=1:((N/M)/SC/Nt)       
        smatcp(1:CP,e,m)=smatifft(SC-CP+1:SC,e,m);
        smatcp(1+CP:SC+CP,e,m)=smatifft(1:SC,e,m);        
    end
end
%paralel to serial
xtx=zeros(Nt,((N/M/Nt)+CP*((N/M)/SC/Nt)));
for m=1:Nt
    tp_smatcp = reshape(smatcp(:,:,m),[],1) ;
    xtx(m,:) = tp_smatcp';% convert matrix to row vector
end

%variance tx
potx=0;
for k =1:Nt
    for m=1:((N/M/Nt)+CP*((N/M)/SC/Nt))
        potx=potx + (abs(xtx(k,m)))^2;
    end
end

px=potx/((N/M/Nt)+CP*((N/M/Nt)/SC));

%channel
%xtxch=conv(xtx,h);
xtxch=H*xtx;

trch=zeros(Nr,(N/M/Nt)+CP*((N/M)/SC/Nt));
for m=1:Nr
    for d=1:((N/M/Nt)+CP*((N/M)/SC/Nt))
        trch(m,d)=xtxch(m,d);
    end
end
%variance after the channel to perform the AWGN
pot=0;
for k=1:Nr
    for m=1:((N/M/Nt)+CP*((N/M)/SC/Nt))
        pot=pot + (abs(trch(k,m)))^2;
    end
end
ps=pot/((N/M)+CP*((N/M)/SC));

for isnr=SNR

var=ps/(8*10^(isnr/10));
%noise
nr=randn(Nr,((N/M/Nt)+CP*((N/M)/SC/Nt)));
nr=sqrt(var)*nr;
ni=randn(Nr,((N/M/Nt)+CP*((N/M)/SC/Nt)));
ni=sqrt(var)*ni;
ni= 1i*ni;
n=nr+ni;
%variance of the noise
potn=0;
for k=1:Nr
for m=1:((N/M/Nt)+CP*((N/M)/SC/Nt))
    potn=potn + (abs(n(Nr,m))^2);
end
end
pn=potn/((N/M)+CP*((N/M)/SC));

r = trch + n ;
%apply MMSE equalization at receiver antennas
xhat = my_MMSE(H,r,pn^2);

%serial to paralel
for m = 1:Nt 
rmatcp(:,:,m) = reshape(xhat(m,:),[SC+CP (N/M/Nt/SC)]);
end
% for d=1:((N/M/Nt)+CP*((N/M)/SC/Nt))
% rmatcp(d)=r(d);
% end
% for m=1:Nr
%     rmat(:,:,m) = rmat_(:,(m-1)*(N/M)/SC/Nr+1:m*(N/M)/SC/Nr);
% end

%remove cyclic prefix
for m=1:Nt
    for e=1:((N/M)/SC/Nt)
        rmat(:,e,m)=rmatcp(1+CP:SC+CP,e,m);
    end
end

%coeficients for  equalizer
HH=fft(h,128);
Ckmmse=zeros(2,SC);

for m=1:Nt
for d=1:SC
    Ckmmse(d)= conj(HH(d))/((abs(HH(d))^2)+pn/px);
    Ckzf(d) = 1/HH(d);
end
end

%fft
MMSE =0;
for m=1:Nt
for e=1:((N/M)/SC/Nt)
    xr=rmat(:,e,m);
    xrfft=fft(xr);
    xrfftscat(:,e,m)=xrfft;
%     for d=1:SC
%         if(MMSE == 1)
%             xrfft(d)= xrfft(d)*(Ckmmse(d)) ;% or Ckzf(d)); %equalization
%         else
%             xrfft(d)= xrfft(d)*Ckzf(d); %equalization
%         end
%     end
    xscat(:,e,m)=xrfft;
    t(:,e)=xrfft;
end
end
tss =[];
for m=1:Nt
    tss = [tss reshape(xscat(:,:,m),[1 (N/M)/Nt])];
end
for d=1:N/M
    ts(d)=t(d);
end

f=zeros(1,N);


%decisor
%convert symbols into bits
[f,Sd] = my_8qam_decisor_circular(tss);
%BER
d=1;
count=0;
for d=1:N
    if x(d)==f(d)
    else count=count+1;
    end
end
BER=count/N
z(kk)=BER;
kk=kk+1;
end
%theoretical BER
yid=SNR;
%z_id=qfunc(sqrt(2*10.^(yid/10)));
z_id =(4/log2(M))*qfunc(sqrt(3*log2(M)*10.^(yid/10))/(M-1));
%BER plots
figure;
y=1:1:length(SNR);
semilogy(yid,z_id, 'r-')
hold on
semilogy(yid,z, 'b*')
hold off;
l=80
figure;
plot(abs(H));
figure;
subplot(131);plot(s,'+')
subplot(132);plot(r(1,1:1000),'o')
subplot(133);plot(xscat(l,:),'*')

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