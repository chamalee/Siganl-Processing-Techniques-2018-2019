clear all;
rng(123)
clc;
%rng(123)
M=3000;    % number of data samples
T=2000;    % number of training symbols
dB=25;     % SNR in dB value

L=20; % length for smoothing(L+1)

LMS = 1;
NLMS = 2;
LMS_TYPE = LMS;

%%%%% QPSK TRANSMISSION %%%%%%%
TxS=round(rand(1,M))*2-1;  
TxS=TxS+sqrt(-1)*(round(rand(1,M))*2-1);


%%%%%%%% Channel %%%%%%%%%%%%%%%%
ChL=5;  % length of the channel(ChL+1)
Ch=randn(1,ChL+1)+sqrt(-1)*randn(1,ChL+1);   
Ch=Ch/norm(Ch);   

% signal filtered by channel
x= filter(Ch,1,TxS);



EqD= round((L+ChL)/2);  %delay for equalization

%%%%%%%%% NOISE %%%%%%%%%%%%%%%%
n=randn(1,M);  
n=n/norm(n)*10^(-dB/20)*norm(x);


%%%%%%%% RECEIVED SIGNAL %%%%%%%%%%%%%%%
y=x+n;                          

%%%%%%%%%%%%%% WRITE LMS ALGORITHM FOR EQUALIZATION HERE %%%%%%%%%%%%%%%%%%
muLMS = 0.01;%mu values found by experiment
muNLMS = 0.5;
if(LMS_TYPE == LMS)
    [c,e,X] = my_lms(muLMS,y,TxS,M,L,T,EqD);
else
    [c,e,X] = my_nlms(muNLMS,y,TxS,M,L,T,EqD);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sb=c'*X;  % recieved symbol estimation
% sb1=sb;
%SER(decision part)
sb1=sb/norm(c);  % normalize the output
sb1=sign(real(sb1))+sqrt(-1)*sign(imag(sb1));  %symbol detection
sb2=sb1-TxS(1:length(sb1));  % error detection
SER=length(find(sb2~=0))/length(sb2); %  SER calculation
disp(SER);
% 
% plot of transmitted symbols
    subplot(2,2,1), 
    plot(TxS,'*');   
    grid,title('Input symbols');  xlabel('real part'),ylabel('imaginary part')
    axis([-2 2 -2 2])
    
% plot of received symbols
    subplot(2,2,2),
    plot(y,'o');
    grid, title('Received samples');  xlabel('real part'), ylabel('imaginary part')
    
% plots of the equalized symbols    
    subplot(2,2,3),
    hold on;
    plot(sb,'o');   
    grid, title('Equalized symbols'), xlabel('real part'), ylabel('imaginary part')

% convergence
    subplot(2,2,4),
    semilogy(abs(e));   
    grid, title('Convergence'), xlabel('n'), ylabel('error signal')