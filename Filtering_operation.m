%% Task 3
clear;
clc; close all;
compFiltering = true;

N=100;
L=3; %number of taps

%% %Comparison of filtering operations
if(compFiltering == true)
    X_n = rand(N,1);
    a_n = [1 rand(1,L)]; % filter coefficents denominator
    b_n =[1]%filter coefficients nominator
    tic;
    y_n_1 = my_filter2(X_n,a_n);%IIR filter implementation
    myfilt_t = toc %time taken by my filter function

    tic;
    y_n_2 = filter(b_n,a_n,X_n); % FIR filter using Matlab built-in function
    matlab_t = toc % time taken by matlab filter function
    %Genrate figures
    subplot(311);stem([1:N]-1,X_n,'fill','r');title('input signal, V1[n]=');grid on;
    subplot(312);stem([1:length(y_n_1)]-1,y_n_1,'fill','r');title('filter output (my filter), Y[n]=');grid on;
    subplot(313);stem([1:length(y_n_2)]-1,y_n_2,'fill','r');title('filter output (built in), Y[n]=');grid on;
    xlabel('time n');
end

%% Generation of X[n]

%Desired signal generation
b_n1 = [1];
a_n1 = [1 0.8458];

%channel as a filter
b_n2 = [1];
a_n2 = [1 -0.9458];

%variances of noise signals
var_v1 = 0.27;
var_v2 =0.1;

%generate the noise signal vectors
v1_n = randn(N,1).*sqrt(var_v1);  % generate v1
v2_n = randn(N,1).*sqrt(var_v2); %generate v2

d_n = my_filter2(v1_n,a_n1);%desired signal from IIR filter output
x_n = my_filter2(d_n,a_n2);%channel output from channel as a IIR filter

u_n = x_n + v2_n;%received signal with noise

%wiener filter coefficients
w_n = [0.8360 -0.7853];%v2 variance 0.1 case
%w_n = [0.9807 -0.9271]; %v2 variance o.01 case

d_n_w = my_filter(u_n,w_n); %wiener filter utput using my_filter function from Lab1

err_2 = (d_n - d_n_w).^2;%squared error

MSE = mean(err_2)  % Mean Squared Error

% plots of required signals
figure;
subplot(311);plot(d_n);title('Desired Signal d(n)');grid on; 
subplot(312);plot(x_n);title('Channel Output x(n)');grid on; 
subplot(313);plot(u_n);title('Distorted Signal u(n)');grid on; 
xlabel('time n');

figure
subplot(211);plot(d_n);title('Desired Signal d(n)');grid on; 
subplot(212);plot(d_n_w);title('Wiener Filter Output d\^(n) ');grid on; 
grid on;xlabel('time n');