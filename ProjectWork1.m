%% SSP Matlab Assignment

%% TASK 1
clc
clear
close all

%Variable initialization
A = 2;
r = 0.5;
N = 100;
var_noise = 0.5;
MC = 100000;
A_estimate_vector = zeros(MC,1);

%signal generation 
time = [0:1:N-1]';
signal = A.*r.^time;

%montecarlo loop
for x = 1:MC
    %observation generation
    white_noise = randn(N,1)*sqrt(var_noise);
    observation = signal + white_noise;
    %optimal estimation
    A_estimate = sum(observation.*r.^time)/sum(r.^(2*time));
    %store result of this iteration
    A_estimate_vector(x) = A_estimate; 
end

%Check whether estimate is unbiased
mean_A_estimate = mean(A_estimate_vector)

%CRLB vs estimate variance
var_A_estimate = var(A_estimate_vector)
CLRB_A_estimate = var_noise/sum(r.^(2*time))

%plots
[A_hist, A_x] = hist(A_estimate_vector,40);

plot(A_x, A_hist/(MC*(A_x(2)-A_x(1))),'ro')
hold
plot(A_x,normpdf(A_x,A,sqrt(CLRB_A_estimate)),'b-')
hold off

%% TASK 2
clc
clear 
close all
%Variable initialization
A = 2;
var_noise = 0.5;
MC = 100000;
rho = -1;
%C = [1 rho; rho 1] * var_noise;
A_estimate_MC = zeros(MC,1);
corr_transformation = [sqrt(var_noise) 0; sqrt(var_noise)*rho sqrt(var_noise)*sqrt(1-rho^2)];

%montecarlo loop

for x = 1:MC
        %noise generation
        white_noise = randn(2,1);%independent standard gaussian noise samples       
        noise = corr_transformation * white_noise; %noise correlation

       %signal generation
       signal = A*ones(2,1) + noise;
       %estimation
       A_estimate = mean(signal);
       %save current loop estimation
       A_estimate_MC(x) = A_estimate;  
end

% mean and variance of simulation and theory
theoretical_mean = A
simulated_mean = mean(A_estimate_MC)

CRLB_A_estimate = (var_noise*(1+rho))/2
var_A_estimate_MC = var(A_estimate_MC)

%plot    
if rho ~=-1 %plot for rho different than -1
    [A_hist, A_x] = hist(A_estimate_MC,40);
    plot(A_x, A_hist/(MC*(A_x(2)-A_x(1))),'ro')
    hold
    plot(A_x,normpdf(A_x,A,sqrt(CRLB_A_estimate)),'b-')
else
    % HISTROGAM
        [A_hist, A_x] = hist(A_estimate_MC,1);
        bar(A_x, A_hist/(MC))
        hold
    % THEORETICAL PDF
        %create x-axis vector
        xxx = (theoretical_mean - 5 : 0.1 : theoretical_mean + 5)';
        %delta
        yyy = zeros(size(xxx,1),1);
        yyy(xxx == theoretical_mean) = 1;
        stem(xxx,yyy,'r-')
end

hold off

%% TASK 3
clc
clear 
close all
%initialization of variables
MC = 100000;
N = 100;
theta = 1;
theta_estimate_ML = zeros(MC,1);
theta_estimate_mean = zeros(MC,1);

for x = 1:MC
   %data generation 
   data = rand(N,1);
   theta_estimate = max(data);
   theta_estimate_ML(x) = theta_estimate;
   theta_estimate_mean(x) = 2 * mean(data); 
end

%plot
[n,x] = hist(theta_estimate_ML,73);
[n1,x1] = hist(theta_estimate_mean,73);
pdf_theo = N*x.^(N-1);
plot(x,n/(sum(n)*(x(2)-x(1))),'r-',x,pdf_theo,'bo')
figure
plot(x1,n1/(sum(n1)*(x1(2)-x1(1))),'gx',x1,normpdf(x1,theta,sqrt(theta/(3*N))),'r-')
