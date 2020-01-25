
%% Task 1
clc
close all
clear
%read data
load data_high_snr
power_data = abs(data);%for this task, the power data is not squared!!! (Rayleigh PDF)

%my_algorithm
%variable initialization
ray_std = zeros(10,1); %Rayleigh standard deviation of current iteration
thresh = zeros(10,1); %Threshold for current iteration
thresh(1) = max(power_data)+1; %first threshold is larger that all signal samples
Nsigma = 4; % Convergence parameter = number of pulse samples / number of noise samples
epsilon = 0.0001;
qqq = 0; %current iteration counter
Nsamples = size(power_data,1); %number of signal samples
Nsamples_out_qqq = Nsamples; 
Nsamples_out_qqq_1 = 0;
plot(timet,20*log10(power_data),'b-')
xlabel('Time')
ylabel('Power of data [dBW] = 10*log10(abs(data)^2)')
title('Threshold for each iteration')
legend('Data samples')
hold

while(Nsamples_out_qqq_1/Nsamples_out_qqq < 1 - epsilon)
        qqq = qqq + 1; %updates iteration counter
        ray_std(qqq+1) = sqrt(var(power_data(power_data < thresh(qqq))));%computes standard deviation of signal samples smaller than previous threshold, assuming a Rayleigh pdf
        thresh(qqq+1) = Nsigma * ray_std(qqq+1)*sqrt(2/(4-pi));%Threshold = Nsigma times Sigma of Rayleigh pdf..... (Sigma of Rayleigh)^2 = 2*ray_std^2/(4-pi)
        Nsamples_out_qqq = sum(power_data < thresh(qqq));%number of samples larger than previous threshold
        Nsamples_out_qqq_1 = sum(power_data < thresh(qqq + 1));%number of samples larger than current threshold  
        plot(timet,ones(size(timet,2),1)*20*log10(thresh(qqq + 1)),'color',rand(1,3))        
end

my_thresh = 20*log10(thresh(qqq+1))

%% Task 3
ttt = thresh(qqq+1); %Task1 algorithm

noise_samples = power_data;
signal_samples = power_data;
noise_samples(noise_samples > ttt) = 0; %samples greater than threshold are set to zero
signal_samples(signal_samples <= ttt) = 0; %samples smaller than threshold are set to zero
figure
plot(timet, 10*log10(noise_samples),'ro') %taking log of samples with zero value, removes them from plot
hold
plot(timet, 10*log10(signal_samples),'bo')%taking log of samples with zero value, removes them from plot
xlabel('Time')
ylabel('Power data (10*log10(data)^2)')
title('Classification')
legend('Noise','Signal')


%% Task 2
clear all
clc
figure()
load('data_high_snr.mat')                   %%%%%% Loading data file
 
y=10*log10(abs(data).^2);
 
power_data = abs(data).^2;
Tcme = 4.6052;
Pr = exp(-Tcme);
M = 1;% Degree of freedom of Chi Squared Distribution
ordered_data = sort(power_data);
I = ceil(0.1 * size(power_data,1));
cumm_sum = sum(ordered_data(1:I-1));

%loop for the test
for k = I : size(ordered_data,1)-1
%   Tk = Tcme/k; %this equation is not right, we need to use fisher function
    Tk=finv(1-Pr, 2*M,2*M*k)/k;
    cumm_sum = cumm_sum + ordered_data(k);
    if ordered_data(k+1) > Tk * cumm_sum
        threshold = ordered_data(k+1);
        break
    end   
end
Threshold = 10*log10(threshold) %getting threshold in dBm

plot(timet,10*log10(abs(data).^2),'b-')
xlabel('Time')
ylabel('Power of data [dBW] = 10*log10(abs(data)^2)')
title('FCME algorith threshold')
hold
plot(timet,ones(size(timet,2),1)*10*log10(threshold),'r-')        
legend('Data samples', 'Threshold')
%% Task 5
clear all
clc
load('data_high_snr.mat') %Loading data file
y=10*log10(abs(data).^2);
plot(timet,y)
[k,v]=findpeaks(y,timet,'MinPeakHeight',-75,'MinPeakDistance',0.0001);
total_diff = 0;
no_of_peaks = size(v,2);
fprintf('No of peaks=  %f.\n',no_of_peaks);
 
 %loop for the summation of gaps between peaks
for n=1:no_of_peaks-1
  total_diff = total_diff + (v(n+1)-v(n));
%   fprintf('Difference=  %f , Diff iteration = %f.\n',total_diff,(v(n+1)-v(n)));
    fprintf(' %f.\n',(v(n+1)-v(n)));
end
 
 %Average pulse duration
pulse_duration = total_diff/size(v,2);

%% Task 6
clc
clear
close all
N = 10000;
noise = randn(N,1);
subplot(1,3,1)
plot(noise)
xlabel('samples')
ylabel('noise')

%pulse characteristics generation
signal = zeros(N,1);
Nsignals = ceil(rand(1,1)*3); %number of signals, from 1 to 3
Wsignals = 100*ones(size(Nsignals,1)); %pulse width=100
Tsignals = floor(rand(Nsignals,1)*(N-100)+50); %signals position inside the sample vector
Asignals = rand(Nsignals,1)*2+1; %signals amplitude, from 1 to 3

for q = 1:Nsignals %signal generation
   W0 = max(Tsignals(q) - Wsignals(q),1);
   W1 = min(Tsignals(q) + Wsignals(q),N);
   signal(W0:W1) = Asignals(q);
end
%Plots
subplot(1,3,2)
plot(signal)
xlabel('samples')
ylabel('signal')
subplot(1,3,3)
plot(signal+noise)
xlabel('samples')
ylabel('signal+noise')

%algorithm
sn = signal+noise;
mean_100 = zeros(N-99,1);
for w = 1:N-99
    mean_100(w) = mean(sn(w:w+99)); %moving mean
end
figure,plot(mean_100)
hold
%threshold plot
plot([1:size(mean_100)],1*ones(size(mean_100,1)),'r-')
