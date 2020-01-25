function [y_n] = my_filter(x_n,h_n)
%% % Definition of variables 
    N = length(x_n); %length of the signal
    L = length(h_n); %number of taps
    y_n = zeros(N,1);%initializing the output signal by zeros
%% % Filter output calculation  
    %implementation of the summation 
    for n = 1:N %iterate over the input signal
      for k = 1:L %iterate over the filter coefficients 
          if(n - k >= 0) % discard the negative indexes
             y_n(n) = y_n(n) + h_n(k)*x_n(n-k+1); %output signal
          end
      end
    end  
end

