%%% Matrix Decomposition Labwork 

function Matrix_Decomposition(varargin)

  
  if isempty(varargin)
    
    disp('using default simulation settings and parameters...')
        
    par.MR = 4;
    par.MT = 4; 
    par.symbols = [ -1-1i,-1+1i, +1-1i,+1+1i ];   
    par.trials = 10000; 
    par.SNRdB_list = 0:5:45; 
    
    par.detector = {'ZF','ZF_QR','ZF_Chol','ZF_LU', 'ZF_LDL','MMSE','MMSE2'};      % define detector(s) to be simulated, you can add more detectors here in this way {'ZF','MMSE'}
  else
      
    disp('use custom simulation settings and parameters...')    
    par = varargin{1}; % only argument is par structure
    
  end

  
  % -- initialization  
  par.Es = mean(abs(par.symbols).^2); 
  par.Q = log2(length(par.symbols)); % number of bits per symbol
  par.bits = de2bi(0:length(par.symbols)-1,par.Q,'left-msb');
  
  % track simulation time
  time_elapsed = 0;
  
  % track decomposition time
 QR_time = 0;
 LU_time = 0;
 LDL_time = 0;
 Cholesky_time = 0;
  
  % initialize result arrays (detector x SNR)
  res.VER = zeros(length(par.detector),length(par.SNRdB_list)); % vector error rate
  res.SER = zeros(length(par.detector),length(par.SNRdB_list)); % symbol error rate
  res.BER = zeros(length(par.detector),length(par.SNRdB_list)); % bit error rate

  % generate random bit stream (antenna x bit x trial)
  bits = randi([0 1],par.MT,par.Q,par.trials);

  % trials loop
  tic
  for t=1:par.trials
  
    % generate transmit symbol
    idx = bi2de(bits(:,:,t),'left-msb')+1;
    s = par.symbols(idx).';
  
    % generate iid Gaussian channel matrix & noise vector
    n = sqrt(0.5)*(randn(par.MR,1)+1i*randn(par.MR,1));
    H = sqrt(0.5)*(randn(par.MR,par.MT)+1i*randn(par.MR,par.MT));
    
    % transmit over noiseless channel (will be used later)
    x = H*s;
  
    % SNR loop
    for k=1:length(par.SNRdB_list)
      
      % compute noise variance (average SNR per receive antenna is: SNR=MT*Es/N0)
      N0 = par.MT*par.Es*10^(-par.SNRdB_list(k)/10);
      
      % transmit data over noisy channel
      y = x+sqrt(N0)*n;
    
      % algorithm loop      
      for d=1:length(par.detector)
          
        switch (par.detector{d}) % select algorithms
          case 'ZF', % zero-forcing detection
            [idxhat,bithat] = ZF(par,H,y);
          case 'ZF_QR', % zero-forcing detection with QR decomposition
            [idxhat,bithat,QR_time] = ZF_QR(par,H,y,QR_time);
          case 'ZF_LU', % zero-forcing detection with LU decomposition
            [idxhat,bithat, LU_time] = ZF_LU(par,H,y,LU_time);
          case 'ZF_Chol', % zero-forcing detection with Cholesky
            [idxhat,bithat, Cholesky_time] = ZF_Chol(par,H,y, Cholesky_time);
          case 'ZF_LDL', % zero-forcing detection with LDL
            [idxhat,bithat, LDL_time] = ZF_LDL(par,H,y, LDL_time);
          case 'MMSE',  % MMSE detection
            [idxhat,bithat] = MMSE(par,H,y,N0);          
          case 'MMSE2', % MMSE detection with QR applied on extended channel matrix
            [idxhat,bithat] = MMSE2(par,H,y,N0); 
          otherwise,
            error('par.detector type not defined.')      
        end

        % -- compute error metrics
        err = (idx~=idxhat);
        res.VER(d,k) = res.VER(d,k) + any(err);
        res.SER(d,k) = res.SER(d,k) + sum(err)/par.MT;    
        res.BER(d,k) = res.BER(d,k) + sum(sum(bits(:,:,t)~=bithat))/(par.MT*par.Q);      
      
      end % algorithm loop
                 
    end % SNR loop    
    
    % keep track of simulation time    
    if toc>10
      time=toc;
      time_elapsed = time_elapsed + time;
      fprintf('estimated remaining simulation time: %3.0f min.\n',time_elapsed*(par.trials/t-1)/60);
      tic
    end      
  
  end % trials loop

  % normalize results
  res.VER = res.VER/par.trials;
  res.SER = res.SER/par.trials;
  res.BER = res.BER/par.trials;
  res.time_elapsed = time_elapsed;
  
  % display decomposition times
  disp(QR_time/par.trials);
  disp(LU_time/par.trials);
  disp(LDL_time/par.trials); 
  disp(Cholesky_time/par.trials);   
    
  % -- show results (generates fairly nice Matlab plot) 
  
  marker_style = {'bo-','rs--','mv-.','kp:','g*-','c>--','yx:'};
  figure(1)
  for d=1:length(par.detector)
    if d==1
      semilogy(par.SNRdB_list,res.BER(d,:),marker_style{d},'LineWidth',2)
      hold on
    else
      semilogy(par.SNRdB_list,res.BER(d,:),marker_style{d},'LineWidth',2)
    end
  end
  hold off
  grid on
  xlabel('average SNR per receive antenna [dB]','FontSize',12)
  ylabel('bit error rate (BER)','FontSize',12)
  axis([min(par.SNRdB_list) max(par.SNRdB_list) 1e-4 1])
  legend(par.detector,'FontSize',12)
  set(gca,'FontSize',12)
  
   figure(2)
  for d=1:length(par.detector)
    if d==1
      semilogy(par.SNRdB_list,res.VER(d,:),marker_style{d},'LineWidth',2)
      hold on
    else
      semilogy(par.SNRdB_list,res.VER(d,:),marker_style{d},'LineWidth',2)
    end
  end
  hold off
  grid on
  xlabel('average SNR per receive antenna [dB]','FontSize',12)
  ylabel('vector error rate (VER)','FontSize',12)
  axis([min(par.SNRdB_list) max(par.SNRdB_list) 1e-4 1])
  legend(par.detector,'FontSize',12)
  set(gca,'FontSize',12)
  
  figure(3)
  for d=1:length(par.detector)
    if d==1
      semilogy(par.SNRdB_list,res.SER(d,:),marker_style{d},'LineWidth',2)
      hold on
    else
      semilogy(par.SNRdB_list,res.SER(d,:),marker_style{d},'LineWidth',2)
    end
  end
  hold off
  grid on
  xlabel('average SNR per receive antenna [dB]','FontSize',12)
  ylabel('symbol error rate (SER)','FontSize',12)
  axis([min(par.SNRdB_list) max(par.SNRdB_list) 1e-4 1])
  legend(par.detector,'FontSize',12)
  set(gca,'FontSize',12)
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Labwork 4 %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% zero-forcing (ZF) detector
function [idxhat,bithat] = ZF(par,H,y)


% Method 1: How Zero-Forcing can be done in one line
%  xhat = inv(H'*H)*H'*y;

% Method 2: We are breaking Method 1 down.
   % Gramian 
   HH = H'*H;
   % Matched Filter
   MF = H'*y;
   % Zero-Forcing
   %xhat = inv(HH)*MF;

  [Q,R] = my_qr(H,0);
  xhat = inv(R)*Q'*y;
 
  [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
  bithat = par.bits(idxhat,:);
end




%% zero-forcing (ZF) detector with QR Decomposition
function [idxhat,bithat, QR_time] = ZF_QR(par,H,y, QR_time)

% Method 2: We are breaking Method 1 down.
   % Gramian 
   HH = H'*H;
   % Matched Filter
   MF = H'*y;
   % Zero-Forcing
   %xhat = inv(HH)*MF;
  tic;
  [Q,R] = my_qr(HH,0); % Decompose the gramian with QR decomposition
  xhat = inv(R)*Q'*MF; 
  QR_time = QR_time + toc;
  
  [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
  bithat = par.bits(idxhat,:);
end


%% zero-forcing (ZF) detector with LU Decomposition
function [idxhat,bithat,LU_time] = ZF_LU(par,H,y,LU_time)


% Decompose the gramian with LU decomposition and solve it
  % Gramian 
   HH = H'*H;
   % Matched Filter
   MF = H'*y;
   tic;
   [L, U] = my_LU(HH);
   % Zero-Forcing
   xhat = inv(U)*inv(L)*MF;
  LU_time = LU_time + toc;
  
  [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
  bithat = par.bits(idxhat,:);
end


%% zero-forcing (ZF) detector with Cholesky Decomposition
function [idxhat,bithat, Cholesky_time] = ZF_Chol(par,H,y, Cholesky_time)


% Decompose the gramian with Cholesky decomposition and solve it
 % Gramian 
   HH = H'*H;
   % Matched Filter
   MF = H'*y;
   tic;
   [L] = my_cholesky(HH);
   % Zero-Forcing
   xhat = inv(L')*inv(L)*MF;
   Cholesky_time = Cholesky_time + toc;

  [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
  bithat = par.bits(idxhat,:);
end


%% zero-forcing (ZF) detector with LDL Decomposition
function [idxhat,bithat, LDL_time] = ZF_LDL(par,H,y,LDL_time)


% Decompose the gramian with LDL decomposition and solve it
% Gramian 
   HH = H'*H;
   % Matched Filter
   MF = H'*y;
   tic;
   [L,D] = my_ldl(HH);
   % Zero-Forcing
   xhat = inv(L')*inv(D)*inv(L)*MF;
   LDL_time = LDL_time + toc;

  [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
  bithat = par.bits(idxhat,:);
end



%% MMSE detector (MMSE)
function [idxhat,bithat] = MMSE(par,H,y,N0)


[m, n] = size(H);
 % Gramian 
HH = H'*H;
% MMSE detector
xhat = inv(HH +N0*eye(n))*H'*y;

  [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
  bithat = par.bits(idxhat,:);  
end


function [idxhat,bithat] = MMSE2(par,H,y,N0)

%%%% Form augmented channel matrix and apply QR on that. Then try to solve
%%%% MMSE detection algorithm. 

[m, n] = size(H);
H=[H;sqrt(N0)*eye(n)];
[Q,R]= my_qr(H,1);

% %%%% MMSE detection algorithm.
y=[y;zeros(n,1)];
xhat = inv(R)*Q'*y;
  
  [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
  bithat = par.bits(idxhat,:);  
end


