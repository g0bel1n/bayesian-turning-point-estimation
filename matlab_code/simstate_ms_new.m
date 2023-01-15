function [Imc, sp, p_filt,pd]= simstate_ms_new(xi,lik)

% simulate multivariate (M) indicators with L states 
% from the posterior distribution proportional to 
% likelihood*unrestricted Markov switching prior

% (all indicators are assumed to start in state j with probability 1/L)

% the procedure automatically recognizes from the third dimesnion of the transition matrix xi
% whether the transition matrices of the different processes are different or not

% Sampling method used:
% 1. multi move sampling (forward filtering - backward sampling)
% 2. if xi is the same for all processes the procedure is vectorized and 
%    samples all components I(t,:) simultaniously, 
%    otherwise there is an additional loop over m=1,...,M

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: 

% xi .. transition matrix  of  a processes with L states

% Comments:
% 1. If xi is an L x L array  (meaning that size(xi,3)=1),
%    then xi is assumed to be the same for all processes.

%    row k of xi corresponds to the conditional probability of I(t,m) given I(t-1,m)=k
%    for arbitrary m


% 2. otherwise xi is a L x L x M with M=size(xi,3)>1,
%    then xi is assumed to be different for the various processes 

%    row k of xi(:,:,m) corresponds to the conditional probability of I(t,m) given I(t-1,m)=k

% lik .. likelihood (L times M times T)
% comments:
% 1. the first index refers to time t
% 2. lik(k,m,t) is the log-likelihood of all data at time t which are relevant for process m,
%    if we assume that the process m takes the value k at time t

% Output:

% Imc ... path of the indicator sampled from the posterior (T times M);
% comments:
% 1.  Imc(t,m) refers to the value of process m at time t

% sp ... sum of the logarithm of all normalizing constants obtained during filtering (1 times M)
% comments:
% 1. sp(1,m) is identical with the log of the marginal likelihood of the conditioning parameters
%    where the m-th indicator process I is integrated out
% 2. sp might give a wrong value if you canceled constants when computing the likelihood
%    (note that even 1/sqrt(2*pi) matters in a setting with a random partion of the data
%    depending on some of the parameters of interest)

% Author: Sylvia Fruehwirth-Schnatter
% Last change: 23. September 2001

T=size(lik,3);
M=size(lik,2);   % number of indicators
L=size(xi,1); % total number of states


% 1.  forward filtering: compute the filter probabilities P(I(t,m)=k|data up to t) 
%     recursively with t running from 2 to T 
%
%     the procedure is vectorized for size(xi,3)=1 in which case the filter probabilities 
%     P(I(t,m)=k|data up to t)  are computed simultaniously for all components 
%     m = 1, ... ,M  and realisations k = 1, ..., L

p_filt=zeros(L,M,T+1);
p_pred=zeros(L,M);

p_filt(1:L,1:M,1)=1/L;
sp=zeros(1,M);


for t=1:T
   if size(xi,3)==1  
     p_pred = xi'*p_filt(:,:,t);   % p_pred is the predicitve density
  elseif size(xi,3)==M
     for m=1:M
        p_pred(:,m) = squeeze(xi(:,:,m))'*p_filt(:,m,t);  
     end     
  else
     'xi and M do not agree in simstae_ms_prior'
  end  
     maxl = max(lik(:,:,t),[],1);
   p =  p_pred.* exp(lik(:,:,t)-ones(L,1)*maxl);
   % p is the non-normalized posterior of I(t,:) 
   % (we did a numerical trick to avoid underflow)
   sp = sp + log(sum(p,1))+maxl;
   p_filt(:,:,t+1) = p./(ones(L,1)*sum(p));  
end   

% 2. backward sampling: sample I(t,:) from the smoothed probabilities P(I(t,m)=k|data up to T) 
%   backwardly with t running from T to 2  

% 2.a. sample  I(T,:) from the filter probabilities P(I(T,m)=k|data up to T)

Imc=ones(T+1,M);
indexL=[1:L]'*ones(1,M); % L times M: [ 1 ... 1; 2 ... 2; L .... L]

p = p_filt(:,:,T+1);
Imc(T+1,:) = sum(cumsum(p) < rand(L,M)) + 1;    
pd =zeros(1,M);
for m=1:M; pd(m) =log(p(Imc(T+1,m),m)); end

for t = (T):-1:1
      if size(xi,3)==1  

      p = p_filt(:,:,t) .* (xi* (ones(L,1)*Imc(t+1,:)==indexL));
   else
      for m=1:M
      p(:,m) = p_filt(:,m,t) .* (xi(:,:,m)* (ones(L,1)*Imc(t+1,m)==[1:L]'));
      end   
      end
      p = p./(ones(L,1)*sum(p));
      c=rand(L,M);
      Imc(t,:) = sum(cumsum(p) < c) + 1;          
%       for m=1:M; pd(m) =pd(m)+log(p(Imc(t,m),m)); end

end

Imc=Imc(2:end,:);

% STILL TO DO:
% richtige berechnung von sp testen