function [xi,cN] = simxi_MS(Imc,c0)

% simulates the transition matrix of multivariate (M) indicators 
% following a Markov switching prior
% from the posterior given a  realisation Imc and prior parameters
% The prior of the leading indicator is adjusted to the implied prior
% derived from the encompassing state

% Sampling method used: Gibbs sampling from the full conditional density
% (the procedure is vectorized and samples all components of xi simultaniously)

% the procedure automatically recognized from the third size of the prior parameter c0
% whether the transition matrices of the different processes is the same or not


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: 

% Imc ... path of the indicator sampled from the posterior (T times M);
% comments:
% 1.  Imc(t,m) refers to the value of process m at time t

% c0 ... parameters of the prior distribution of xi

% Comments:
% 1. If c0 is an L x L array  (meaning that size(c0,3)=1),
%    then xi is assumed to be the same for all processes 
%    all conditional distributions xi(k,.) have the same Dirichlet prior 
%    with parameter D(c0(k,1),...,c0(k,L))

% 2. otherwise c0 is a L x L x M with M=size(c0,3)>1,
%    then xi is assumed to be different for the various processes 
%    all conditional distributions xi(k,..m) have a Dirichlet prior 
%    with parameter D(c0(k,1,m),...,c0(k,L,m))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output:

% xi .. defines the transition matrix xi 
% comments:
% 1. if size(c0,3)=1 xi is L times L 
% 2. if size(c0,3)=M xi is L times L times M


% cN ... parameters of the posterior distribution of xi

% Comments:
% 1.  cN has the same dimension as c0
% 2.  If xi are assumed to be different for the various processes 
%    each conditional distribution xi(k,:,m) has a Dirichlet posterior
%    with parameters D(cN(k,1,m),...,cN(k,L,m))
% 3  If xi is assumed to be the same for all processes 
%    all  conditional distribution xi(k,.) have a Dirichlet posterior
%    Beta(cN(k,1),...,cN(k,L))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjusted: Sylvia Kaufmann
% Last change: May 2005
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M=size(Imc,2);
L=size(c0,1);
T=size(Imc,1);

if size(c0,3)==1
   %  xi is assumed to be the same for all processes 
   %  posterior of xi is a beta(a,b)-distribution with following parameters:
   
   
   cN =  c0;
   for m=1:M
      cN = cN + countrns(Imc(:,m)',L);
   end
   gam = gamrnd(cN,1);
   xi = gam ./ (sum(gam,2)*ones(1,L));
   
   
elseif size(c0,3)==M
   %  xi is assumed to be different for the various processes 
   
   xi=zeros(L,L,M);
   cN=zeros(L,L,M);

   
   for m=1:M
      cN(:,:,m) = countrns(Imc(:,m)',L) + c0(:,:,m);
      gam = gamrnd(cN(:,:,m),1);
      xi(:,:,m) = gam ./ (sum(gam,2)*ones(1,L));
   end
else
   'xi and M do not agree in simxi_ms'
end  

