function lh = multloglikeli_MS2_re(y1,sig,nu);
% log likelihood of y1 

% y1^(0.5) is normally distributed with unit variance: y1=(y-mu)'Sigma^-1(y-mu)
% sig nxnstxK log of determinant |Sigma| that enters likelihood evaluation of group-specific state indicator k, k=1,...,K
% y1 nxnstxK  (n observations in each of nst states conditional on group membership, with K group-specific state indicators
% nu 1xK vector, number of units in group k

% Sylvia Kaufmann, 16.5.2002

nst=size(y1,2);
K=size(y1,3);
nu=reshape(nu,1,1,K);
n = size(y1,1);
flg=zeros(n,nst,K);
fln = -0.5*(nu(ones(n,1),ones(1,nst),:)*log(2*pi)+sig);%sig time dependent
flg = -0.5*y1;
lh = fln+flg;

