function lik = full_lik_re(stat1,stat2,S,ZG,ZMS,Q,seps,K);

%function to compute the log-likelihood of all data relevant for the group-specific state process at time t

% stat1 = y - E(y) in state 1 conditional on group membership S_i
% stat2 = y - E(y) in state 2 conditional on group membership S_i
% ZG contains group-specific random effects
% ZMS contains the switching effects
% Q is the covariance matrix of the random effects
% S     1 x N vector of group-indicator
% seps  variance of the observation equation
% K     number of groups

% output
% lik   2 x K x T  matrix, lik(j,k,t) is the log-likelihood of all data in group k at time t if the process is in state j

% Sylvia Kaufmann, 16.5.2002

groups=[1:K];
nst=2;
N=size(S,2);
T=size(stat1,1);
lik=zeros(nst,K,T);
sepsscale=zeros(T,nst,K);
res=zeros(T,nst,K);

Z1=cat(2,ZG,-ZMS);
Z2=cat(2,ZG,zeros(T,size(ZMS,2),N));

if K>1
    D = S(ones(K,1),:)==kron(groups',ones(1,N));
elseif K==1;
    D= ones(1,N)==1;
end


for t=1:T
    for k=1:K
        sig=permute(Z1(t,:,D(k,:)),[3 2 1])*Q(:,:,k)*permute(Z1(t,:,D(k,:)),[3 2 1])'+diag(seps(D(k,:)));
        res(t,1,k)=stat1(t,D(k,:))*inv(sig)*stat1(t,D(k,:))';
        sepsscale(t,1,k)=log(det(sig));
        sig=permute(Z2(t,:,D(k,:)),[3 2 1])*Q(:,:,k)*permute(Z2(t,:,D(k,:)),[3 2 1])'+diag(seps(D(k,:)));
        res(t,2,k)=stat2(t,D(k,:))*inv(sig)*stat2(t,D(k,:))';
        sepsscale(t,2,k)=log(det(sig));
    end
end

% stat1=stat1./sqrt(seps(ones(T,1),:));
% stat2=stat2./sqrt(seps(ones(T,1),:));

% if K>1
%     D = permute([S(ones(K,1),:)==kron(groups',ones(1,N))],[3 2 1]);
% elseif K==1;
%     D= ones(1,N);
% end
% res1 = sum(stat1(:,:,ones(1,K)).*(stat1(:,:,ones(1,K)).*D(ones(T,1),:,:)),2);
% res2 = sum(stat2(:,:,ones(1,K)).*(stat2(:,:,ones(1,K)).*D(ones(T,1),:,:)),2);
% nu = squeeze(sum(D,2));
% for k=1:K
%     sepsscale(1,1,k)=sum(log(seps(D(1,:,k)==1)));
% end

lh = multloglikeli_MS2_re(res,sepsscale,sum(D,2));

for j=1:K;
    lik(:,j,:)=reshape(lh(:,:,j)',2,1,T);
end
    