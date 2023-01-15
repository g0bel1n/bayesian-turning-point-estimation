function [S, sp, Vinv] = Simstateexmix_enc_logit(y,gam,ZG,ZMS,Zlogit,IMS,LIMSind,pool,Q,alpha,alphaH,indexMIX,indexFIX,indexMS,seps,lamb,dd,sp)

% S .. (1 times N) vector of group-number

K = size(IMS,1);
T = size(y,1);
d = size(ZG,2);
dMS = size(ZMS,2);
N = size(y,2);
% if size(IMS,1)<K;
%     IMS=IMS(ones(K,1),:);
% end
IMStr=IMS';
sepsi=seps./lamb;

Igr=(1-pool);


m = zeros(T,K+Igr,N);
llh = zeros(K+Igr,N);
%ZG = squeeze(ZG);
%ZMS = squeeze(ZMS);
Vinvk = zeros(T,T,N,K);
fv =zeros(N,K+Igr);
eps = zeros(T,K+Igr,N);

llh1 = -.5*T*(log(2*pi));

alpha1=alpha';
alpha2= reshape(alpha1(indexMIX),1,dd,1,K);

alpha2(1,dd+1:d,1,1:K) = alpha1(1,indexFIX,1,ones(1,K));
alpha2=alpha2(ones(T,1),:,ones(N,1),:);
if ~pool
    alpha2=cat(4,alpha2,permute(alphaH(:,1:d,ones(T,1)),[3,2,1]));
end
m1 = permute([sum(ZG(:,:,:,ones(1,K+Igr)).*alpha2,2)],[1 4 3 2]);
%append alphaH to alpha2
%append LIMSind to IMStr
if dMS>0;
    alpha3=reshape(alpha1(indexMS),1,dMS,1,K);
    alpha3=alpha3(ones(T,1),:,ones(N,1),:);
    IMStr=IMStr(:,:,ones(N,1));
    if ~pool
        alpha3=cat(4,alpha3,permute(alphaH(:,d+1:end,ones(T,1)),[3,2,1]));
        IMStr=cat(2,IMStr,permute(LIMSind,[2 3 1]));
    end
    m1= m1 + (permute([sum(ZMS(:,:,:,ones(1,K+Igr)).*alpha3,2)],[1 4 3 2])).*(IMStr-1);
end

K=K+Igr;
ytr= reshape(y,T,1,N);
eps = ytr(:,ones(1,K),:) - m1;

if any(Q)
    XX=ZG(:,1:dd);
    XXMS=XX;

    for k=1:K
        for j=1:N;
            if dMS>0
                XXMS=[XX ZMS.*(IMStr(:,k,j)-1)];
            end
            if k<=K-Igr
                Vinvk(:,:,j,k) = inv(XXMS(:,:,j)*(squeeze(Q(:,:,k)))*XXMS(:,:,j)' + sepsi(j)*eye(T));
            elseif k>K-Igr
                Vinvk(:,:,j,k) = sepsi(j)*eye(T);
            end
            llh(k,j) = llh1 +.5*(log(det(Vinvk(:,:,j,k))))- .5*eps(:,k,j)'*Vinvk(:,:,j,k)*eps(:,k,j);
        end
    end

else
    Vinv=(1./sepsi)*eye(T);
    llh = llh1 -.5*(T*log(sepsi(ones(K,1),:)) + permute(sum(eps.^2,1),[2 3 1])./sepsi(ones(K,1),:));

end

maxl = max(llh);
lh = exp(llh - maxl(ones(1,K),:));

eta=exp(Zlogit*gam);
eta=eta./kron(ones(1,K),sum(eta,2));
p = eta'.*lh;
sp = sp + sum(log(sum(p,1))+maxl,2);

p = p ./ kron(ones(K,1),sum(p,1));
fv = p';

rnd = rand(N,1);
S = (sum(cumsum(fv,2) < rnd(:,ones(K,1)),2) + 1)';      % sampling

if any(Q)
    Vinv = Vinvk(:,:,:,S);   % Vinv conditional on S
end

