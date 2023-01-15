function  [gam,acc] = simtrans_logit_mh(Zlogit,gam,prgamm,prgaminf,postgamm,postgaminf,postgamcov,cholpostgamcov,D,k0);
  
K = size(gam,2);
dg=size(gam,1);
if size(D,1)<K
    D=[D;1-sum(D,1)];
end
% 1. Sample hidden propensities based on gam old

lamlog=Zlogit*gam;
lam=exp(lamlog); % predictor

 U=rand(size(lam)).*(1-D')+D'.*ones(size(lam));
 ymin= -log(rand(size(lam,1),1))./sum(lam,2);
ystar= ymin(:,ones(1,K))- log(U)./lam;
ystarlog=log(ystar);

% 3. Sample new parameter of the logit model
% index=[1:K];
%  k0=index(all(gam==0,1)==1);k0=k0(1);   % baseline category k0 
%  k0=3;
baseline=0; %baseline=0: sample also baseline category from unidentified model
baseline=1; %baseline=1: fix baseline category k0 (no permutation sampling : 1);

if baseline
    indexk=[1:k0-1 k0+1:K];
else
    indexk=index;
end

for k=indexk
    sr=sqrt(pi^2/6);
    Zk=-Zlogit./sr;
    postgamk=postgamcov*(Zk'*ystarlog(:,k)./sr+postgamm);
    gamneu=postgamk+cholpostgamcov*randn(dg,1);
    gamold= gam(:,k);

    priorratio=-0.5*((gamneu-prgamm)'*prgaminf*(gamneu-prgamm)-(gamold-prgamm)'*prgaminf*(gamold-prgamm));
     qratio=0.5*((gamneu-postgamk)'*postgaminf*(gamneu-postgamk)-(gamold-postgamk)'*postgaminf*(gamold-postgamk));
    
      lamlogneu=Zlogit*gamneu;
    lamneu=exp(lamlogneu); % predictor
     likratio=  sum(lamlogneu-ystar(:,k).*lamneu,1)- sum(lamlog(:,k)-ystar(:,k).*lam(:,k),1);
    acc=(log(rand)<(likratio+priorratio+qratio));
    if acc
        gam(:,k)=gamneu;
    else
       gam(:,k)=gamold;
   end 
end    




gam=gam-gam(:,k0*ones(1,size(gam,2)));


