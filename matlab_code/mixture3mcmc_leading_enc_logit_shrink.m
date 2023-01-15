function [alphaHmc,alphamc,bmc,Qmc,Qinvmc,Qinvdet,postQS,postQnu,Smc,LIMSmc,sepsmc,lambmc,postseps,postlamb,gammc,acc,etamc,rhomc,etaMSmc,etapostmc,etapostMSmc,anmc,ancholmc,dlo_miss] = ... 
    mixture3mcmc_leading_enc_logit_shrink(y,ZG,Zlogit,bas_cat,mr,sr,ZMS,miss,dlo,dlo_out,indexb,indexn,index0,isp_lag_dlo,lag_dir,maxlag,IMS,LIMSind,pool,S,ind_cont,ind_lead,Istar,rho,...
    pralm,pralinf,pralmH,pralinfH,prQnu,prQS,Q,prseps,seps,prlamb,lamb,unit_spec_var,K,eta,etaMS,etaMSind,alpha,alphaH,gam,var_logit,prgamm,prgaminf,perm,restMS,e0,e0MS,e0MSind,frac,M,it0,rth,itpost,post,ieff_restr);

%[alphamc,bmc,Qmc,Qinvmc,Qinvdet,postQS,postQnu,Smc,IMSmc,sepsmc,postseps,etamc,etaMSmc,etapostmc,etapostMSmc,anmc,ancholmc,dlo_miss] = mixture3mcmc_neu(y,ZG,ZMS,miss,dlo,dlo_out,indexb,indexn,index0,isp_lag_dlo,IMS,pralm,pralinf,prQnu,prQS,Q,prseps,seps,K,eta,etaMS,alpha,perm,e0,e0MS,frac,M,it0,it0cov,it0b)
% latent class mixed effect modell

% die Matrix ZG besteht aus einem variablen und einem fixen Teil
% ZG(:,1:dd,:) .. variabel .. X2
% ZG(:,dd+1:end,:) .. fix .. X1
% die Matrix ZMS enthaelt spalten von ZG, welche ueber die Zeit einen gruppenspezifischen, wechselnden (switchenden) Einfluss haben
% ZMS(:,1:dMS,:) 


% Modell: y = X1*alpha + X2*D(1)*betak(1) + .. + X2*D(K)*betak(K) + X3*D(K)betaR(K)*(I_t^K,-1) + seps  fuer alle i (Index weggelassen)
% mit K .. Anzahl an Klassen
% D(k) .. Dummy fuer Klassenzugehörigkeit
% seps - N(0,sigmaeps^2*I)
% b(i) - N(0,Q(i)); Q(i) = Q(:,:,S(i)); S(i).. gibt Gruppenzug an, 

% input: 
%    - observations and design matrices:

%        y (T times N)  observations
%           T .. number of repetions
%           N .. number of objects

%        ZG (T times dd+1:end times N) design matrix for the fixed effects
%           d-dd .. dimension of  parameter alpha of fixed effects

%        ZG (T times 1:dd times N) design matrix for random effects
%           dd .. dimension of  parameter betak of random effects

%        ZMS (T times 1:dMS times N) design matrix for switching effects
%           dMS .. dimension of  parameter betaR of switching effects

%        IMS (K times T) matrix of group-specific state variables, wenn size(IMS,1)==1, dann gleiches IMS fuer alle Gruppen (mit gruppenspezifischen switching effects)
%           

%    - prior information:

%        pralm  (d+dMS times 1)/(d+dd*(k-1)+dMS*K) ... normal prior distribution of fixed effect:  mean 
%        pralinf (d+dMS times d+dMS)/(d+dd*(k-1)+dMS*K times d+dd*(k-1)+dMS*K) ... normal prior distribution of fixed effect: information matrix 
%                                                                                                                                                                       (inverse of covariance)
%                        bei dim d noch aufblasen in Programm entspechend auskommentieren

%
%        prQ  .. inverted Wishart prior distribution of Q, equal for all groups
%             prQnu ... degrees of freedom (scalar)
%             prQS  ... scale matrix (dd times dd)
%        prseps (1 times 2).. inverted gamma prior distribution of seps
%             1,1 ... degrees of freedom
%             1,2 ... scale parameter

%
%        e0 (1 times K) .. Dirichlet prior for eta
%        e0MS  Dirichlet prior for etaMS (2 times 2 times K) if group-specific transition probabilities,
%                                        (2 times 2) if group-independent transition probabilities

%    - Size and starting value for MCMC:

%        M ... size of MCMC-sample
%        it0 ... burn-in 
%        it0b ... M-it0b values are saved for Smc and IMSmc
%        it0cov ... M-it0cov values are saved for ancm and ancholmc
%        sp=M-it0;
%        spcov=fix(M-it0cov);
%        spb=min(100,M-it0b);


%        Starting value:

%          seps (scalar) variance of the observation equation
%          Q (dd times dd times K) symmetric covariance matrix of the random effects
%          alpha (d+dd*(K-1)+dMS*K) fixed effects and group specific effects
%          eta (K) weights
%          etaMS (2 times 2 times K) if group-specific transition probabilities,
%                (2 times 2) if group-independent transition probabilities
%          IMS (K times T) group-specific state variable

%    - Model:
%          perm ... to define the restriction for the permutation sampler
%          K ... to define the number of groups

% output:
%        MCMC-sample:
%          alphamc (sp times d+dd*(k-1)+dMS*K)  group specific, fixed and group-specific switching effects
%          sepsmc (sp times 1) variance of the observation equation
%          etamc (sp times K) weigths
%          etaMSmc (size(e0MS) times sp) transition probabilities
%          Smc (spb times N) group indicator
%          IMSmc (K times N times spb) state variable

%        Parameters of conditional densities used for estimation (if post==1): 
%          postQS(spcov times dd*(dd+1)/2 times K) scale matrix of inverted  Wishart posterior of Q
%          postseps(spcov)  scale parameter of inverted gamma posterior densities of seps
%          etapostmc(spcov times K)
%          anmc(spcov times d+dd*(k-1)+dMS*K)
%          ancolmc(spcov times d+dd*(k-1)+dMS*K times d+dd*(k-1)+dMS*K)


acc=0;
d = size(ZG,2); % davon sind die ersten dd Spalten variabel
dMS=size(ZMS,2); %zahl der spalten, die in der Zeit switchen
dd = size(Q,1)-dMS;
gamd=size(Zlogit,2); % dimension von gamd
S_logit=any(var_logit);
DStruc=size(Istar,1)>1;
enc_cols= [find(rho==1) find(rho==2)];

N = size(y,2);
T = size(y,1);
lag_dlo=size(isp_lag_dlo,2);
% nobs_nu=(lag_dir-lag_dlo)*(lag_dir>lag_dlo);
nobs_nu=maxlag-lag_dlo;
dlo_miss=zeros(size(dlo));

% [r1,c1] = find(Istar==1);

dlo_m=median(dlo);dlo_s=iqr(dlo);      
nst=max(sum(etaMS>0,1));
states = 1:nst;
groups = [1:K]';
etapostMS=zeros(size(etaMS));
etapostMSind=zeros(nst,nst,N*(~pool)+pool);
sp=length(it0+1:rth:M);
sps=[it0+1:rth:M];
jj=0;

IMS_enc=zeros(1,size(IMS,2));
if K>1
    for j=1:nst;
        for i=1:nst;
            IMS_enc= IMS_enc+ and((IMS(1,:)+1)==j,(IMS(2,:)+1)==i)*((j-1)*j+i);
        end
    end
end

Z=zeros(T,dd*(K-1)+d+dMS*K,N);      % aufgeblasene Matrix aus ZG: erste dd*K Spalten fuer die K Gruppen,letzte (d-K*dd) Spalten fuer die fixen Effekte
alphamc=zeros(sp,dd*(K-1)+d+dMS*K);  % dimension entprechend zu Z passend
alphaHmc=zeros(N*(~pool),d+dMS,sp);
sepsmc=zeros(sp,1);
lambmc=lamb;
if unit_spec_var
    lambmc=zeros(sp,N);
end
gammc=zeros(sp,gamd,max(max(S)+any(S==0),K));
etamc=zeros(sp,max(max(S)+any(S==0),K),N*S_logit+~S_logit);

etaMSmc=[]; etapostMSmc=[]; rhomc=[];
if and(dMS>0,K>1) 
    etaMSmc=zeros(4,4,size(etaMS,3),sp);
    rhomc=zeros(sp*DStruc,size(rho,2));

    if post
        etapostMSmc=zeros(4,4,size(etaMS,3),itpost);
%         likMSmc=zeros(nst,K,size(IMS,2),spcov);
    end
elseif and(dMS>0,K==1) 
    etaMSmc=zeros(2,2,1,sp);
    rhomc=[];

    if post
        etapostMSmc=zeros(2,2,1,itpost);
%         likMSmc=zeros(nst,K,size(IMS,2),spcov);
    end
end
if dMS>0
    LIMSmc=(zeros(size(IMS,1),size(IMS,2),sp)==1);
%     LIMSind=(zeros(N,size(IMS,2))==0);
    
    index_col=[1:size(ZG,2)];
    index_perm2=zeros(dMS,1);
    for iMS=1:dMS
        incol=index_col(all(all((ZG-ZMS(:,iMS*ones(size(ZG,2),1),:))==0,3)==1,1));
        %%incol=1;'TEST prior in mixture3mcmc_neu'
        if isempty(incol) 
            ['warning: column ' num2str(iMS) ' of ZMS not contained in ZG']
        else   
            index_perm2(iMS)=incol;
        end  
    end   
    % match group specific and fixed parameters
    iiMIX=index_perm2(index_perm2<=dd);
    iiFIX=index_perm2(index_perm2>dd)+(K-1)*dd;
    iiNS= find(1-any([kron(ones(1,length(iiMIX)),[1:dd]')-iiMIX(:,ones(1,dd))']==0,2));

    %search_colMS; % bestimmt index_perm2, gibt an in welchen spalten von ZG sich die switching variablen befinden;
end

shilf=1;
it0mmc=0;
postseps=[];
postlamb=[];
etapostmc=[];
anmc=[];
ancholmc=[];

if post
   it0mmc=simuni(M-it0,itpost)+it0;
    
    postseps=zeros(itpost,2);
    if unit_spec_var
        postlamb=zeros(itpost,2,N);
    end
    etapostmc=zeros(itpost,K+(K>2));
    anmc=zeros(itpost,dd*(K-1)+d+dMS*K);
    ancholmc=zeros(itpost,dd*(K-1)+d+dMS*K,dd*(K-1)+d+dMS*K);
end

Qinv=zeros(size(Q));
detQ=zeros(1,K);
Qq=zeros(size(qincol(Q(:,:,1))));
Qqinv=Qq;
QS=zeros(dd+dMS,dd+dMS,K);
Qnu=0;
if any(Q)
    spQ=sp; spcovQ=itpost; 
 else  
    spQ=0; spcovQ=0; 
end
% 
Qmc=zeros(spQ,(dd+dMS)*(dd+dMS+1)/2,K);             
Qinvmc=zeros(spQ,(dd+dMS)*(dd+dMS+1)/2,K);             
Qinvdet=zeros(spQ,K);             
% 
postQS=[];
postQnu=[];
if post
   postQS=zeros(spcovQ,(dd+dMS)*(dd+dMS+1)/2,K);
   postQnu=zeros(spcovQ,K);
end
% 
bmc=zeros(spQ,N,dd+dMS);
if or(var_logit==3,var_logit==4) %'MH based on a normal proposal derived from the non-linear regression model'
    sr=sqrt(pi^2/6);   mr=-0.5772;
    Zk=-Zlogit./sr;   ym=(-mr*ones(size(Zlogit,1),1))./sr;
    postgaminf=prgaminf+Zk'*Zk;
    postgamcov=inv(postgaminf);
    postgamm=Zk'*ym+prgaminf*prgamm;
    cholpostgamcov=chol( postgamcov)';
end

eps = zeros(T,N); 
%S=zeros(1,N); %defined in daten_read_b
Smc=(zeros(K,N,sp)==1);
D=zeros(K,N)==1;
Vinv = zeros(T,T,N); 
Kk = zeros(dd+dMS,T);
%B = zeros(dd+dMS,dd+dMS,N);
b = zeros(N,dd+dMS);
yb = zeros(T,N);
y_mf=zeros(T,N);
y_MS=zeros(T,N);

sp=0;  % zur Berechnung der Likelihood

indexG = kron(ones(1,K),[1:dd]); %reshape(index1(ones(1,K),:)',1,K*dd);
indexK = kron([1:K]',ones(dd,1))'; %reshape(states(:,ones(1,dd))',1,K*dd);
indexR = kron(ones(1,K),[1:dMS]);
indexM = kron([1:K]',ones(dMS,1))';
indexMIX=[1:dd*K];
indexFIX=[K*dd+1:d+(K-1)*dd];
indexMS=[(d+(K-1)*dd+1):(d+(K-1)*dd+dMS*K)];

permQ_mat=diag([ones(dd,1);-ones(dMS,1)]);
permQ_mat(1:dMS,dd+1:end)=-eye(dMS);
permal_mat=diag([ones(dd*(K-1)+d,1);-ones(dMS*K,1)]);
permal_mat(1:dd*K,end-dMS*K+1:end)=kron(eye(K),[-eye(dMS);zeros(dd-dMS,dMS)]);
alpha=alpha';

for m=1:M
    
    if any(m==[1:1000:M]) 
        m
    end
    
    ANinv = pralinf;
    a = pralinf*pralm;
    
  
    % 3. sampling the group indicator
    
    if dMS>0;
        if size(IMS,1)<K
            IMStr=IMS(ones(K,1),:)';
        else
            IMStr=IMS';
        end
    end

    if K>1
    S([ind_cont ind_lead])= [find(rho==1)*ones(1,length(ind_cont)) find(rho==2)*ones(1,length(ind_lead))];
    end
    
    for j=1:N
        if K>1
            if j~=[ind_cont ind_lead];
                if and(~S_logit,size(S,2)>2)
                    [S(j),sp, Vinv(:,:,j)] = Simstateexmix_enc(y(:,j),eta,ZG(:,:,j),ZMS(:,:,j),IMS,LIMSind(j*(~pool)==j,:),Q,pool,alpha,alphaH(j*(~pool)==j,:),indexMIX,indexFIX,indexMS,seps,lamb(j),dd,sp);
                    %adjust sp for the first time series for the marginal likelihood estimation
                elseif S_logit
                    [S(j),sp,Vinv(:,:,j)] = Simstateexmix_enc_logit(y(:,j),gam,ZG(:,:,j),ZMS(:,:,j),Zlogit(j,:),IMS,LIMSind(j*(~pool)==j,:),pool,Q,alpha,alphaH(j*(~pool)==j,:),indexMIX,indexFIX,indexMS,seps,lamb(j),dd,sp);
                end
            elseif isempty([ind_cont ind_lead]);
                if and(~S_logit,size(S,2)>2)
                    [S(j),sp, Vinv(:,:,j)] = Simstateexmix_enc(y(:,j),eta,ZG(:,:,j),ZMS(:,:,j),IMS,LIMSind(j*(~pool)==j,:),Q,pool,alpha,alphaH(j*(~pool)==j,:,:),indexMIX,indexFIX,indexMS,seps,lamb(j),dd,sp);
                    %adjust sp for the first time series for the marginal likelihood estimation
                elseif S_logit
                    [S(j),sp,Vinv(:,:,j)] = Simstateexmix_enc_logit(y(:,j),gam,ZG(:,:,j),ZMS(:,:,j),Zlogit(j,:),IMS,LIMSind(j*(~pool)==j,:),pool,Q,alpha,alphaH(j*(~pool)==j,:),indexMIX,indexFIX,indexMS,seps,lamb(j),dd,sp);
                end
            elseif or(any(j==[ind_cont ind_lead]),size(S,2)<=2)
                Vinv(:,:,j)=inv([ZG(:,1:dd,j) ZMS(:,:,j).*(IMS(S(j)*ones(dMS,1),:)'-1)]*Q(:,:,S(j))*[ZG(:,1:dd,j) ZMS(:,:,j).*(IMS(S(j)*ones(dMS,1),:)'-1)]' + (seps./lamb(j))*eye(T));
            end
%             if size(S,2)<=2
%                 Vinv(:,:,j)=inv([ZG(:,1:dd,j) ZMS(:,:,j).*(IMS(S(j)*ones(dMS,1),:)'-1)]*Q(:,:,S(j))*[ZG(:,1:dd,j) ZMS(:,:,j).*(IMS(S(j)*ones(dMS,1),:)'-1)]' + (seps./lamb(j))*eye(T));
%             end
        else
            S(j)=1;
            Vinv(:,:,j)=inv([ZG(:,1:dd,j) ZMS(:,:,j).*(IMS(S(j)*ones(dMS,1),:)'-1)]*Q(:,:,S(j))*[ZG(:,1:dd,j) ZMS(:,:,j).*(IMS(S(j)*ones(dMS,1),:)'-1)]' + (seps./lamb(j))*eye(T));
        end

        D(:,j) = S(ones(K,1),j)==groups; %group dummy

        %   1a.  Sampling the fixed and group-specific (non-random) effects from the marginal model

        % constructing the new Designmatrix Z

        Z(:,[dd*K+1:dd*(K-1)+d],j) = ZG(:,[dd+1:d],j);
        Z(:,[1:dd*K],j) = ZG(:,indexG,j);
        D2 = kron(ones(T,1),D(indexK,j)');
        %D2 = D1(ones(T,1),:);
        Z(:,[1:dd*K],j) = squeeze(Z(:,[1:dd*K],j)).*D2;
        if dMS>0
            D2=kron(ones(T,1),D(indexM,j)').*(IMStr(:,indexM)-1);
            Z(:,[dd*(K-1)+d+1:end],j)=ZMS(:,indexR,j).*D2;
        end


        ANinv = ANinv + (squeeze(Z(:,:,j)))'*squeeze(Vinv(:,:,j))*squeeze(Z(:,:,j));
        a = a + (squeeze(Z(:,:,j)))'*squeeze(Vinv(:,:,j))*squeeze(y(:,j));

    end
    
    % schur decomposition   - informationsmatrix
    
    [U, Td]=schur(ANinv);
    Td=diag(1./(abs(diag(real(Td))).^.5));
    U=real(U);
    %   ancov=U*Td.^2*U';
    an = U*Td.^2*U'*a;
    anchol = U*Td;
    alpha_old=alpha;
    alpha = an + anchol*randn(dd*(K-1)+d+dMS*K,1);
    alphaiMS=reshape(alpha(indexMS),dMS,K);

    %restricted sampling: state 1 is below-average state a_st >0
   
    if ~perm
        nrest=1:(1+(K>1))*(1-DStruc)+K*DStruc;
        if any(alphaiMS(restMS,nrest)<0)
            %               alpha = an + anchol*randn(dd*(K-1)+d+dMS*K,1);
            alpha=alpha_old;
            alphaiMS=reshape(alpha(indexMS),dMS,K);

            m
            %[int2str(m) 'alpha']
        end
    end

    
    %3.a sampling the parameters of individual processes
%     for j=find(S>K)
if ~pool
    for j=[1:N]
        ANinvH=pralinfH;
        aH=pralinfH*pralmH;

        Xtr=ZG(:,:,j);

         if (dMS>0);
             Xtr=[Xtr [ZMS(:,:,j).*[LIMSind(j*ones(dMS,1),:)'-1]]]; %X-transpose
             end

        XX = Xtr'*Vinv(:,:,j)*Xtr;
        Xy = Xtr'*Vinv(:,:,j)*y(:,j);
        
        ANinvH = ANinvH + XX;        
        aH = inv(ANinvH)*(aH + Xy);
        alphaH(j,:)=[aH + chol(inv(ANinvH))'*randn(d+dMS,1)]';
        
        alphaHMS=alphaH(j,d+1:end);
        if ~perm
        while any(alphaHMS(restMS)<0)
              alphaH(j,:)=[aH + chol(inv(ANinvH))'*randn(d+dMS,1)]';
              alphaHMS=alphaH(j,d+1:end);
              %[int2str(m) 'alphaH' int2str(j)]              
        end
    end
    end
end
    
   %  4. sampling the group probability
   if K>1
       if ~S_logit
           if length(e0)>K
               [eta,etapost] = simtransex([D;S>K],e0);
           elseif length(e0)==K
               [eta,etapost] = simtransex(D,e0);
           end
       elseif S_logit
           if var_logit==1
               [gam] = simtrans_logit(Zlogit,gam,prgamm,prgaminf,D);
           elseif var_logit==2
               [gam,mr,sr] = simtrans_logit_var2(Zlogit,gam,prgamm,prgaminf,D,mr,sr);
           elseif var_logit==3
               [gam,accm] = simtrans_logit_mh(Zlogit,gam,prgamm,prgaminf,postgamm,postgaminf,postgamcov,cholpostgamcov,D,bas_cat);
               if  m > it0,           acc=acc+accm;       end
           elseif var_logit==4
               [gam,accm] = simtrans_logit_mhrw(Zlogit,gam,prgamm,prgaminf,postgamm,postgaminf,postgamcov,cholpostgamcov,D);
               if  m > it0,           acc=acc+accm;       end

           end
       end
   end
   %alph=permute(alpha([(S(ones(dd,1),:)'-1)*dd+kron(ones(N,1),[1:dd]) kron(ones(N,1),indexFIX) reshape([indexMIX([(S(ones(dMS,1),:)'-1)*dMS+kron(ones(N,1),[1:dMS])])],N,dMS)]),[3 2 1])
   alph=[];
   alph(:,:,find(S<=K))=permute(alpha([kron(ones(dd,1),(S(find(S<=K))-1))'*dd+kron(ones(sum(S<=K),1),[1:dd]) kron(ones(sum(S<=K),1),indexFIX)]),[3 2 1]);
   if ~pool
       alph(:,:,find(S>K))= permute(alphaH(S>K,1:d),[3 2 1]);
   end
    if size(alph,1)>1
        alph=permute(alph,[3 1 2]);
    end
%     if N==1
%         alph=reshape(alph,1,size(alph,3),1);
%     end
    y_mf = squeeze(sum(ZG.*alph(ones(T,1),:,:),2));
    if dMS>0 
        alph=[];
        IMSmat(:,:,find(S<=K))=reshape(IMStr(:,S(find(S<=K))),T,1,sum(S<=K));

        alph(:,:,find(S<=K))=permute(alpha([indexMS([kron(ones(dMS,1),(S(find(S<=K))-1))'*dMS+kron(ones(sum(S<=K),1),[1:dMS])])]),[3 2 1]);
        if ~pool
            IMSmat(:,:,find(S>K))=permute(LIMSind(find(S>K),:),[2 3 1]);
            alph(:,:,find(S>K))=permute(alphaH(S>K,d+1:end),[3 2 1]);
        end
        if size(alph,1)>1
            alph=permute(alph,[3 1 2]);
        end

        y_MS = squeeze(sum( (ZMS.*(IMSmat(:,ones(1,dMS),:)-1)).*alph(ones(T,1),:,:),2));
        eps = y - y_mf - y_MS;
    else
        eps = y - y_mf;
    end
    
    if any(Q)

        for k=1:K,

            [U, Td]=schur(squeeze(Q(:,:,k)));
            Q2(:,:,k)=real(U)*diag(max(diag(real(Td)),0).^.5);

        end

        for j = 1:N
            if S(j)<=K
                ZQ=[ZG(:,1:dd,j) ZMS(:,:,j).*(IMS(S(j)*ones(1,dMS),:)'-1)];
                Kk = Q(:,:,S(j))*ZQ'*Vinv(:,:,j);
                C=ZQ*Q2(:,:,S(j))/(seps/lamb(j))^.5;
                b(j,:) = (Kk*eps(:,j))' + randn(1,dd+dMS)*chol(inv(C'*C+eye(dd+dMS)))*Q2(:,:,(S(j)))';
                if ieff_restr
                    bMS=b(j,d+1:end)+alpha(dd*(K-1)+d+([1:dMS]+(S(j)-1)*dMS));
                    while bMS(restMS)<0
                        b(j,:) = (Kk*eps(:,j))' + randn(1,dd+dMS)*chol(inv(C'*C+eye(dd+dMS)))*Q2(:,:,(S(j)))';
                        bMS=b(j,d+1:end)+alpha(dd*(K-1)+d+([1:dMS]+(S(j)-1)*dMS));
                    end
                end
            else
                b(j,:)=zeros(1,dd+dMS);
                ZQ=zeros(1,dd+dMS);
        end

        eps(:,j) = eps(:,j) - ZQ*b(j,:)';

    end
    %eps = eps1;

    %        2a. covariance of the random effects

    for k=1:K,

        Qnu(1,k)=prQnu + 0.5*(sum(D(k,:),2));
        QS(:,:,k)=prQS + 0.5*b((D(k,:)==1)',:)'*b((D(k,:)==1)',:); %warum nicht -beta^k?
        [Q(:,:,k) Qinv(:,:,k) detQ(1,k)] = raninvwi_neu(Qnu(1,k),QS(:,:,k));

    end

    else

        for k=1:K

            Q(:,:,k) = zeros(dd+dMS,dd+dMS);

        end
end

    
    %  2b. variances of observation equation   
    
    posts = prseps + .5*frac*[N*T  sum(sum(eps.^2,1).*lamb,2)];
    % 'TEST VARIANCE: prior=post'     posts=prseps;
    %%%%%
    seps = 1./gamrnd(posts(1,1),1./posts(1,2));
    
    % sample unit-specific lambda's
    if unit_spec_var
        postl = prlamb(ones(N,1),:)./2 + 0.5*[T*ones(N,1) (sum(eps.^2,1)./seps)'];
        lamb = gamrnd(postl(:,1),1./postl(:,2))';
    end
    
    
    % Sample MS-Indicator, transition probabilities 
    
    if dMS>0
        Dstat1=permute(D,[3 1 2]);
        ystat1 = squeeze(sum((-ZMS(:,indexR,:).*Dstat1(ones(size(ZMS,1),1),kron(1:K,ones(1,dMS)),:))...
            .*permute(alpha([indexMS],ones(size(ZMS,1),1),ones(size(ZMS,3),1)),[2 1 3]),2));
        lik = full_lik_re([y-y_mf-ystat1],[y-y_mf],S,ZG(:,1:dd,:),ZMS,Q,seps./lamb,size(IMS,1));

        %which lead-contemporaneous grouping is more probable
        %         istar = sim_istar(lik,Istar([r1(c1==S(ind_cont(1)))],:),etaMS,xia);
        %simulate initially the state indicator independently to get an
        %initial grouping of the series
        %         if m<it00
        %             IMS = simstate_ms_new(etaMS,lik)'-1;
        %         else
        if K>1;
            if DStruc
                rho = sim_istar(lik,Istar,etaMS,etaMSind);
                if any(m==[1:20:M])
                    [m rho]
                end
                enc_cols= [find(rho==1) find(rho==2)];
            end
            

            [IMS(enc_cols,:),IMS_enc]= sim_ecompIMS(etaMS,lik(:,enc_cols,:));
            [etaMS,etapostMS] = simxi_ms(IMS_enc',e0MS);

            %fill in implied group-specific transition probabilities
            etaMSind(:,:,find(rho==1))=[(etaMS(2,2)+1)/2 1-((etaMS(2,2)+1)/2); 1-((etaMS(3,3)+1)/2) (etaMS(3,3)+1)/2];
            etaMSind(:,:,find(rho==2))=[(etaMS(1,1)+1)/2 1-((etaMS(1,1)+1)/2); 1-((etaMS(4,4)+1)/2) (etaMS(4,4)+1)/2];
        end

        if ~pool
        %     Sample MS-indicator of individual series
            Dstat1=permute(S<=(K+1),[3 1 2]);
            ystat2 = squeeze(sum(ZG.*permute(alphaH(:,1:d,ones(1,T)),[3 2 1]),2));
            ystat1 = squeeze(sum((-ZMS.*Dstat1(ones(size(ZMS,1),1),ones(1,dMS),:))...
                .*permute(alphaH(:,d+1:end,ones(1,T)),[3 2 1]),2));
            %         lik = full_lik([y-y_stat2-ystat1],[y-y_stat2],cumsum(S>K).*(S>K),seps./lamb,sum(S>K));
            lik = full_lik([y-ystat2-ystat1],[y-ystat2],[1:N],seps./lamb,N);
            LIMSind = (simstate_ms_new(etaMSind,lik)'-1)==1;
            [etaMSind,etapostMSind] = simxi_ms(LIMSind'+1,e0MSind);
        elseif pool
            if K>1
                IMS(find(rho==0),:) = simstate_ms_new(etaMSind(:,:,find(rho==0)),lik(:,find(rho==0),:))'-1;
                [etaMSind(:,:,find(rho==0)),etapostMSind(:,:,find(rho==0))] = simxi_ms(IMS(find(rho==0),:)'+1,e0MSind(:,:,find(rho==0)));
            elseif K==1
                IMS(1,:) = simstate_ms_new(etaMSind,lik(:,1,:))'-1;
                [etaMSind,etapostMS] = simxi_ms(IMS(1,:)+1,e0MSind);
                etaMS=etaMSind;
            end

        end
            
    end

    if miss==1
        % sample missing values
        Tall=size(dlo,1)-nobs_nu;
        dlo_use=dlo(1+nobs_nu:end,:);
        for i=1:N
            ib=indexb(i);
            imiss=dlo_out(1+nobs_nu:end,ib)'==1;
            if any(imiss)
                if S(i)<=K
                betai=[alpha((S(i)-1)*dd+[1:dd]);alpha(indexFIX);zeros(dMS,1)];
                elseif S(i)>K
                    betai=[alphaH(i,1:d)';zeros(dMS,1)];
                end
                betai=betai(:,ones(1,T));
                if dMS>0
                    if S(i)<=K
                        betai(d+1:d+dMS,:)=betai(d+1:d+dMS,:)+alpha(indexMS([(S(i)-1)*dMS+[1:dMS]]),ones(1,T)).*(permute([IMS(S(i),:,ones(dMS,1))],[3,2,1])-1);
                    elseif S(i)>K
                        betai(d+1:d+dMS,:)=betai(d+1:end,:)+alphaH(i*ones(1,T),d+1:end)'.*(permute([LIMSind(i,:,ones(dMS,1))],[3,2,1])-1);
                    end
                end

                if lag_dlo>0
                    isp_nolag=[1:(isp_lag_dlo(1)-1) (isp_lag_dlo(end)+1):size(ZG,2)];
                    %construct Bi
                    Bi=reshape([betai(isp_lag_dlo(end:-1:1),:)' -ones(T,1) zeros(T)]',T*(Tall+1),1);
                else
                    isp_nolag=[1:size(ZG,2)];
                    Bi=reshape([-ones(T,1) zeros(T)]',T*(Tall+1),1); %richtig?
                end   
                Bi=reshape(Bi(1:end-T),Tall,T)';
                Vst=inv(chol(inv(Vinv(:,:,i))))';
                Bi=Vst*Bi;
                ZGst=Vst*ZG(:,:,i);
                ZMSst=Vst*ZMS(:,:,i);
                izero=all(Bi(:,imiss)==0,2);
                Bitilde=Bi(~izero,imiss);
                dc=Bi(~izero,~imiss)*dlo_use(~imiss,ib)+sum(squeeze(ZG(~izero,isp_nolag,i)).*betai(isp_nolag,~izero)',2);
                if dMS>0; %THIS WAS NOT IN THE ORIGINAL PROGRAMS!
                    dc=dc+sum(squeeze(ZMS(~izero,:,i)).*betai(d+1:end,~izero)',2);
                end
                C0=(1.34/5)^2*diag((1./dlo_s(ib*ones(1,size(Bitilde,2)))).^2);
                Sy=inv(Bitilde'*Bitilde+C0);                                           
                my=Sy*(-Bitilde'*dc+C0*dlo_m(ib*ones(1,size(Bitilde,2)))');
                
                ytilde=my+chol(Sy)'*randn(size(my,1),1);
                dlo_use(imiss,ib)=ytilde;
            end   
            
        end   
        dlo(1+nobs_nu:end,:)=dlo_use;
        % impute missing values
        y=dlo(indexn,indexb);
        if lag_dlo>0
            for i1=1:lag_dlo
                ZG(1:size(indexn,1),isp_lag_dlo(i1),1:size(indexb,2))=dlo(indexn-i1,indexb);
            end
        end   
    end   
    

% abspeichern    
    if and(perm,~DStruc)
        %si = randperm(K); 
        si=[1:K];
        siMS = randperm(1+(dMS>0));
        %siMS=[1:2];  %keine permutation 
    elseif  and(~perm,DStruc)
        si =randperm(K);
        siMS=[1:2];
    else
        si=[1:K];
        siMS=[1:2];
    end
    
    %CHANGE THE FOLLOWING IF FIXED EFFECTS ARE ALSO SWITCHING
    alphai=reshape(alpha(indexMIX),dd,K);
    alphaiMS=reshape(alpha(indexMS),dMS,K);
    
    if and(dMS>0,siMS(1)==2)
        if K>1
        siMSenc=[4:-1:1];
        else
            siMSenc=siMS;
        end
        
        %CHANGE THE FOLLOWING IF FIXED EFFECTS ARE ALSO SWITCHING
        alphai(iiMIX,:,1)=alphai(iiMIX,:)-alphaiMS;
        %alpha(iiFIX)=alpha(iiFIX)-alpha(indexMS(???));  %CHOOSE WHICH STATE VARIABLE IS RELEVANT FOR FIXED EFFECTS
        alphaiMS=-alphaiMS;
        IMS = 1-IMS;
        etaMS=etaMS(siMSenc,siMSenc,:);
%         xia=xia(siMS,siMS,:);
        alpha(indexMIX,1) = reshape(alphai,dd*K,1);
        alpha(indexMS,1) = reshape(alphaiMS,dMS*K,1);
        %        lik=lik(siMS,:,:);
        if ~pool
            alphaH(:,iiMIX)=alphaH(:,iiMIX)-alphaH(:,d+1:end);
            alphaH(:,d+1:end)=-alphaH(:,d+1:end);
            LIMSind = 1-LIMSind;
        end
        etaMSind=etaMSind(siMS,siMS,:);
        if any(Q)
            b(:,iiMIX)=b(:,iiMIX)-b(:,dd+1:end);
            b(:,dd+1:end)=-b(:,dd+1:end);
            for k=1:K
                Q(:,:,k)=permQ_mat*Q(:,:,k)*permQ_mat';
                Qinv(:,:,k)=inv(permQ_mat')*Qinv(:,:,k)*inv(permQ_mat);
            end
        end
        pralm=permal_mat*pralm;
    end

    if any(si~=[1:K])
        [si1,ssort]=sort(si);
        S = ssort(S);
        %CHANGE THE FOLLOWING IF FIXED EFFECTS ARE ALSO SWITCHING
        alpha(indexMIX,1) = reshape(alphai(:,si),dd*K,1);
        if dMS>0 alpha(indexMS,1) = reshape(alphaiMS(:,si),dMS*K,1); end
        
        lik=lik(:,si,:);
        if any(Q)
        Q = Q(:,:,si);
        Qinv=Qinv(:,:,si);
        detQ=detQ(1,si);
        end
        if size(IMS,1)>1;
            IMS = IMS(si,:);
        end
        if and(size(etaMSind,3)>1,pool);
            etaMSind=etaMSind(:,:,si);
        end
        if S_logit
            gam=gam(:,si);
        elseif ~S_logit
            eta = eta(si);
        end
        if DStruc
            rho=rho(si);
        end
    end

    
    % 6. abspeichern
    if any(m==sps) 
        jj=jj+1;
        alphamc(jj,:)=alpha';
        if ~pool
        alphaHmc(:,:,jj)=alphaH;
        end
        if any(Q)
            bmc(jj,:,:)=permute(b,[3 1 2]);
            for k=1:K,
                Qq(:,k)=qincol(Q(:,:,k));
                Qqinv(:,k)=qincol(Qinv(:,:,k));
            end
            Qinvdet(jj,:)=detQ;
            Qmc(jj,:,:) = Qq;
            Qinvmc(jj,:,:) = Qqinv;
        end
        sepsmc(jj)=seps;
%         istarmc(m-it0,:)=istar;
        if unit_spec_var
            lambmc(jj,:)=lamb;
        end
        if S_logit
            gammc(jj,:,:)=gam;
            etag=exp(Zlogit*gam);
            etamc(jj,:,:)=permute(etag./kron(ones(1,size(etag,2)),sum(etag,2)),[3 2 1]);
        elseif ~S_logit
            etamc(jj,:)=eta;
        end
        %         xiamc(:,:,m-it0)=xia;
        if dMS>0 
            etaMSmc(:,:,:,jj)=etaMS; 
            LIMSmc(:,:,jj)=(IMS==1); 
        end
        Smc(:,:,jj)= S(ones(K,1),:)==groups(:,ones(N,1)); 
        if DStruc
            rhomc(jj,:)=rho;
        end
        dlo_miss=dlo_miss+dlo;
    end
    
    if any(it0mmc==m)
        shilf=find(it0mmc==m);
        sn=sum(it0mmc==m);
        anmc(shilf,:)=an(:,ones(1,sn))';
        ancholmc(shilf,:,:)=anchol(ones(1,sn),:,:);
        postseps(shilf,:)=posts(ones(1,sn),:);
        if unit_spec_var
            postlamb(shilf,:,:)=permute(postl(:,:,ones(1,sn)),[3 2 1]);
        end
        etapostmc(shilf,:)=etapost(ones(1,sn),:);
        if dMS>0 
            etapostMSmc(:,:,:,shilf)=etapostMS(:,:,:,ones(1,sn)); 
%            likMSmc(:,:,:,shilf)=lik;
        end
        if any(Q)
            for k=1:K,
                QS1(:,k)=qincol(QS(:,:,k));
            end
            postQS(shilf,:,:)=permute(QS1(:,:,ones(1,sn)),[3 1 2]);
            postQnu(shilf,:)=Qnu(ones(1,sn),:);
        end
    end
    
    
end 

dlo_miss=dlo_miss/sp;
if or(var_logit==3,var_logit==4)
    acc=acc/(M-it0);
else
    acc=1;
end    
