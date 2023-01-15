startg=0;  startu =   0;startn =  0;    %% old

randn('state',startn);
rand('state',startu);

dd=ill; %dimension der gruppenspezifischen Effekte, die ersten dd spalten von Z
d=size(Z,2); %dimension der matrix Z
dMS=size(ZMS,2); %dimension der switchenden matrix ZMS (wieviele variablen switchen)

pral_mat=[0 1; -1 1]; %transformation of parameterization [m_1-m_2;m_1] to [m_1,m_2]
prQ_mat=diag([zeros(dMS,1);ones(dd,1)]);
prQ_mat(1:dMS,end-dMS+1:end)=eye(dMS);
prQ_mat(end-dMS+1:end,1:dMS)=-eye(dMS);

Q0 = zeros(dd+dMS,dd+dMS,K); %startwert für Q (varianz für random effects oder shrinking within groups model), wenn 0 pooling within groups
Q=eye(dd+dMS)*0.0;  %to design prior hyperparameters of Q if there is shrinkage  
if shrink
    Q=eye(dd+dMS)*0.025;
end
Q=prQ_mat*Q*prQ_mat';
N=size(y,2);

% wahl der prior
prQnu=(dd+dMS+2)*2;
prQS=Q*(prQnu-(dd+dMS))/2;
prseps=[1 1]; %

if unit_spec_var
    prlamb=[8 8]; % lambda_i
else
    prlamb=[1000  1000];
end

prm=pral_mat*prm;
[pralm, pralinf] = run_lmc_MS_priorbeta_neu(d,dd,ill_sw,K,dMS,lag_dlo,prm,prm_end,scale,scale_auto,scale_dum,auto_switch,auto_grspec);
 
pralmH=pralm(dd*(K-1)+1:dd*(K-1)+d+dMS);
pralinfH=pralinf(1:dd,1:dd);
pralinfH(end+1:end+dMS,:)=pralinf(dd*K+1:dd*K+dMS,1:dd);
pralinfH(:,end+1:end+dMS)=[zeros(dd,dMS);pralinf(dd*K+1:dd*K+dMS,dd*K+1:dd*K+dMS)];

if S_logit
    gamdim=size(Zlogit,2);
    prgamm=zeros(gamdim,1);
    prgaminf=0.05*eye(gamdim);
    e0=[];
else
    e0 = ones(1,max(max(S)+(min(S)==0),K));
    prgamm=[];
    prgaminf=[];
end


%startwerte
if true==0
    seps0=0.1; %startwert für die varianz
    lamb0=ones(1,size(y,2)); %startwert für lambda_i
    alpha_sk=[0.5*ones(dd,1);zeros(d-dd,1)];
    alpha0=pralm;
    Q0=Q(:,:,ones(1,K));
    rho0=[1:1+(K>1)];
    if K>2
        rho0=[rho0 zeros(1,K-2)];
    end

    if ~pool
        alpha0H=kron(ones(N,1),pralmH');
    elseif pool
        alpha0H=[];
    end
    
    if S_logit
        gam0=zeros(gamdim,max(max(S)+(min(S)==0),K));
        gam0(2:3,1:2)=eye(2);
        %zum testen: conditional sampling

        lamlog=Zlogit*gam0;
        lam=exp(lamlog); % predictor
        D=(rand(size(lam))<0.5)';
        U=rand(size(lam)).*(1-D')+D'.*ones(size(lam));
        ymin= -log(rand(size(lam,1),1))./sum(lam,2);
        ystar= ymin(:,ones(1,max(max(S)+(min(S)==0),K)))- log(U)./lam;
        ystarlog=log(ystar);
        [mr,sr,r]=draw_comps_mix(ystarlog+lamlog);
        eta0=[];
        % random starting value: [mr,sr,r]=draw_comps_mix(zeros(size(Zlogit,1),K));

    else
          eta0=ones(1,max(max(S)+(min(S)==0),K))/max(max(S)+(min(S)==0),K);
        gam0=0;
        mr=[];
        sr=[];
    end

    eta0MS=[1];e0MS=[1];
    eta0MSind=[1];e0MSind=[1];
    %     xia0=[1];a0=[1];
    if dMS>0

        eta0MS=[0.7 0.3 0 0;0 0.3 0 0.7;0.75 0 0.25 0;0 0 0.25 0.75];
        e0MS=[5 2 0 0;0 3 0 7;7 0 3 0; 0 0 3 7];

        eta0MSind=[0.7 0.3;0.3 0.7];e0MSind=[2 1; 1 2];
        if and(eta_different,pool)
            eta0MSind=eta0MSind(:,:,ones(1,K));
            e0MSind=e0MSind(:,:,ones(1,K));
        end
        if and(eta_different,~pool)
            eta0MSind=eta0MSind(:,:,ones(1,size(y,2)));
            e0MSind=e0MSind(:,:,ones(1,size(y,2)));
        end
    end

elseif true==2
    load(['start_values'])
    seps0=seps; %startwert für die varianz
    lamb0=lamb(indexb); %startwert für lambda_i
    alpha0=[alpha(1:d*K) alpha(d*3+1:d*3+dMS*K)];
       rho0=[1:1+(K>1)];
    if K>2
        rho0=[rho0 zeros(1,K-2)];
    end
    if shrink
        for k=1:K
            Q0(:,:,k)=Qinmatr(Q(1,:,k)');
        end
    end
    Q=Q0(:,:,1);
    
    IMS=IMS_est(1:K,1:size(IMS,2));
    
    groups=kron(ones(1,N),[1:K]');
    if K>1;
        S=groups(IS(1:K,indexb))';
    elseif K==1
        S=ones(1,N);
    end

    if ~pool
        alpha0H=kron(ones(N,1),pralmH');
    elseif pool
        alpha0H=[];
    end
    
    if S_logit
        gam0=gam;
        %zum testen: conditional sampling

        lamlog=Zlogit*gam0;
        lam=exp(lamlog); % predictor
        D=(rand(size(lam))<0.5)';
        U=rand(size(lam)).*(1-D')+D'.*ones(size(lam));
        ymin= -log(rand(size(lam,1),1))./sum(lam,2);
        ystar= ymin(:,ones(1,max(max(S)+(min(S)==0),K)))- log(U)./lam;
        ystarlog=log(ystar);
        [mr,sr,r]=draw_comps_mix(ystarlog+lamlog);
        eta0=[];
        % random starting value: [mr,sr,r]=draw_comps_mix(zeros(size(Zlogit,1),K));

    else
          eta0=ones(1,max(max(S)+(min(S)==0),K))/max(max(S)+(min(S)==0),K);
        gam0=0;
        mr=[];
        sr=[];
    end

    eta0MS=[1];e0MS=[1];
    eta0MSind=[1];e0MSind=[1];
    %     xia0=[1];a0=[1];
    if dMS>0
        eta0MS=etaMS;
        e0MS=[5 2 0 0;0 3 0 7;7 0 3 0; 0 0 3 7];

        eta0MSind=[0.7 0.3;0.3 0.7];e0MSind=[2 1; 1 2];
        if and(eta_different,pool)
            eta0MSind=eta0MSind(:,:,ones(1,K));
            e0MSind=e0MSind(:,:,ones(1,K));
        end
        if and(eta_different,~pool)
            eta0MSind=eta0MSind(:,:,ones(1,size(y,2)));
            e0MSind=e0MSind(:,:,ones(1,size(y,2)));
        end
    end

elseif true==1 %do not use, is from a previous program version
%     seps0=parms_sim(1,end);
%     lamb0=parms_sim(1,end)./parms_sim(S,end)'; %startwert für lambda_i
%     alpha0=[reshape(parms_sim(:,1:2)',dd*K,1)' parms_sim(:,end-1)']';
%    
% 
%     eta0=sum(S(ones(K,1),:)'==kron(ones(size(S,2),1),[1:K]),1)./size(S,2);
%     e0 = ones(1,K);
%     istar0=[1 2];
%     eta0MS=etaMS_sim;e0MSind=[2 1; 1 2];
%     Q0=Q(:,:,ones(1,K));
%     %eta0MS(:,:,2)=etaMS(:,:,1);
%     eta_different=size(eta0MS,3)>1; %flag for different transition probabilities of each group-specific I_t, 0=equal trans. probs. for group-specific I_t
%     xia0=xia_sim;a0=[1 2;2 1];
%     %e0MS and a0 would imply prior parameters [5 1;1 5] for the leading group's
%     %trans.probs. if imputed into trans.matrix of encompassing state
%     %[2 1 0 0;0 1 0 2;2 0 1 0;0 0 1 2];
%     %this matrix implies a prior D[4 2;2 4] for the contemporaneous
%     %group
%     % corrected to:
%     e0lead=[5 2;2 5];
%     e0cont=[4 2;2 4];
%     e0MS=cat(3,e0cont,e0lead);
%     enc_ind =  [1 1 0 0; 0 1 0 1; 1 0 1 0; 0 0 1 1]; %non-zero values of the encompassing state's trans. matrix (needed for simulation of lead durations
%     if eta_different
%         e0MS=cat(3,e0MS,e0MSind(:,:,ones(1,K-2)));
%         %             e0MS=e0MSind(:,:,ones(1,K));
%         %             eta0MS=eta0MS(:,:,ones(1,K));
%     end

end

%M=max-iteration; it0=burn-in; it0cov=M-number of retained posterior moments (equally spaced from it0 onwards);
%it0b=M-number of retained IMS and S-indicators (the last simulated values);
%L=the length of the q-sample for model likelihood estimation

if true==2
    M=9000;it0=1000;itpost=1000;L=1000;rth=1;%(not implemented)
else
    M=13000;it0=8000;itpost=1000;L=1000;rth=1;%(not implemented)
end
%  M=10;it0=2;itpost=8;L=8;
% 
% % perm=1;
% restMS=[];
perm=0;
restMS=1; %which parameter to restrict

% if isempty(ind_cont)
%     runmcmc3mixture_neu2;
% elseif ~isempty(ind_cont)
    Runmcmc3mixture_enc_logit
%end
% runmcmc3mix_ml_sepsind_MS;

clear Z y;