% fuehrt zu random-permutation sampler im nachhinein eine entsprechende
% Sortierung/Permutation durch

% fuer alphamc, etamc, etapostmc, anmc, ancholmc,(postmc fehlt)

%eta_restr=0;
%group=1;  sort beta^G, (2) sort beta^R, (3) sort beta^G-beta^R, (4) sort etaMS_jj
% indexMSsort=[1 1 1 1]; %as many elements as groups
% indexsort=[2];
M=size(alphamc,1);
% mn=size(anmc,1);
nst=max(1,sum(etaMSmc(1,:,1,1)>0));
dMS=size(ZMS,2);
dd=size(Q0,2)-dMS;
MSmc=M;
indexMIX=[1:dd*K];
indexFIX=[K*dd+1:d+(K-1)*dd];
indexMS=[(d+(K-1)*dd+1):(d+(K-1)*dd+dMS*K)];
%neu anfang
ZG=Z;
search_colMS
iiMIX=index_perm2(index_perm2<=dd);
iiFIX=index_perm2(index_perm2>dd);
iiNS= find(1-any([kron(ones(1,length(iiMIX)),[1:dd]')-iiMIX(:,ones(1,dd))']==0,2));
permQ_mat=diag([ones(dd,1);-ones(dMS,1)]);
permQ_mat(1:dMS,dd+1:end)=-eye(dMS);
for i = 1:M

    %alphai = reshape(alphamc(i,indexMIX),dd,K);
    alphai=reshape(alphamc(i,indexMIX),dd,K);
    alphaiMS=reshape(alphamc(i,indexMS),dMS,K);
    if any(Q)
    bi=permute(bmc(i,:,:),[2 3 1]);
    for k=1:K;
        Q(:,:,k)=qinmatr(Qmc(i,:,k)');
        Qinv(:,:,k)=qinmatr(Qinvmc(i,:,k)');
    end
    end
    if nst>1
        etai=etaMSmc(:,:,1,i);
    end

    alpha=alphamc(i,:)';
    if dMS>0
        for k=1;
            if (alpha(indexMS((k-1)*dMS+indexMSsort(k)))<0)*(1-eta_restr)+(etai(1,1)>etai(4,4))*eta_restr
                si_eta=[4 3 2 1];
                %                 if etai(1,1)>etai(2,2)
                %                             if alphai(indexMSsort(k),k,1)<alphai(indexMSsort(k),k,1)-alphai(indexMSsort(k),k,2)
                % match group specific and fixed parameters
                %                 iiMIX=index_perm2(index_perm2<=dd);
                %iiFIX=index_perm2(index_perm2>dd);
                alphai(iiMIX,:)=alphai(iiMIX,:)-alphaiMS;
                %alpha(indexFIX(iiFIX))=alpha(indexFIX(iiFIX))-alpha(indexMS(iiFIX));
                alphaiMS=-alphaiMS;
                alphaHmc(:,iiMIX,i)=alphaHmc(:,iiMIX,i)-alphaHmc(:,d+1:end,i);
                alphaHmc(:,d+1:end,i)=-alphaHmc(:,d+1:end,i);
                if any(Q)
                    bi(:,iiMIX)=bi(:,iiMIX)-bi(:,dd+1:end);
                    bi(:,dd+1:end)=-bi(:,dd+1:end);
                    for k=1:K
                        Q(:,:,k)=permQ_mat*Q(:,:,k)*permQ_mat';
                        Qinv(:,:,k)=inv(permQ_mat')*Qinv(:,:,k)*inv(permQ_mat);
                        Qmc(i,:,k)=qincol(Q(:,:,k))';
                        Qinvmc(i,:,k)=qincol(Qinv(:,:,k))';
                    end
                end

                LIMSmc(:,:,i) = (1-LIMSmc(:,:,i))==1;
                etaMSmc(:,:,1,i)=etaMSmc(si_eta,si_eta,1,i);
%                 if pool
%                     LIMSmc(3,:,i) = (1-LIMSmc(3,:,i))==1;
%                 end
            end
        end
        alphamc(i,indexMIX)=reshape(alphai,dd*K,1)';
        alphamc(i,indexMS)=reshape(alphaiMS,dMS*K,1)';
        if any(Q)
            bmc(i,:,:)=permute(bi,[3 1 2]);
        end
    end

    % Group specific parameters
        alpha1 = alphai;
        alpha1MS=alphaiMS;
        
        if group
            
            if group==1
                [alph,si] = sort(alpha1(indexsort,:)');
            elseif group==5
                [rho,si] = sort(rhomc(i,:)');
            end

            Smc(:,:,i) = Smc(si,:,i);
            %CHANGE THE FOLLOWING IF FIXED EFFECTS ARE ALSO SWITCHING
            
            alphamc(i,indexMIX) = reshape(alphai(:,si),dd*K,1)';
            if dMS>0
                alphamc(i,indexMS) = reshape(alphaiMS(:,si),dMS*K,1)';
                LIMSmc(:,:,i) = LIMSmc(si,:,i);
            end

            if any(Q)
                Q = Q(:,:,si);
                Qinv=Qinv(:,:,si);
                
                for k=1:K,
                    Qq(:,k)=qincol(Q(:,:,k));
                    Qqinv(:,k)=qincol(Qinv(:,:,k));
                end
                Qinvdet(i,:)=Qinvdet(i,si);
                Qmc(i,:,:) = Qq;
                Qinvmc(i,:,:) = Qqinv;
            end
            
            if S_logit
                gammc(i,:,:)=gammc(i,:,si);
                
            elseif ~S_logit
                etamc(i,:) = etamc(i,si);
            end
            if DStruc
                rhomc(i,:)=rhomc(i,si);
            end
        end
end
    
iiMIX_cols=reshape(iiMIX(:,ones(K,1))+  kron(dd(ones(length(iiMIX),1)),[0:K-1]),length(iiMIX)*K,1)';
am_MIX=reshape(mean(alphamc(:,indexMIX)),dd,K);
am_FIX=mean(alphamc(:,indexFIX));
t_MIX=reshape(mean(alphamc(:,indexMIX)),dd,K)./reshape(std(alphamc(:,indexMIX)),dd,K);
t_FIX=mean(alphamc(:,indexFIX))./std(alphamc(:,indexFIX));
sd_MIX=reshape(std(alphamc(:,indexMIX)),dd,K);
sd_FIX=std(alphamc(:,indexFIX));
if any(Q)
    Qm=mean(Q,1);
    sd_Q=std(Q,[],1);
    t_Q=Qm./sd_Q;
    [ac, hm, sd, med, ki, ineff] = autocovneu(reshape(Qmc,M,size(Qmc,2)*K),10);
    ki_Q=reshape(ki',2,size(Qmc,2),K);
end

% if or(auto_switch,auto_grspec)
%     indexLAG=reshape(isp_lag_dlo(ones(K,1),:)'+  kron(dd(ones(lag_dlo,1)),[0:K-1]),lag_dlo*K,1)';
%     lag_sum=permute(sum(reshape(alphamc(:,indexLAG)',lag_dlo,K,M),1),[3 2 1]);
%     am_unc= mean(alphamc(:,indexMIX(1:dd:end))./ (1-lag_sum),1);
%     lag_sum_rec=lag_sum;
%     if auto_switch
%         lag_sum_rec=permute(sum(reshape((alphamc(:,indexLAG)-alphamc(:,indexMS(indexLAG)))',lag_dlo,K,M),1),[3 2 1]);
%     end
%     am_unc_rec= mean((alphamc(:,indexMIX(1:dd:end))-alphamc(:,indexMS(indexMIX(1:dd:end))))./ (1-lag_sum_rec),1);
%
%     [ac, hm, sd, med, ki, ineff] = autocovneu([alphamc(:,indexMIX(1:dd:end))./(1-lag_sum)...
%                     (alphamc(:,indexMIX(1:dd:end))-alphamc(:,indexMS(indexMIX(1:dd:end))))./ (1-lag_sum_rec)],10);
%     ki_unc= reshape(ki,1,nst*K,2);
% end
if nst>1

    t_MS=reshape(mean(alphamc(:,indexMS))./std(alphamc(:,indexMS)),dMS,K);
    am_MS=reshape(mean(alphamc(:,indexMS)),dMS,K);
    sd_MS=reshape(std(alphamc(:,indexMS)),dMS,K);

    am_REC=reshape([mean([reshape(alphamc(:,iiMIX_cols),M,dMS,K)-reshape(alphamc(:,indexMS),M,dMS,K)],1)],dMS,K);
    t_REC=reshape([(mean([reshape(alphamc(:,iiMIX_cols),M,dMS,K)-reshape(alphamc(:,indexMS),M,dMS,K)],1))./(std([reshape(alphamc(:,iiMIX_cols),M,dMS,K)-reshape(alphamc(:,indexMS),M,dMS,K)],[],1))],dMS,K);
    sd_REC=reshape(std([reshape(alphamc(:,iiMIX_cols),M,dMS,K)-reshape(alphamc(:,indexMS),M,dMS,K)],[],1),dMS,K);
    [ac, hm, sd, med, ki, ineff] = autocovneu([alphamc alphamc(:,iiMIX_cols)-alphamc(:,indexMS)],10);
    ki_MIX=reshape(ki(indexMIX,:),dd,K,2);
    ki_REC=reshape(ki(indexMS(end)+1:end,:),dMS,K,2);
    ki_MS=reshape(ki(indexMS,:),dMS,K,2);

    if or(auto_switch,auto_grspec)
        indexLAG=reshape(isp_lag_dlo(ones(K,1),:)'+  kron(dd(ones(lag_dlo,1)),[0:K-1]),lag_dlo*K,1)';
        lag_sum=permute(sum(reshape(alphamc(:,indexLAG)',lag_dlo,K,M),1),[3 2 1]);
        am_unc= mean(alphamc(:,indexMIX(1:dd:end))./ (1-lag_sum),1);
        if auto_switch
            lag_sum_rec=permute(sum(reshape((alphamc(:,indexLAG)-alphamc(:,indexMS(indexLAG)))',lag_dlo,K,M),1),[3 2 1]);
            am_unc_rec= mean((alphamc(:,indexMIX(1:dd:end))-alphamc(:,indexMS(indexMIX(1:dd:end))))./ (1-lag_sum_rec),1);

            [ac, hm, sd, med, ki, ineff] = autocovneu([alphamc(:,indexMIX(1:dd:end))./(1-lag_sum)...
                (alphamc(:,indexMIX(1:dd:end))-alphamc(:,indexMS(indexMIX(1:dd:end))))./ (1-lag_sum_rec)],10);
        end
        if auto_grspec
            lag_sum_rec=lag_sum;


            am_unc_rec= mean((alphamc(:,indexMIX(1:dd:end))-alphamc(:,indexMS))./ (1-lag_sum_rec),1);

            [ac, hm, sd, med, ki, ineff] = autocovneu([alphamc(:,indexMIX(1:dd:end))./(1-lag_sum)...
                (alphamc(:,indexMIX(1:dd:end))-alphamc(:,indexMS))./ (1-lag_sum_rec)],10);
        end
        ki_unc= reshape(ki,1,nst*K,2);
    end

    if and(~auto_switch,~auto_grspec)
        indexLAG=isp_lag_dlo;
        lag_sum=permute(sum(alphamc(:,indexLAG,ones(1,K)),2),[1 3 2]);
        am_unc= mean(alphamc(:,indexMIX(1:dd:end))./ (1-lag_sum),1);
        am_unc_rec= mean((alphamc(:,indexMIX(1:dd:end))-alphamc(:,indexMS))./ (1-lag_sum),1);

        [ac, hm, sd, med, ki, ineff] = autocovneu([alphamc(:,indexMIX(1:dd:end))./(1-lag_sum)...
            (alphamc(:,indexMIX(1:dd:end))-alphamc(:,indexMS))./ (1-lag_sum)],10);
        ki_unc= reshape(ki,1,nst*K,2);
    end

end

if nst==1
    [ac, hm, sd, med, ki, ineff] = autocovneu([alphamc],10);
    ki_MIX=reshape(ki(indexMIX,:),dd,K,2);

    if or(auto_switch,auto_grspec)
        indexLAG=reshape(isp_lag_dlo(ones(K,1),:)'+  kron(dd(ones(lag_dlo,1)),[0:K-1]),lag_dlo*K,1)';
        lag_sum=permute(sum(reshape(alphamc(:,indexLAG)',lag_dlo,K,M),1),[3 2 1]);
        am_unc= mean(alphamc(:,indexMIX(1:dd:end))./ (1-lag_sum),1);

        [ac, hm, sd, med, ki, ineff] = autocovneu([alphamc(:,indexMIX(1:dd:end))./(1-lag_sum)],10);
        ki_unc= reshape(ki,1,nst*K,2);
    end

    if and(~auto_switch,~auto_grspec)
        indexLAG=isp_lag_dlo;
        lag_sum=permute(sum(alphamc(:,indexLAG,ones(1,K)),2),[1 3 2]);
        am_unc= mean(alphamc(:,indexMIX(1:dd:end))./ (1-lag_sum),1);

        [ac, hm, sd, med, ki, ineff] = autocovneu([alphamc(:,indexMIX(1:dd:end))./(1-lag_sum)],10);
        ki_unc= reshape(ki,1,nst*K,2);
    end

end


if nst>1;
    etadiag=[];
    etadiag_enc=[];
    for j=1:size(etaMSmc,1)
        etadiag=[etadiag permute(etaMSmc(j,j,:,:),[4 3 1 2])];
    end
    [ac, hm, sd, med, ki2, ineff] = autocovneu(etadiag,10);
    ki_eta=ki2';
%     ki_eta=reshape(ki2(1:size(etaMSmc,1),:),1,size(etaMSmc,1));
%     for j=2:nst
%         ki_eta=[ki_eta; reshape(ki2((j-1)*size(etaMSmc,3)+1:j*size(etaMSmc,3),:),1,size(etaMSmc,3),2)];
%     end

end
% am_MIXSUM=squeeze(mean(sum(reshape(alphamc(:,indexMIX),M,dd,K),2),1))
% sd_MIXSUM=squeeze(std(sum(reshape(alphamc(:,indexMIX),M,dd,K),2)))
% t_MIXSUM=squeeze(mean(sum(reshape(alphamc(:,indexMIX),M,dd,K),2),1))./squeeze(std(sum(reshape(alphamc(:,indexMIX),M,dd,K),2)))
% am_MSSUM=squeeze(mean(sum([reshape(alphamc(:,indexMIX),M,dd,K)-alphamc(:,indexMS,ones(1,K))],2),1))
% sd_MSSUM=squeeze(std(sum([reshape(alphamc(:,indexMIX),M,dd,K)-alphamc(:,indexMS,ones(1,K))],2)))
% t_MSSUM=squeeze(mean(sum([reshape(alphamc(:,indexMIX),M,dd,K)-alphamc(:,indexMS,ones(1,K))],2),1))./squeeze(std(sum([reshape(alphamc(:,indexMIX),M,dd,K)-alphamc(:,indexMS,ones(1,K))],2)))

cov_MIXMS=cov([alphamc(:,indexMIX) alphamc(:,indexMS)]);
if nst>1
    m_eta=mean(etaMSmc,4);
    s_eta=sort(etaMSmc,4);
end

