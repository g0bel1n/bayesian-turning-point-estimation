function [pralm, pralinf] = run_lmc_MS_priorbeta_neu(d,dd,ill_sw,K,dMS,lag_dlo,prm,prm_end,scale,scale_auto,scale_dum,auto_switch,auto_grspec)

% scale ... scale of the prior covariance
if and(~auto_switch,~auto_grspec)
    n_dum=d-dd-lag_dlo;
    scale=scale*ones(1,dd);
    scale=kron(eye(K*(1+(dMS>0))),diag(1./scale));
elseif or(auto_switch,auto_grspec)
    n_dum=d-dd;   
    scale=[scale*ones(1,dd-lag_dlo) scale_auto*ones(1,lag_dlo)];
    scale=kron(eye(K*(1+(dMS>0))),diag(ones(1,dd)./scale));
    if dMS>0
        ind_sel= reshape(ill_sw(ones(1,K),:)'+ kron(ones(dMS,1),[0:K-1]*dd),dMS*K,1)'; %index to select the appropriate elements if 
        %not all group-specific parameters are switching
        scale=scale([ind_sel dd*K+1:2*dd*K],[ind_sel dd*K+1:2*dd*K]);
    end
end
    

pralm=[kron(ones(K,1),[prm(1);prm_end*ones(dd-1,1)]);zeros(d-dd,1);kron(ones(K,1),[prm(2);zeros(dMS-1,1)])];
if dMS>0
    A=kron(ones(2,1),kron(eye(K),eye(dd)));
    A=[A [kron(-eye(K),eye(dd));zeros(K*dd)]];
    if or(auto_switch,auto_grspec)
        A=A( [ind_sel dd*K+1:2*dd*K], [1:dd*K dd*K+ind_sel]); %select the rows and columns if not all group-specific parms are switching
    end
    A=inv(A'*A)*A'; %generalized inverse
    G_inf=inv(A*scale*A');%information matrix of group- and state-specific coeff's
elseif dMS==0
    G_inf=inv(scale);
end

%insert information on group-independent coeff's                    
pralinf=eye(dd*(K-1)+d+dMS*K);
pralinf([1:dd*K dd*(K-1)+d+1:end],1:dd*K)= G_inf(:,1:dd*K);
pralinf([1:dd*K dd*(K-1)+d+1:end],dd*(K-1)+d+1:end)=G_inf(:,dd*K+1:end);
if and(~auto_switch,~auto_grspec)
pralinf(dd*K+1:dd*K+lag_dlo,dd*K+1:dd*K+lag_dlo)=pralinf(dd*K+1:dd*K+lag_dlo,dd*K+1:dd*K+lag_dlo)*scale_auto;
pralinf(dd*K+lag_dlo+1:dd*K+lag_dlo+n_dum,dd*K+lag_dlo+1:dd*K+lag_dlo+n_dum)= ... 
    pralinf(dd*K+lag_dlo+1:dd*K+lag_dlo+n_dum,dd*K+lag_dlo+1:dd*K+lag_dlo+n_dum)*scale_dum;
elseif or(auto_switch,auto_grspec)
    pralinf(dd*K+1:dd*K+n_dum,dd*K+1:dd*K+n_dum)= ... 
    pralinf(dd*K+1:dd*K+n_dum,dd*K+1:dd*K+n_dum)*scale_dum;
end
