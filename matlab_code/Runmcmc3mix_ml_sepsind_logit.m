% runmcmc3mixture_ml
fid=fopen('bf_canova.out','a');

% latent class mixed effect modell - berechnen bayes factor aus vorhandenen simulationen
% (workspace aktiv)
            
if ~S_logit
    [mllogbs, sdrel,loglikq,loglikmc,priorq,priormc,qq,qmc,mu,Qsim,Qsiminv,detQsiminv,qmall,qqall] = mlbsfullddmixpr5_sepsind(alphamc,sepsmc,etamc,Qmc,Qinvmc,Qinvdet,prQnu,prQS,postQnu,postQS,e0,etapostmc,prseps,postseps,prlambda,anmc,ancholmc,pralm,pralinf,y,Z,dd);
else
    if perm~=100
        [mllogbs, sdrel,loglikq,loglikmc,priorq,priormc,qq,qmc,mu,Qsim,Qsiminv,detQsiminv,qmall,qqall] = mlbsfullddmixpr5_sepsind_logit(alphamc,sepsmc,gammc,Zlogit,Qmc,Qinvmc,Qinvdet,prQnu,prQS,postQnu,postQS,prgamm,prgaminf,prseps,postseps,prlambda,anmc,ancholmc,pralm,pralinf,y,Z,dd);
    else
        ' BAYES FACTOR Logit auf permutation sampling ERWEITERN'
    end   
end

if perm~=100 
    mllogbs=mllogbs+gammaln(K+1);
end    
ergmiran=[mllogbs(1,1) sdrel(1) mllogbs(1,1)-5*sdrel(1) mllogbs(1,1)+5*sdrel(1)] %importance
ergriran=[mllogbs(1,2) sdrel(2) mllogbs(1,2)-5*sdrel(2) mllogbs(1,2)+5*sdrel(2)] %reciprocal imp.
ergbsran=[mllogbs(end,1) sdrel(3) mllogbs(end,1)-5*sdrel(3) mllogbs(end,1)+5*sdrel(3)]  %bridge

if and(~not_switch,~auto_switch)
    t=['model_K' int2str(K) 'end' int2str(lag_dlo)  '_canova_base          '];
elseif and(~not_switch,auto_switch)
    t=['model_K' int2str(K) 'end' int2str(lag_dlo)  '_canova_autoswit      ' ];
elseif and(not_switch,~auto_switch)
    t=['model_K' int2str(K) 'end' int2str(lag_dlo)  '_canova_nonswit       ' ];
elseif and(not_switch,auto_switch)
    t=['model_K' int2str(K) 'end' int2str(lag_dlo)  '_canova_nonswit_autoswit ' ];
end

            if ~S_logit
                t=[t  '_prlambda' num2str(prlambda) '  '];
            else
                 t=[t  '_prlambda' num2str(prlambda) '_Slogit'  num2str(var_logit) '  '];
            end   

t=[t 'imp ' num2str(ergmiran,'%6.4f %6.4f %6.4f') ' r.imp ' num2str(ergriran,'%6.4f %6.4f %6.4f') ' bsmi ' num2str(ergbsran,'%6.4f %6.4f %6.4f') '\n'];

fprintf(fid,t);
fclose(fid)

clear *post* loglik 
clear qnorm loglikq loglikmc priorq priormc qq qmc qmall qqall%ergmiran ergriran ergbsran
clear mllogbs sdrel
clear functions
