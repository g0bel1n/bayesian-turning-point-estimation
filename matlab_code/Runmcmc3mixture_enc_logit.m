% runmcmc3mixture


frac=1;

[alphaHmc,alphamc,bmc,Qmc,Qinvmc,Qinvdet,postQS,postQnu,Smc,LIMSmc,sepsmc,lambmc,postseps,postlamb,gammc,acc,etamc,rhomc,etaMSmc,etapostmc,etapostMSmc,anmc,ancholmc,dlo_miss] = ...
    mixture3mcmc_leading_enc_logit_shrink(y,Z,Zlogit,bas_cat,mr,sr,ZMS,miss,dlo,dlo_out,indexb,indexn,index0,isp_lag_dlo,lag_dir,maxlag,IMS,LIMSind,pool,S,ind_cont,ind_lead,Istar,rho0,...
    pralm,pralinf,pralmH,pralinfH,prQnu,prQS,Q0,prseps,seps0,prlamb,lamb0,unit_spec_var,K,eta0,eta0MS,eta0MSind,alpha0,alpha0H,gam0,var_logit,prgamm,prgaminf,perm,restMS,e0,e0MS,e0MSind,frac,M,it0,rth,itpost,post,ieff_restr);

mixture3mcmc_leading_enc_logit_shrink(y,Z,Zlogit,bas_cat,mr,sr,ZMS,miss,dlo,dlo_out,indexb,indexn,index0,isp_lag_dlo,lag_dir,maxlag,IMS,LIMSind,pool,S,ind_cont,ind_lead,Istar,rho0,...
    pralm,pralinf,pralmH,pralinfH,prQnu,prQS,Q0,prseps,seps0,prlamb,lamb0,unit_spec_var,K,eta0,eta0MS,eta0MSind,alpha0,alpha0H,gam0,var_logit,prgamm,prgaminf,perm,restMS,e0,e0MS,e0MSind,frac,M,it0,rth,itpost,post,ieff_restr);