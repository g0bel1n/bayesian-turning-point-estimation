%program to plot some interesting output for various estimated models
starting_date= 2003
ending_date = 2022.75

span=19;
rest_diff=0;
pool=1;
shrink=1;
sg_change=1;
unit_spec_var=1;
perm=0;
File='logit';
DSTruc=0;
S_Logit=1;
Pclas=1;
for cal_beg=[starting_date];
    cal_end=cal_beg+span;
    for cal_end_Est=[ending_date];
    for K=[3];
        for lag_dlo=[2]; ['lag=' int2str(lag_dlo)] %lag endogenous variable
            lag_dir=0; %lag exogenous variable

            if Pclas
                load('-mat',[File int2str(S_Logit) '_dyn' int2str(DSTruc) '_shr' int2str(shrink) '_var' int2str(unit_spec_var) '_perm' int2str(perm) '_sg' int2str(sg_change) '_K' int2str(K) 'end' int2str(lag_dlo) '_ex' int2str(lag_dir) '_grspecstd_' num2str(cal_end_Est,'%5.2f') '_' int2str(cal_end) '_iden']);
            elseif ~Pclas
                load('-mat',[File int2str(S_Logit) '_pclas' int2str(Pclas) '_dyn' int2str(DSTruc) '_shr' int2str(shrink) '_var' int2str(unit_spec_var) '_perm' int2str(perm) '_sg' int2str(sg_change) '_K' int2str(K) 'end' int2str(lag_dlo) '_ex' int2str(lag_dir) '_grspecstd_' num2str(cal_end_Est,'%5.2f') '_' int2str(cal_end) '_iden']);
            end

            eval=1;
            Q=Q0;
            run_model4_endswit_logit;

            plotrd;
            plotrdetaMS_ranperm;
            if any(Q)
            plotrdQ
            end
            plotprob;
            if any(Q)
                plotbetai
            end
            if S_Logit
                plotzlogit
            end
             plotseps
            clear functions
            % ctr=ctr+1;
        end
    end
    end
end
