%program to identify various estimated models

starting_date = 2003;
ending_date = 2022.75
span=19;
pool=1;
shrink=1;
rest_diff=0;
unit_spec_var=1;
perm=0;
sg_change=1;
File='logit';
DSTruc=0;
S_Logit=1;
Pclas=1;
indexMSsort=[1 1 1];indexsort=1; %which switching parameter ( in each group) is used to identify the state
                                   % which group-specific parameter
                                   % identifies the cluster
ind_enc=[1 2]; %dynamic structure identification: first contemporaneous, second leading
for cal_beg=[starting_date];
    cal_end=cal_beg+span;
    for cal_end_Est=[ending_date];
    for K=3;

        for lag_dlo=[2]; ['lag=' int2str(lag_dlo)] %lag endogenous variable
            add=0;
            group=0;
            eta_restr=0;

            lag_dir=0; %lag exogenous variable

            if Pclas
                load('-mat',[File int2str(S_Logit) '_dyn' int2str(DSTruc) '_shr' int2str(shrink) '_var' int2str(unit_spec_var) '_perm' int2str(perm) '_sg' int2str(sg_change) '_K' int2str(K) 'end' int2str(lag_dlo) '_ex' int2str(lag_dir) '_grspecstd_' num2str(cal_end_Est,'%5.2f%') '_' int2str(cal_end) ]);
            elseif ~Pclas
                load('-mat',[File int2str(S_Logit) '_pclas' int2str(Pclas) '_dyn' int2str(DSTruc) '_shr' int2str(shrink) '_var' int2str(unit_spec_var) '_perm' int2str(perm) '_sg' int2str(sg_change) '_K' int2str(K) 'end' int2str(lag_dlo) '_ex' int2str(lag_dir) '_grspecstd_' num2str(cal_end_Est,'%5.2f') '_' int2str(cal_end)]);
            end

            %             load(['logit_test_' int2str(cal_beg) '_' int2str(cal_end)])
            Q=Q0;
            eval=1;
            file=File;
            run_model4_endswit_logit;
            %             run_model_test_endswit
            permut_MS_enc;
            clear group eta_restr

            if Pclas
                save('-mat',[File int2str(S_Logit)  '_dyn' int2str(DSTruc) '_shr' int2str(shrink) '_var' int2str(unit_spec_var) '_perm' int2str(perm) '_sg' int2str(sg_change) '_K' int2str(K) 'end' int2str(lag_dlo) '_ex' int2str(lag_dir) '_grspecstd_' num2str(cal_end_est,'%5.2f') '_' int2str(cal_end) '_iden' ]);
            elseif ~Pclas
                save('-mat',[File int2str(S_logit) '_pclas' int2str(Pclas) '_dyn' int2str(DStruc) '_shr' int2str(shrink) '_var' int2str(unit_spec_var) '_perm' int2str(perm) '_sg' int2str(sg_change) '_K' int2str(K) 'end' int2str(lag_dlo) '_ex' int2str(lag_dir) '_grspecstd_' num2str(cal_end_est,'%5.2f') '_' int2str(cal_end) '_iden']);
            end
            clear functions
            %             clear ifig
        end
    end
    end
end
