starting_date= 2003;

clear *mc;
index0=[];
file='logit'; % to obtain results for PDS: 'logit'; PP: 'pp',; BDS: 'bds'; UNI: 'uni'
rest_diff=0; %flag to include confidence indicators in levels
ieff_restr=0; %flag to impose state identification restriction on unit-specific effects
post=0; %flag to save the parameters of posteriors

%
not_switch=0; %flag to restrict model to non-switching
auto_switch=0; %flag to make autoregressive parameters group- and state-specific (needed for prior definition)
%(i.e. all group-specific parameters are also switching)
auto_grspec=1; %flag to make autoregressive parameters group-specific but not state-specific (needed for prior definition and for building Z)
%(i.e. not all group-specific parameters are switching)
pool=1; %flag to pool the series in the third group
start_S=0; %flag to start with fixed grouping start_S=1;
preclas=1; %flag to indicate whether some series are pre-classified into coincident or leading group
S_logit=1;% also set bas_cat and var_logit!; flag to predict the group indictaor using a logit model: S_logit=1
DStruc=0; %flag to estimate the dynamic structure in the panel -- if DStruc=1 then use pool=1, S_logit=0 in that case;
bas_cat=[3]; % set baseline group (no dependence on covariates) for logit prior
var_logit=[3];
%var_logit=3; % zur zeit verl�sslichste variante
%var_logit=2;
%var_logit=3;
%var_logit=4;

eta_different=1; %flag for different transition probabilities of each group-specific I_t, 0=equal trans. probs. for group-specific I_t

unit_spec_var=1; %flag to model unit specific error variance e_it~N(0,sigma^2/lambda_i)
sg_change=1; %flag which multiplies series by -1 which are negatively correlated with GDP
dum_out=1; %flag to model dummies as (outlying) missing values
shrink=1;

true=0; %flag to start with true values of I_t and of parameters (if 1), to start with general values (if 0),
% to start with mean value of initial estimation (if 2);

span=19;
maxlag=4; %maximum endogenous or exogenous variables

for cal_beg=[starting_date];
    cal_end=cal_beg+span;
    for cal_end_est=2022.75;%[2000.00:0.25:2006.75]; % use the second range for out-of-sample forecasting exercise
    for K=[3];['K=' int2str(K)] %number of groups
        for lag_dlo=[2]; ['lag=' int2str(lag_dlo)] %lag endogenous variable
            lag_dir=0; %lag exogenous variable

            run_model4_endswit_logit;
            %test_run_model4;

            prm=[-0.25;0.25]; %mean of group-specific parameter if I_t==0 (recession), I_t==1 (boom) ;
            prm_end=0;  %mean for autocorrelation and other group-specific variables
            scale=0.45; % prior inf. for group-specific parms gleich wie sylvia (Identit�tsmatrix), scale kleiner--info kleiner
            scale_auto=4; % prior inf. for lagged endogenous variables
            scale_dum=0.04; % prior inf. for dummy variables

            tic
            run_lcm_MS_logit;
            toc
            if preclas
                save('-mat',[file int2str(S_logit) '_dyn' int2str(DStruc) '_shr' int2str(shrink) '_var' int2str(unit_spec_var) '_perm' int2str(perm) '_sg' int2str(sg_change) '_K' int2str(K) 'end' int2str(lag_dlo) '_ex' int2str(lag_dir) '_grspecstd_' num2str(cal_end_est,'%5.2f') '_' int2str(cal_end)] );
            elseif ~preclas
                save('-mat',[file int2str(S_logit) '_pclas' int2str(preclas) '_dyn' int2str(DStruc) '_shr' int2str(shrink) '_var' int2str(unit_spec_var) '_perm' int2str(perm) '_sg' int2str(sg_change) '_K' int2str(K) 'end' int2str(lag_dlo) '_ex' int2str(lag_dir) '_grspecstd_' num2str(cal_end_est,'%5.2f') '_' int2str(cal_end)]);
            end

            clear functions
            % ctr=ctr+1;
        end
        end
    end
end

clear all




% %SWITCHING, AUTOREGRESSIVE PARMS NOT SWITCHING (NOT GROUP-SPECIFIC)
% GDP=1; %flag to choose between GDP and IP
% IP=0;
% ctries=['AUT';'BEL';'DEU';'ESP';'FIN';'FRA';'GRC';'IRL';'ITA';'NLD';'PRT';'DNK';'GBR';'SWE';'CHE';'NOR';'CAN';'JPN';'USA'];
%
% % ctry=[1:12];   %ctry=[3:17 31:34];
% % ctry=[3:22];
% ctry=[3:18 20:22];
% % ctry=[34]
%
% clear *mc;
% index0=[];
% % gr_fix=0; %flag to fix the groups a priori NOT FULLY IMPLEMENTED
% %
% not_switch=0; %flag to restrict model to non-switching  FULLY IMPLEMENTED by Sylvia F-S, Juli 2003
% auto_switch=1; %flag to make autoregressive parameters group- and state-specific (needed for prior definition)
% %(i.e. all group-specific parameters are also switching)
% auto_grspec=0; %flag to make autoregressive parameters group-specific but not state-specific (needed for prior definition and for building Z)
% %(i.e. not all group-specific parameters are switching)
%
% eta_different=1; %flag for different transition probabilities of each group-specific I_t, 0=equal trans. probs. for group-specific I_t
% unit_spec_var=1; %flag to model unit specific error variance e_it~N(0,sigma^2/lambda_i)
% dum_out=0; %flag to model dummies as (outlying) missing values
%
% true=0; %flag to start with true values of I_t and of parameters (if 1), to start with general values (if 0)
% % cal_beg=[1980 1985 1990 1995];
% % cal_end=[1986:2003];
% span=22;
% maxlag=4; %maximum endogenous or exogenous variables
% for cal_beg=[1980];
%     cal_end=cal_beg+span;
%     for K=[5];['K=' int2str(K)]
%         for lag_dlo=[1:4]; ['lag=' int2str(lag_dlo)] %lag endogenous variable
%             lag_dir=0; %lag exogenous variable
%
%             run_model4;
%             %test_run_model4;
%
%             prm=0.0; %mean of group-specific parameter if I_t==1;
%             scale=0.45; % prior inf. for group-specific parms gleich wie sylvia (Identit�tsmatrix), scale kleiner--info kleiner
%             scale_auto=4; % prior inf. for lagged endogenous variables
%             scale_dum=0.04; % prior inf. for dummy variables
%
%             tic
%             run_lcm_MS;
%             toc
%             % save(['model_MS_end' int2str(lag_dlo) '_ex' int2str(lag_dir) '_long_' ctries(ctr,:) '_78']);
%             save(['model_K' int2str(K) 'end' int2str(lag_dlo) '_ex' int2str(lag_dir) '_GDP_autopool_' int2str(cal_beg) '_' int2str(cal_end) '_dmspan']);
%             clear functions
%             % ctr=ctr+1;
%         end
%     end
% end
%
% clear all
% pack
%
% % NOT SWITCHING
% GDP=1; %flag to choose between GDP and IP
% IP=0;
% ctries=['AUT';'BEL';'DEU';'ESP';'FIN';'FRA';'GRC';'IRL';'ITA';'NLD';'PRT';'DNK';'GBR';'SWE';'CHE';'NOR';'CAN';'JPN';'USA'];
% % ctry=[1:12];   %ctry=[3:17 31:34];
% ctry=[3:18 20:22];
% % ctry=[34]
%
% clear *mc;
% index0=[];
% % gr_fix=0; %flag to fix the groups a priori NOT FULLY IMPLEMENTED
% %
% not_switch=1; %flag to restrict model to non-switching  FULLY IMPLEMENTED by Sylvia F-S, Juli 2003
% auto_switch=0; %flag to make autoregressive parameters group- and state-specific (needed for prior definition)
% %(i.e. all group-specific parameters are also switching)
% auto_grspec=0; %flag to make autoregressive parameters group-specific but not state-specific (needed for prior definition and for building Z)
% %(i.e. not all group-specific parameters are switching)
%
% eta_different=0; %flag for different transition probabilities of each group-specific I_t, 0=equal trans. probs. for group-specific I_t
% unit_spec_var=1; %flag to model unit specific error variance e_it~N(0,sigma^2/lambda_i)
% dum_out=0; %flag to model dummies as (outlying) missing values
%
% true=0; %flag to start with true values of I_t and of parameters (if 1), to start with general values (if 0)
% % cal_beg=[1980 1985 1990 1995];
% % cal_end=[1986:2003];
% span=22;
% maxlag=4; %maximum endogenous or exogenous variables
%
% for cal_beg=[1980];
%     cal_end=cal_beg+span;
%     for K=[5];['K=' int2str(K)]
%         for lag_dlo=[1:4]; ['lag=' int2str(lag_dlo)] %lag endogenous variable
%             lag_dir=0; %lag exogenous variable
%
%             run_model4;
%             %test_run_model4;
%
%             prm=0.0; %mean of group-specific parameter if I_t==1;
%             scale=0.45; % prior inf. for group-specific parms gleich wie sylvia (Identit�tsmatrix), scale kleiner--info kleiner
%             scale_auto=4; % prior inf. for lagged endogenous variables
%             scale_dum=0.04; % prior inf. for dummy variables
%
%             tic
%             run_lcm_MS;
%             toc
%             % save(['model_MS_end' int2str(lag_dlo) '_ex' int2str(lag_dir) '_long_' ctries(ctr,:) '_78']);
%             save(['model_K' int2str(K) 'end' int2str(lag_dlo) '_ex' int2str(lag_dir) '_GDP_notswit_' int2str(cal_beg) '_' int2str(cal_end) '_dmspan']);
%             clear functions
%             % ctr=ctr+1;
%         end
%     end
% end
%
% clear all
% pack
%
%
% % NOT-SWITCHING, AUTO-REGRESSIVE PARMS GROUP-SPECIFIC
% GDP=1; %flag to choose between GDP and IP
% IP=0;
% ctries=['AUT';'BEL';'DEU';'ESP';'FIN';'FRA';'GRC';'IRL';'ITA';'NLD';'PRT';'DNK';'GBR';'SWE';'CHE';'NOR';'CAN';'JPN';'USA'];
% % ctry=[1:12];   %ctry=[3:17 31:34];
% ctry=[3:18 20:22];
% % ctry=[34]
%
% clear *mc;
% index0=[];
% % gr_fix=0; %flag to fix the groups a priori NOT FULLY IMPLEMENTED
% %
% not_switch=1; %flag to restrict model to non-switching  FULLY IMPLEMENTED by Sylvia F-S, Juli 2003
% auto_switch=1; %flag to make autoregressive parameters group- and state-specific (needed for prior definition)
% %(i.e. all group-specific parameters are also switching)
% auto_grspec=0; %flag to make autoregressive parameters group-specific but not state-specific (needed for prior definition and for building Z)
% %(i.e. not all group-specific parameters are switching)
%
%
% eta_different=0; %flag for different transition probabilities of each group-specific I_t, 0=equal trans. probs. for group-specific I_t
% unit_spec_var=1; %flag to model unit specific error variance e_it~N(0,sigma^2/lambda_i)
% dum_out=0; %flag to model dummies as (outlying) missing values
%
% true=0; %flag to start with true values of I_t and of parameters (if 1), to start with general values (if 0)
% % cal_beg=[1980 1985 1990 1995];
% % cal_end=[1986:2003];
% span=22;
% maxlag=4; %maximum endogenous or exogenous variables
%
% for cal_beg=[1980];
%     cal_end=cal_beg+span;
%     for K=[5];['K=' int2str(K)]
%         for lag_dlo=[1:4]; ['lag=' int2str(lag_dlo)] %lag endogenous variable
%             lag_dir=0; %lag exogenous variable
%
%             run_model4_endswit;
%             %test_run_model4;
%
%             prm=0.0; %mean of group-specific parameter if I_t==1;
%             scale=0.45; % prior inf. for group-specific parms gleich wie sylvia (Identit�tsmatrix), scale kleiner--info kleiner
%             scale_auto=4; % prior inf. for lagged endogenous variables
%             scale_dum=0.04; % prior inf. for dummy variables
%
%             tic
%             run_lcm_MS;
%             toc
%             % save(['model_MS_end' int2str(lag_dlo) '_ex' int2str(lag_dir) '_long_' ctries(ctr,:) '_78']);
%             save(['model_K' int2str(K) 'end' int2str(lag_dlo) '_ex' int2str(lag_dir) '_GDP_nonswitgrspec_' int2str(cal_beg) '_' int2str(cal_end) '_dmspan']);
%             clear functions
%             % ctr=ctr+1;
%         end
%     end
% end
%
% clear all
% pack
