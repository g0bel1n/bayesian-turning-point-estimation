daten_read_france
% daten_read_test;
%  K=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aufbau von y und der Matrix Z, Zuerst variabel, dann fix

%  lag_dlo=4; %lag endogenous variable
%  lag_dir=0; %lag exogenous variable

lag_y=max(lag_dir,lag_dlo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auswahl der cross-section unit

index=[1:nbank];iout=(index~=index);

% Ausreisserbehandlung: iout(i)==1 exclude time series i

%outmerge=1;  %exclude time series containing mergers
outmerge=0;
if outmerge==1; iout=any(dlo_out==1); end

miss=1; %treat mergers as missing values
%miss=0;
if miss==1
    dlo_obs=dlo; % actual observations;
    dlo_m=median(dlo);
    for i=find(any(dlo_out))
            dlo(dlo_out(:,i)==1,i) = dlo_m(i);
            %substitute median as starting value for missing observations at mergers
    end
end

dlo_max=max(dlo);dlo_min=min(dlo);
outy=0;
%outy=1;scale=5;%exclude time series with outliers
if outy==1; dlo_m=median(dlo);dlo_s=iqr(dlo);iout=or(dlo_max>dlo_m+scale/1.34*dlo_s,dlo_min<dlo_m-scale/1.34*dlo_s) ;%'max > scale*(iqr/1.34)'
end

%liq_sh=0; %flag (1) to eliminate banks with liquidity shares above 100%
%exclude_ind; %ev. exclude some units
%iout(1,ind)=1;

plot_out=0;
if plot_out==1
    iplot=0;ifig=20;nplot=12;nline=3;xmin=-50;xmax=100;
    for i=index(iout)
        iplot=iplot+1;
        if iplot==1;ifig=ifig+1;figure(ifig);end
        subplot(nplot/nline,nline,iplot);
        plot(dlo(:,i));
        %  axis([0 nobs+2  xmin xmax]);
        alpha=0.3;
        axis([0 nobs+2  dlo_min(i)/(1-alpha) dlo_max(i)*(1+alpha)]);
        legend(num2str(BLZ(i)));
        if iplot==nplot;iplot=0;end
    end
end


indexb=index(~iout);%Auswahl der units
% indexn=[(1+lag_y):nobs]';%Ber�cksichtigen von autokorrelation in y
indexn=[(1+maxlag):nobs]';%Ber�cksichtigen von autokorrelation in y oder von verz�gerten endogenen variablen


y=dlo(indexn,indexb);
cal=cal(indexn);

Z_dlolag=[];
if lag_dlo>0
    for i1=1:lag_dlo
        Z_add(1:size(indexn,1),1,1:size(indexb,2))=dlo(indexn-i1,indexb);
        Z_dlolag=cat(2,Z_dlolag,Z_add);
    end;clear i1
end

Zvar=[];

% for iv=[1 ]; % 1 .. zinssatz
%    %Zvar=cat(2,Zvar,Zbasis(indexn,iv,indexb));
% if lag_dir>0
%    for i1=1:lag_dir
%       Z_add(1:size(indexn,1),1,1:size(indexb,2))=Zbasis(indexn-i1,iv,indexb);
%       Zvar=cat(2,Zvar,Z_add);
%    end;clear i1
% end
% end; clear iv;

%Zvar=cat(2,c(indexn,:,indexb),Z_dlolag); %autoregressive parameters group specific
%Z=Zvar;
Zvar=cat(2,c(indexn,:,indexb),Z_dlolag); %group-specific autoregressive parameters
Z=cat(2,Zvar);

if lag_dlo>0
    ill=size(Zvar,2); %gibt die letzte der spalten von Z an, in denen sich die gruppenspezifischen effekte befinden
    isp_lag_dlo=[2:ill];
    %gibt an in welchen spalten von Z sich die zeitverz�gerten Variablen befinden
else
    ill=size(Zvar,2);
    isp_lag_dlo=[];
end
isp_size=0;   %gibt an in welcher spalte von Z sich die size befindet - nicht definiert

if S_logit
    Zlogit=[ones(size(indexb))' Zlogit(indexb,:)]; %read in daten_read_cont0604.m
%   Zlogit=[ones(size(indexb))' squeeze(Z_dlolag(1,1,:))];
   %Zlogit=[ones(size(indexb))'];  % test  1: should be equal to dirichlet
                                %sampling
   %Zlogit=[ones(size(indexb))' rand(size(indexb))'];  % test  2: should
   %be insignificant

else
    Zlogit=[];
end

% ZMS enth�lt jene spalten von Z, deren Parameter sich nach einem M_S_modell �ndern k�nnen

if not_switch
    ill_sw=[];
    ZMS=[];

else
    if auto_switch
        ill_sw=[1:ill]; %gibt an, welche spalten von Z switchen
        ZMS=Z(:,ill_sw,:);
    elseif auto_grspec
        ill_sw=[1]; %gibt an, welche spalten von Z switchen
        ZMS=Z(:,ill_sw,:);
    end
end

if length(index0)>0;
    IMS=probMS;
else
    if ~DStruc
        if K==1
            IMS=[mean(y,2)>mean(mean(y))]';
            Istar=[1];
        elseif K>1
            IMS=[mean(y(:,ind_cont),2)>mean(mean(y(:,ind_cont),2))]';
            IMS=[IMS ; [mean(y(:,ind_lead),2)>mean(mean(y(:,ind_lead),2))]'];
            Istar=[1 2 zeros(1,K-2)];
        end
    elseif DStruc
        IMS=[mean(y,2)>mean(mean(y))]';
        IMS=[IMS ; IMS(3:end) IMS(ones(1,2)*end)];
        Istar=make_mat_istar(K);
    end
    if ~pool
        LIMSind=[y > kron(ones(length(indexn),1),mean(y,1))]';
    elseif pool
        IMS=[IMS ; IMS(2*ones((K-2),1),:)];
        LIMSind=[];
    end

    if not_switch
        IMS=IMS*0+1; %uncomment if only one state
    end

end
%starting value for MS-indicator: IMS(t+1) = 1, if mean of industrial production across countries_t>0, 0 otherwise
% IMS=IMS(ones(K,1),:); % one state variable for each group

%defining S
S=zeros(1,length(indexb));
S([ind_cont ind_lead])=[ones(1,length(ind_cont)) 2*ones(1,length(ind_lead))];
% if gr_fix
%     group_fix
%     S=gr_ind;
% end



clear ind_prod Zbasis Zdum Zdumli Zwirt Zvar Z_dlolag Z_add  c;%startwerte
