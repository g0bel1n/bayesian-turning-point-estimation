% loading our data file
[data ser]=xlsread('data.xlsx');

% Zlogit=load('corr_0604.asc');

ser=ser(1,2:end); % string array of data series
if ~strcmp(file,'uni');
    %PDS
    if strcmp(file,'logit');
        ser_excl={'price_products'
            }; %string of series to exclude
        ser_sel = ser(find(1-sum(strcmp(ser(ones(size(ser_excl,1),1),:),ser_excl(:,ones(1,size(ser,2)))),1))); % string array of included series
        ser_ind=find(1-sum(strcmp(ser(ones(size(ser_excl,1),1),:),ser_excl(:,ones(1,size(ser,2)))),1)); %choose which countries to include, index array
        if preclas
            ser_cont={'YER';
                'MTR';
                'PCR';
                'capital_fixe';
                'XTR';
             
                'price_food';
                'price_energy';
                'rate_production';
                
                'past_trend_workforce';
                'price_servies';
                'building_industry_past';
                'rent';
                'unemployment_rate';
                'past_service_activities';
                 'buisness_climate'}; %series falling into the contemporaneous group
            ser_lead={'expected_service_activities';
                'trend_services';
                'manufacturing_index';
                  'building_industry_expected';
                'building_industry_overall';
                }; %series falling into the leading group
            ind_cont=find(sum(strcmp(ser_sel(ones(size(ser_cont,1),1),:),ser_cont(:,ones(1,size(ser_sel,2)))),1)); %index of included series defining the contemporaneous group
            ind_lead=find(sum(strcmp(ser_sel(ones(size(ser_lead,1),1),:),ser_lead(:,ones(1,size(ser_sel,2)))),1)); %index of included series defining the leading group
        elseif ~preclas
            ser_cont=[]; %series falling into the contemporaneous group
            ser_lead=[]; %series falling into the leading group
            ind_cont=[]; %index of included series defining the contemporaneous group
            ind_lead=[];
        end
        %PP
    elseif strcmp(file,'pp');
        ser_excl={'price_products'}
        ser_sel = ser(find(1-sum(strcmp(ser(ones(size(ser_excl,1),1),:),ser_excl(:,ones(1,size(ser,2)))),1))); % string array of included series
        ser_ind=find(1-sum(strcmp(ser(ones(size(ser_excl,1),1),:),ser_excl(:,ones(1,size(ser,2)))),1)); %choose which countries to include, index array
        ser_cont=[];%{'YER';'PCR';'ITR';'MTR';'XTR'}; %series falling into the contemporaneous group
        ser_lead=[];%{'manuf-oscd';'KTPROL';'QTAUF';'EINDSE'}; %series falling into the leading group
        ind_cont=[];
        ind_lead=[];
        %BDS
    elseif strcmp(file,'bds');
        ser_excl={'price_products'};
        ser_sel = ser(find(1-sum(strcmp(ser(ones(size(ser_excl,1),1),:),ser_excl(:,ones(1,size(ser,2)))),1))); % string array of included series
        ser_ind=find(1-sum(strcmp(ser(ones(size(ser_excl,1),1),:),ser_excl(:,ones(1,size(ser,2)))),1)); %choose which countries to include, index array
        ser_cont={'PIB'};%'PCR';'ITR';'MTR';'XTR'}; %series falling into the contemporaneous group
        ser_lead={'manufacturing_index'};%'KTPROL';'QTAUF';'EINDSE'}; %series falling into the leading group
        ind_cont=find(sum(strcmp(ser_sel(ones(size(ser_cont,1),1),:),ser_cont(:,ones(1,size(ser_sel,2)))),1)); %index of included series defining the contemporaneous group
        ind_lead=find(sum(strcmp(ser_sel(ones(size(ser_lead,1),1),:),ser_lead(:,ones(1,size(ser_sel,2)))),1)); %index of included series defining the leading group

    end
    ser_diff={}; %series falling into the contemporaneous group
    ser_lev={'HICP-IG';'HICP';'YER';'MTR';'PCR';'ITR';'XTR';'ser-capa';'manuf-tppa';'manuf-ossk';'bat-apa';'bat-epa';...
        'bat-tuc';'HICP-FO';'HICP-E';'CLIMAT-FR';'CLIMAT-FR-EMPL';'IPI-CZ';...
        'ser-capre';'ser-dem';'manuf-pgp';'manuf-tppre';'manuf-oscde';'manuf-oscd';'bat-apre';'bat-jcc'};
    ind_lev=[];
    ind_levex=ones(1,size(ser,2))==1;

elseif strcmp(file,'uni') % Supprimer surement
    ser_excl={'HICP-IG';'HICP'};
    ser_sel = ser(find(1-sum(strcmp(ser(ones(size(ser_excl,1),1),:),ser_excl(:,ones(1,size(ser,2)))),1))); % string array of included series
    ser_ind=find(1-sum(strcmp(ser(ones(size(ser_excl,1),1),:),ser_excl(:,ones(1,size(ser,2)))),1)); %choose which countries to include, index array
    ser_cont={'YER'};%series falling into the contemporaneous group
    ser_lead=[]; %series falling into the leading group
    ind_cont=find(sum(strcmp(ser_sel(ones(size(ser_cont,1),1),:),ser_cont(:,ones(1,size(ser_sel,2)))),1));
    ind_lead=[]; %index of included series defining the leading group
    ser_diff={};
    ind_lev=[];
    ind_levex=ones(1,size(ser,2))==1;
end


cat_dum=data(1,:);
cat_dum=cat_dum(ser_ind);
data=data(2:end,:);
data_raw=data;


nobs=size(data,1);
cend = 2022.75;
cal = [sort(cend-[1:nobs]/4)+1/4]';
% cal_beg=cal(1);cal_end=cal(end);


% lrmean_smpl=find((cal>=cal_beg).*(cal<cal_end+1)); %choose sample period to compute the long-run mean of the series
obs_per=find((cal>=cal_beg).*(cal<=cal_end_est)); %choose sample period, contains the respective rows of original matrix
% lrmean_smpl=find((cal>=1980).*(cal<2003)); %choose sample period to compute the long-run mean of the series
% obs_per=find((cal>=1980).*(cal<2003)); %choose sample period, contains the respective rows of original matrix
cal=cal(obs_per);

nbank=length(ser_ind);
S=zeros(1,nbank);


% dlo=(log(data(2:end,:))-log(data(1:end-1,:)))*100; %endogenous variable, GDP quarterly growth rate
dlo=data(:,ser_ind);
data_raw=data_raw(:,ser_ind);
% Zlogit=Zlogit(ser_ind,:);

%dlo=dlo(obs_per(1:end),:); %ev. end-1 if data is differenced on line 51
%gdpyearlyrate=gdpyearlyrate(obs_per(1:end-3));
dlo_unadj=dlo;
nobs=size(dlo,1);
% cal=cal(2:end); %ev. uncomment if data is differenced on line 51
% dlo_obs=dlo;



nan_ret=find(1-((mean(dlo,1)>-Inf).*(mean(dlo,1)<Inf))); %find series with NaN  
n_ret=zeros(nobs,nbank);
for j=nan_ret;
   n_ret([find(1-((dlo(:,j)>-Inf).*(dlo(:,j)<Inf)))],j)=1; %construct matrix of missing values for NaN
end

dlo_out=(1-(1-n_ret).*(1))==1;

%compute mean and std without outliers
for j=1:nbank
    lrmean(j)=mean(dlo(dlo_out(:,j)==0,j));
    lrstd(j)=std(dlo(dlo_out(:,j)==0,j));
%    dlo([find(dlo_out(:,j))],j)=mean(dlo(dlo_out(:,j)==0,j));
end

% lrmean=mean(dlo); %ev. end-1 if data is differenced on line 51
% lrstd=std(dlo,0,1);

%Commentaires AQLT
%if ~unit_spec_var
%    dlo=(dlo-kron(ones(nobs,1),lrmean));%./kron(ones(nobs,1),lrstd);
%elseif and(unit_spec_var,nbank>2)
%    dlo=(dlo-kron(ones(nobs,1),lrmean))./kron(ones(nobs,1),lrstd);
%elseif and(unit_spec_var,nbank==2)
%    dlo=dlo-kron(ones(nobs,1),lrmean);
%end    
dlo_obs=dlo;

if sg_change
    % [sg] = max_conc_sign(dlo,ser_sel,{'YER'},4); %change sign for countercyclical variables according to the sign of largest concordance
    
    [sg] = neg_corr_sign(dlo,ser_sel,{'PIB'}); %change sign for series significantly negatively correlated with GDP
    %dlo=dlo.*sg(ones(nobs,1),:);%Commentaires AQLT

    %dlo_obs=dlo_obs.*sg(ones(nobs,1),:);
    %dlo_unadj=dlo_unadj.*sg(ones(nobs,1),:);
    %data_raw=data_raw.*sg(ones(size(data_raw,1),1),:);
    if K>1
    [sg] = neg_corr_sign(dlo,ser_sel,{'manufacturing_index'}); %change sign for series significantly negatively correlated with manuf-oscd
    %dlo=dlo.*sg(ones(nobs,1),:);%Commentaires AQLT
    %dlo_obs=dlo_obs.*sg(ones(nobs,1),:);
    %dlo_unadj=dlo_unadj.*sg(ones(nobs,1),:);
    %data_raw=data_raw.*sg(ones(size(data_raw,1),1),:);
    end
end

if S_logit
    corr_c={'PIB'}; %correlation with GDP
    c_ind=find(sum(strcmp(ser_sel(ones(size(corr_c,1),1),:),corr_c(:,ones(1,size(ser_sel,2)))),1));
    corr_l={'manufacturing_index'}; %correlation with orders
    l_ind=find(sum(strcmp(ser_sel(ones(size(corr_l,1),1),:),corr_l(:,ones(1,size(ser_sel,2)))),1));
    [R P]=corrcoef(dlo);
    Zlogit=[R(c_ind,:)' R(l_ind,:)'];
end



% Zbasis=zeros(nobs,[],nbank);
% Zbasis(:,1,:)=reshape(lenddata_b(:,9),nobs,nbank)*100; % zinssatz dir
% Zbasis(:,2,:)=reshape(lenddata_b(:,3),nobs,nbank); %Gr?sse si

% Zdum=zeros(nobs,4,nbank);%dummies 1-3: saisonals 4: structural break in der liquidit?t ab 4.quartal 1995 
% for i1=1:4;
%    Zdum(:,i1,:)=reshape(lenddata_b(:,4+i1),nobs,nbank); 
% end;clear i1;   

% Zwirt=zeros(nobs,2,nbank); % wirtschaftsdaten 1: dp86 (inflation) 2: dyr (wachstum bip)
% for i1=1:2;
%    Zwirt(:,i1,:)=reshape(lenddata_b(:,9+i1),nobs,nbank)*100; 
% end;clear i1;   
% dyr=lenddata_b(1:nobs,9+2)*100;
c=ones(size(dlo,1),1,size(dlo,2)); % constante 



clear outl n_ret 