function [sig] = neg_corr_sign(data,ser_str,corr_ser)

T=size(data,1);
N=size(data,2);

[R P]=corrcoef(data);

corr_ind=find(sum(strcmp(ser_str(ones(size(corr_ser,1),1),:),corr_ser(:,ones(1,size(ser_str,2)))),1));

sig  = 1-(((R(corr_ind,:)<0).*(P(corr_ind,:)<=0.05)))+(((R(corr_ind,:)<0).*(P(corr_ind,:)<=0.05))).*(-1);

