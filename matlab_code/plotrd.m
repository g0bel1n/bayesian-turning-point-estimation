dMS=size(ZMS,2);
dd = size(Q,1)-dMS;
str=['r';'k';'b';'y';'g';'c'];
leg=['group 1'; 'group 2';'group 3';'group 4'];

indexG = kron(ones(1,K),[1:dd]); %reshape(index1(ones(1,K),:)',1,K*dd);
%index2 = [1:K]';
indexK = kron([1:K]',ones(dd,1))'; %reshape(states(:,ones(1,dd))',1,K*dd);
indexR = kron(ones(1,K),[1:dMS]);
indexM = kron([1:K]',ones(dMS,1))';
indexMIX=[1:dd*K];
indexFIX=[K*dd+1:d+(K-1)*dd];
indexMS=[(d+(K-1)*dd+1):(d+(K-1)*dd+dMS*K)];

nh=4;%number of horizental plots
nv=fix((dd+nh-1)/nh);
plotrdalpha;
title(['K' int2str(K) ', ' int2str(cal_beg) '-' int2str(cal_end)])
if dMS>0
    if size(etamc,2)==1
        ifig=ifig-1;
    end
    %str([1 2])=str([2 1]);
    plotrdalpha_MS;
    title(['K' int2str(K) ', ' int2str(cal_beg) '-' int2str(cal_end)])
    subplot(nv,nh,i1);
    
    plotrdalphaR;
    title(['K' int2str(K) ', ' int2str(cal_beg) '-' int2str(cal_end)])
    
end
nvfix=fix(((d-dd)+nh-1)/nh);
% plotrdalpha_fix;
% if size(etamc,2)>1
%   plotrdeta;
% end  
% if dMS>0
% plotrdetaMS;
  title(['K' int2str(K) ', ' int2str(cal_beg) '-' int2str(cal_end)])

% end

% ifig=ifig+1;
% figure(ifig);
% for k=1:K
%     [xa,fa]=dichte(istarmc(indexmc,k));
%     plot(xa,fa,str(k));hold on;
% end
%   title(['i^*, K' int2str(K) ', ' int2str(cal_beg) '-' int2str(cal_end)])

