K = size(etamc,2);
M = size(etamc,1);
dMS=size(ZMS,2);
dd=size(Q,1)-dMS;
str=['k- ';'k--'];%'b';'y';'g';'c'};
set(0,'DefaultLineLineWidth',1.5)
Qdim=size(Q,1);

index1 = [1:dd];
index = reshape(index1(ones(1,K),:)',1,K*dd);
states = [1:K]';
index2 = reshape(states(:,ones(1,dd))',1,K*dd);
indexMIX=[1:dd*K];
indexFIX=[K*dd+1:d+(K-1)*dd];
indexMS=[(d+(K-1)*dd+1):(d+(K-1)*dd+dMS*K)];

am=mean(alphamc,1);
betai=squeeze(mean(bmc,1));
for k=1:K
    betai(prob(k,:)==1,:)=betai(find(prob(k,:)),:)+am(ones(sum(prob(k,:),2),1),[dd*(k-1)+[1:dd] indexMS(dMS*(k-1)+[1:dMS])]);
end
amc=permute(reshape(alphamc(:,indexMIX),M,dd,K),[1 3 2]);
amc=cat(3,amc,permute(reshape(alphamc(:,indexMS),M,dMS,K),[1 3 2]));
D=states(:,ones(1,size(prob1,2)));
D=D(prob1(1:K,:)==kron(ones(K,1),max(prob1,[],1)))'; %dem maximum zuordnen

bimc=mean(bmc+amc(:,D,:),1);

uc_am=amc(:,:,1)./(1-sum(amc(:,:,2:3),3));
uc_am=cat(3,uc_am,(amc(:,:,1)-amc(:,:,end))./(1-sum(amc(:,:,2:3),3)));
uc_bi=(bmc(:,:,1)+amc(:,D,1))./(1-sum(bmc(:,:,2:3)+amc(:,D,2:3),3));
uc_bi=mean(cat(3,uc_bi,(bmc(:,:,1)+amc(:,D,1)-bmc(:,:,end)-amc(:,D,end))./(1-sum(bmc(:,:,2:3)+amc(:,D,2:3),3))),1);
uc_mbi=mean(bmc(:,:,1)+amc(:,D,1))./mean(1-sum(bmc(:,:,2:3)+amc(:,D,2:3),3));
uc_mbi=cat(3,uc_mbi,mean(bmc(:,:,1)+amc(:,D,1)-bmc(:,:,end)-amc(:,D,end))./mean(1-sum(bmc(:,:,2:3)+amc(:,D,2:3),3)));
stat_bi=(sum(bmc(:,:,2:3)+amc(:,D,2:3),3));

% for j=1:Qdim
%     figure(gcf+1)
%     for k=1:K;
%         [xb,fb]=dichte(betai(prob(k,:)==1,j));
%         if j<=dd
%             [xa,fa]=dichte(alphamc(:,(k-1)*dd+j));
%         elseif  j>dd
%             [xa,fa]=dichte(alphamc(:,indexMS((k-1)*dMS+j-dd)));
%         end
%         plot(xa,fa,str(k),xb,fb,str(k));
%         hold on
%     end
%     title(['beta_i vs. beta_k, parameter ' int2str(j)])
% end
% 
% for j=1:Qdim
%     figure(gcf+1)
%     for k=1:K;
%         xb=betai(prob(k,:)==1,j);
%         if j<=dd
%             xa=am(ones(1,size(xb,1)),(k-1)*dd+j);
%         elseif  j>dd
%             xa=am(ones(1,size(xb,1)),indexMS((k-1)*dMS+j-dd));
%         end
%         scatter(xa,xb,4);
%         hold on
%     end
%     title(['beta_i vs. beta_k, parameter ' int2str(j)])
% end

% for j=1:Qdim
%     figure(gcf+1)
%     for k=1:K-(K>2);
%         xb=betai(prob(k,:)==1,j);
%         xbb=bimc(1,prob(k,:)==1,j);
%         if j<=dd
%             [xa,fa]=dichte(alphamc(:,(k-1)*dd+j));
%             xam=am(ones(1,size(xb,1)),(k-1)*dd+j);
%         elseif  j>dd
%             [xa,fa]=dichte(alphamc(:,indexMS((k-1)*dMS+j-dd)));
%             xam=am(ones(1,size(xb,1)),indexMS((k-1)*dMS+j-dd));
%         end
%         scatter(xb,xam,4,'o',str(k))
%         hold on
%         scatter(xbb,xam,4,'*',str(k))
%         hold on
%         plot(xa,fa,str(k));
%         hold on
%     end
% title(['beta_i vs. beta_k, parameter ' int2str(j)])
% for k=1:K-(K>2);
%     figure(gcf+1)
%     for l=find(prob(k,:)==1);
%     [xa,fa]=dichte(bmc(:,l,j)+amc(:,k,j));
%     plot(xa,fa,str(k));
% %     if or(j==2,j==3)
% %         [j l]
% %         pause
% %     end
%     hold on
%     end
%     title(['distribution of parameter ' int2str(j) ' for group ' int2str(k)])
% end

% end

for j=1:nst;
gh=   figure(double(gcf)+1)
    for k=1:K-(K>2);
        xb=uc_mbi(1,prob(k,:)==1,j);
            [xa,fa]=dichte(uc_am(:,k,j));
            xam=kron(ones(1,size(xb,2)),2*k*abs(mean(uc_am(:,k,j))));
        scatter(xb,xam,4,'o',str(k))
        hold on
        plot(xa,fa,str(k,:));
        set(gca,'FontSize',14)
        hold on
    end
title(['unconditional mean: state ' int2str(3-j)])
print(gh,'-deps',[file '_uncmeanst' int2str(j)]);
end

gh=   figure(double(gcf)+1)
    for k=1:K-(K>2);
        xb=mean(stat_bi(:,prob(k,:)==1),1);
            [xa,fa]=dichte(sum(amc(:,k,2:3),3));
            xam=kron(ones(1,size(xb,2)),2*abs(mean(sum(amc(:,k,2:3),3))));
        scatter(xb,xam,4,'o',str(k))
        hold on
        plot(xa,fa,str(k,:));
        set(gca,'FontSize',14)
        hold on
    end
title(['sum of autoregressive coefficients '])
print(gh,'-deps',[file '_sumauto']);


% for k=1:K-(K>2);
%     figure(gcf+1)
%     for j=find(prob(k,:)==1);
%     [xa,fa]=dichte(stat_bi(:,j));
%     plot(xa,fa,str(k));
%     hold on
%     end
%     title(['distribution of sum of autoregressive coefficients for group ' int2str(k)])
% end




