
MC=size(gammc,1);
N=size(Zlogit,1);
nZ=size(Zlogit,2);
K=size(gammc,3);
prob_prior=zeros(N,K,MC);
prob_post=zeros(N,K,MC);

for j=1:MC
    gamp=[prgamm(:,ones(1,K-1))+inv(prgaminf)*randn(nZ,K-1) zeros(nZ,1)];
    prob_prior(:,:,j)=exp(Zlogit*gamp)./kron(ones(1,K),sum(exp(Zlogit*gamp),2));
    prob_post(:,:,j)=exp(Zlogit*permute(gammc(j,:,:),[2 3 1]))./kron(ones(1,K),sum(exp(Zlogit*permute(gammc(j,:,:),[2 3 1])),2));
end

figure(double(gcf)+1);
for j=1:K
    subplot(K,1,j)
bar(1:N,[mean(prob_post(:,j,:),3) mean(prob_prior(:,j,:),3)])
if j==1
    title('prior versus posterior group probabilities')
end
end


Z_logit_b=[-1.5:0.01:1.5]';
PP=length(Z_logit_b);
Z_logit_pp=[ones(PP,1) Z_logit_b zeros(PP,1)];
prob_prior_pp=zeros(PP,K,MC);
for j=1:MC
    prob_prior_pp(:,:,j)=exp(Z_logit_pp*permute(gammc(j,:,:),[2 3 1]))./kron(ones(1,K),sum(exp(Z_logit_pp*permute(gammc(j,:,:),[2 3 1])),2));
%    prob_post(:,:,j)=exp(Zlogit*permute(gammc(j,:,:),[2 3 1]))./kron(ones(1,K),sum(exp(Zlogit*permute(gammc(j,:,:),[2 3 1])),2));
end
[Zlogit_GDP,Zind_GDP]=sort(Zlogit(:,2));
figure(double(gcf)+1);
subplot(2,1,1)
plot(Z_logit_pp(:,2),mean(prob_prior_pp(:,1,:),3))
title('P(S_i=1)|\cdot, depending on corr with GDP')
subplot(2,1,2)
bar(1:N,mean(prob_post(Zind_GDP,1,:),3),0.1)

Z_logit_pp=[ones(PP,1) zeros(PP,1) Z_logit_b];
prob_prior_pp=zeros(PP,K,MC);
for j=1:MC
    prob_prior_pp(:,:,j)=exp(Z_logit_pp*permute(gammc(j,:,:),[2 3 1]))./kron(ones(1,K),sum(exp(Z_logit_pp*permute(gammc(j,:,:),[2 3 1])),2));
%    prob_post(:,:,j)=exp(Zlogit*permute(gammc(j,:,:),[2 3 1]))./kron(ones(1,K),sum(exp(Zlogit*permute(gammc(j,:,:),[2 3 1])),2));
end
[Zlogit_KTAUF,Zind_KTAUF]=sort(Zlogit(:,3));
figure(double(gcf)+1);
subplot(2,1,1)
plot(Z_logit_pp(:,3),mean(prob_prior_pp(:,2,:),3))
title('P(S_i=2)|\cdot, depending on corr with KTAUF')
subplot(2,1,2)
bar(1:N,mean(prob_post(Zind_KTAUF,2,:),3),0.1)

[Z1,Z2]=meshgrid(-1.0:0.1:1.0);
P1p=0;P2p=0;
for j=1:MC
%P1=exp(mean(gammc(j,1,1),1)+Z1*mean(gammc(:,2,1),1)+Z2*mean(gammc(:,3,1),1));
%P2=exp(mean(gammc(j,1,2),1)+Z1*mean(gammc(:,2,2),1)+Z2*mean(gammc(:,3,2),1));
P1=exp(gammc(j,1,1)+Z1*gammc(j,2,1)+Z2*gammc(j,3,1));
P2=exp(gammc(j,1,2)+Z1*gammc(j,2,2)+Z2*gammc(j,3,2));

P1p=P1p+P1./(1+P1+P2);
P2p=P2p+P2./(1+P1+P2);
end
P1p=P1p./MC;
P2p=P2p./MC;
% figure(gcf+1)
% subplot(1,2,1)
% mesh(Z1,Z2,P1p);
% xlabel('corr with GDP')
% ylabel('corr with KTAUF')
% zlabel('P(S_i=1|Z_i)')
% subplot(1,2,2)
% mesh(Z1,Z2,P2p);
% xlabel('corr with GDP')
% ylabel('corr with KTAUF')
% zlabel('P(S_i=2|Z_i)')

figure(double(gcf)+1);
%subplot(1,2,1)
%plot3(Z1,Z2,P1p);
mesh(Z1,Z2,P1p);
xlabel('corr with GDP')
ylabel('corr with KTAUF')
zlabel('P(S_i=1|Z_i)')
hold on
%plot3(Zlogit_GDP,Zlogit(Zind_GDP,3),mean(prob_post(Zind_GDP,2,:),3),'.','MarkerSize',15)
%pp=cat(3,[Zlogit_GDP,Zlogit(Zind_GDP,3),zeros(N,1)],[Zlogit_GDP,Zlogit(Zind_GDP,3),mean(prob_post(Zind_GDP,1,:),3)]);
pp=cat(3,[Zlogit(:,2),Zlogit(:,3),zeros(N,1)],[Zlogit(:,2),Zlogit(:,3),mean(Smc(1,:,:),3)']);
pp1=permute(pp,[3 2 1]);
for j=1:size(pp1,3)
    plot3(pp1(:,1,j),pp1(:,2,j),pp1(:,3,j),'b-h','MarkerSize',5)
    hold on
end
hidden off

figure(double(gcf)+1);
%subplot(1,2,2)
%plot3(Z1,Z2,P2p);
mesh(Z1,Z2,P2p);
xlabel('corr with GDP')
ylabel('corr with KTAUF')
zlabel('P(S_i=2|Z_i)')
hold on
%plot3(Zlogit_GDP,Zlogit(Zind_GDP,3),mean(prob_post(Zind_GDP,2,:),3),'.','MarkerSize',15)
pp=cat(3,[Zlogit(:,2),Zlogit(:,3),zeros(N,1)],[Zlogit(:,2),Zlogit(:,3),mean(Smc(2,:,:),3)']);
pp1=permute(pp,[3 2 1]);
for j=1:size(pp1,3)
    plot3(pp1(:,1,j),pp1(:,2,j),pp1(:,3,j),'b-h','MarkerSize',5)
    hold on
end
hidden off






