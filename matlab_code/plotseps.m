str=['r';'k';'b';'y';'g';'c'];
set(0,'DefaultLineLineWidth',1.5)

if unit_spec_var
    sepsm=kron(ones(N,1),mean(sepsmc,1));
    sepsi=mean(sepsmc(:,ones(1,N))./lambmc,1);
    [ssepsi si]=sort(sepsi);
end


gh=figure(double(gcf)+1);
[xa,fa]=dichte(sepsmc);
if unit_spec_var
    scatter(sepsi',sepsm*10,4,'o',str(k))
    hold on
%     for j=[1:N];
%         text(ssepsi(j),sepsm(1)+0.15*(N-j+1),ctries(si(j),:),'VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',12)
%     end
end
plot(xa,fa,str(k));
set(gca,'FontSize',14)
%title(['marginal of \sigma^2, \sigma^2/\lambda_i'])
print(gh,'-deps',[file '_sepsi']);



