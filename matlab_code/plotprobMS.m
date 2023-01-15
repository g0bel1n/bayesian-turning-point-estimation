%load testdata_itdiff
%load testdata_itequal
%load testdata_big
%[ss si]=sort(parms(:,1));
%IMS=IMS(si,:);
rec=2; %1 means below-average growth, 2 means recession
nst=2;
it=size(LIMSmc,2);
KIMS=size(LIMSmc,1);

%probMS=zeros(nst,it,K);
probMS=1-mean(LIMSmc,3);
if unit_spec_var
    y_gr=dlo(end-it+1:end,:)./ kron(ones(it,1),sqrt(mean(sepsmc(:,ones(1,size(lambmc,2)))./lambmc,1)));
elseif ~unit_spec_var
    y_gr=dlo(end-it+1:end,:)./ kron(ones(it,1),sqrt(mean(sepsmc(:,ones(1,size(lambmc,2))),1)));
end
  
% y_gr=dlo(end-it+1:end,:);
time = 1:it;
% c = 2003;
c=cal_end_est+1/4;
%c=it+1
cal = sort(c-time/4)';
%cal = sort(c-time)';
% if rec==1;
%     dyrrec1=(dyr(end-it+2:end)<mean(dyr(end-it+6:end)))*2;%
%     dyrrec2=(dyr(end-it+2:end)<mean(dyr(end-it+6:end)))*(-2);
% end
% if rec==2;
%    dyrrec1=((cal>1991.5).*(cal<1993.25))*2;
%    dyrrec2=((cal>1991.5).*(cal<1993.25))*(-2);
% end   
% CEPR_dat=(cal>1980).*(cal<=1982.5)+(cal>1992).*(cal<=1993.5);
% NBER_dat=(cal>1980).*(cal<=1980.5)+(cal>1981.5).*(cal<=1982.75)+(cal>1990.5).*(cal<=1991)+(cal>2001).*(cal<=2001.75)

% for j=1:KIMS;
%     ifig=ifig+1;
%     figure(ifig);
%     
%     for k=1:nst
%         probMS(k,:,j)=sum(permute(IMSmc(j,:,:)+1,[3 2 1])==k*ones(size(IMSmc,3),size(IMSmc,2)),1)/size(IMSmc,3);
%         %probMS(k,:)=sum((2-IMSmc)==k*ones(size(IMSmc)))/size(IMSmc,1);
%         subplot(nst+1,1,k);
%         bar(cal,probMS(k,:,j));
% %        hold on;
% %        plot(cal,IMS(j,end-it+1:end));
%         axis([cal(1) cal(end) 0 1]);
%         if k==1 title(['I_t=0, group ' int2str(j)]);end  
%         if k==2 title(['I_t=1, group ' int2str(j)]);end  
%     end
%     probMS(:,:,1)=1-probMS(:,:,1);
%     % subplot(nst+1,1,nst+1);
%     % bar(cal(end-length(dyrrec1)+1:end),dyrrec1,'y');
%     % %axis([cal(1) cal(end) -2 2]);
%     % hold on;
%     % bar(cal(end-length(dyrrec2)+1:end),dyrrec2,'y');
%     % hold on;
%     % plot(cal(2:end),dyr(end-it+2:end),'k-',cal(2:end),zeros(it-1,1),'k-');
%     % hold on;
%     % plot(cal(2:end),Zbasis(end-it+2:end,1,1),'k--',cal(2:end),zeros(it-1,1),'k-');
%     % axis([cal(1+lag_y) cal(end) -2 2]);
%     % hold on;
%     
%     %axis([cal(1) cal(end) -2 2]);
%     
%     
%     
%     
%     ifig=ifig+1;
%     figure(ifig);
%     probm=max(probMS(:,:,j));
%     for k=1:nst
%         subplot(nst+1,1,k);
%         bar(cal,probMS(k,:,j)==probm);
% %        hold on;
% %        plot(cal,IMS(j,end-it+1:end));
%         axis([cal(1) cal(end) 0 1]);
%         if k==1 title(['I_t=0, group ' int2str(j)]);end  
%         if k==2 title(['I_t=1, group ' int2str(j)]);end  
%     end
% end
% subplot(nst+1,1,nst+1);
% bar(cal(end-length(dyrrec1)+1:end),dyrrec1,'y');
% %axis([cal(1) cal(end) -2 2]);
% hold on;
% bar(cal(end-length(dyrrec2)+1:end),dyrrec2,'y');
% hold on;
% plot(cal(2:end),dyr(end-it+2:end),'k-',cal(2:end),zeros(it-1,1),'k-');
% axis([cal(1) cal(end) -2 2]);
%axis([cal(1) cal(end) -2 2]);

% ifig=ifig+1;
% figure(ifig)
% for k=1:K
%     subplot(K,1,k);
%     bar(cal,probMS(1,:,k));
%     hold on
%     plot(cal,ones(length(cal),1)*0.5);
%     axis([cal(1) cal(end) 0 1]);
%     if k==1 title(['I_t=0']);end  
%     
% end

figure(double(gcf)+1)
for k=1:(KIMS-(K>2)*(1-DStruc));
    subplot(KIMS-(K>2)*(1-DStruc),1,k);
%     maxk=max(max(max(y_gr(:,prob(k,:)==1)))+0.1,1);
%     maxk=1;
%     mink=min(min(y_gr(:,prob(k,:)==1)))-0.1;
    bar(cal,probMS(k,:),0.2);
    hold on
%     if or(k==1,k==3);
%         bar(cal,[(probMS(1,:,k)'.*CEPR_dat.*(prob(k,end)==0))+(probMS(1,:,k)'.*NBER_dat.*(prob(k,end)==1))],'m');
%     elseif k==2
%         bar(cal,[probMS(1,:,k)'.*NBER_dat],'m');
%     end
    hold on;
    if ~DStruc
        ylabel(['S_i=' int2str(k) ],'FontSize',12)
    elseif DStruc
        ylabel(['\rho=' num2str(mean(rhomc(:,k)),'%f5.2')])
    end
    ax1=gca;
    set(ax1,'XLim',[cal(1) cal(end)],'XTick',[cal_beg:2:cal_end],'YLim',[0 1],'YTick',[0:0.5:1])
%     ax2=axes('Position',get(ax1,'Position'),'XAxisLocation','top','YAxisLocation','right','Color','none');
%     set(ax2,'XLim',[cal(1) cal(end)],'XTickLabel',' ','YLim',...
%         [min(min(y_gr(:,prob(k,:)==1)))-0.5 max(max(y_gr(:,prob(k,:)==1)))+0.5])
% %     %     set(ax2,'XLim',[cal(1) cal(end)],'XTickLabel',' ','YLim',[-5 5])
% %     %     if k==1 set(ax2,'XLim',[cal(1) cal(end)],'YLim',[-5 5]); end
% %     %     if find(y_gr(:,prob(k,:)==1))
%     if any(prob(k,:))
%         line(cal,y_gr(:,prob(k,:)==1),'Parent',ax2,'Color','k')
%         hold on
%         line(cal,y_gr(:,1),'Parent',ax2,'Color','r')
% %       
% %         hold on
%     end
% %     line(cal,zeros(length(cal),1),'Parent',ax2,'Color','k')
%     line(cal,ones(length(cal),1)*mean(mean(y_gr(:,prob(k,:)==1))),'Parent',ax2,'Color','k')

    if k==1 title(['I_t=1']);end  
    
end

% if k>1
%     pr1=[cal/1000 squeeze(probMS(1,:,:))>0.5]
%     pr2=squeeze(probMS(1,:,:))
% elseif k==1
%     pr1=[cal/1000 (probMS(1,:))'>0.5]
%     pr2=(probMS(1,:))'
% end
