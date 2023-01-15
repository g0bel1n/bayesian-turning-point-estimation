%str=['r';'k';'b';'y';'g';'c'];
str=['k- ';'k- ';'b  ';'y  ';'g  ';'k--'];
leg=['group 1'; 'group 2';'group 3';'group 4'];
qdim=size(Q,1);
M=size(Qmc,1);
if M< 500
    indexmc=[1:M];
else
    ispace=fix(M/750);indexmc=[1:ispace:M];
end
Qsim=zeros(size(Qmc(:,:,1)));
for j=1:M
    Q=raninvwi_neu(prQnu,prQS);
    Qsim(j,:)=qincol(Q)';
end
qno=0;
figure(double(gcf)+1);
for j=1:size(Q,1)
    for l=1:j
        subplot(qdim,qdim,l+qdim*(j-1));
        qno=qno+1;
        [xa,fa]=dichte(Qsim(indexmc,qno));
        %plot(xa,fa,str(end));hold on;
        plot(xa,fa,'k--');hold on;
        for k=1:K-(K>2)
            [xa,fa]=dichte(Qmc(indexmc,qno,k));
            plot(xa,fa,str(k));hold on;
        end
        if (l+qdim*(j-1))==(qdim^2-fix(qdim/2))
            legend(['prior  ';leg(1:k,:)]);
        end
    end
    if j==1
        title(['Q of group-specific parameters, ' int2str(cal_beg) '-' int2str(cal_end)]);
    end


end



% var=[];
% iplot=0;
% if dMS>0;
%     figure(gcf+1);
%     if fix(dMS/2)>=1;
%         for i=1:fix(dMS/2)
%             i1=sum(1:1+dd);i2=sum(1:2+dd);ih=4;iplot=iplot+1;
%             if dMS>2
%                 subplot(fix((fix(dMS/2)+ih-1)/ih),ih,iplot);
%             end
%             for k=1:K
%                 %scatter(alphamc(:,dd*K+i1),etamc(:,i),1.,str(i));
%                 scatter(Qmc(indexmc,i1,k),Qmc(indexmc,i2,k),1.,str(k));
%                 %var=[var alphamc(indexmc,dd*(k-1)+i1),alphamc(indexmc,dd*(k-1)+i2)];
%                 %save scatG_K2.asc
%                 xlabel(['Q_{' num2str(i1) num2str(i1) '}']);
%                 ylabel(['Q_{' num2str(i2) num2str(i2) '}']);
%                 hold on;
%             end
%             if i==1 title('Q of switching parameters'); end
%         end
%     else
%         i1=sum(1:iiMIX(1));i2=sum(1:1+dd);ih=4;iplot=iplot+1;
%         for k=1:K
%             %scatter(alphamc(:,dd*K+i1),etamc(:,i),1.,str(i));
%             scatter(Qmc(indexmc,i1,k),Qmc(indexmc,i2,k),1.,str(k));
%             %var=[var alphamc(indexmc,dd*(k-1)+i1),alphamc(indexmc,dd*(k-1)+i2)];
%             %save scatG_K2.asc
%             xlabel(['Q_{' num2str(i1) num2str(i1) '}']);
%             ylabel(['Q_{' num2str(i2) num2str(i2) '}']);
%             hold on;
%         end
%         title('Q of switching parameters'); 
%     end
% end
