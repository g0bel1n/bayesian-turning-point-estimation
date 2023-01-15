nst=size(etaMSmc,2);
Keta=size(etaMSmc,3);
if M< 500
    indexmc=[1:M];
else
    ispace=fix(M/500);indexmc=[1:ispace:M];
end   

str=['r';'k';'b';'y';'g';'c'];
leg=['group 1'; 'group 2';'group 3';'group 4'];
xlab={'\xi_{11}';'\alpha';'\beta';'\xi_{22}'};

figure(double(gcf)+1);
for i=1:nst
%     for j=1:nst
%         if i==j   
            subplot(2,2,i)  
            [xa,fa]=dichte(squeeze(etaMSmc(i,i,1,indexmc)));
            plot(xa,fa,str(k));hold on;
            xlabel(char(xlab(i,1)));
%         else
%             subplot(2,1,2)   
%             for k=1:Keta
%                 [xa,fa]=dichte(squeeze(etaMSmc(i,j,k,indexmc)));
%                 plot(xa,fa,str(k));hold on;
%             end
%             xlabel(['\eta_{i,j}']);
%         end
%     end
end
