figure(double(gcf)+1);

if M< 500
   indexmc=[1:M];
else
   ispace=fix(M/500);indexmc=[1:ispace:M];
end   


for i1=1:dMS
   subplot(nv,nh,i1);
   for k=1:K-(K>2)
   [xa,fa]=dichte(alphamc(indexmc,indexMS((k-1)*dMS+i1)));
   plot(xa,fa,str(k));hold on;
   xlabel(['\beta^R_{' num2str(i1) '}']);
   if i1==dd
       axis([-2 2 min(fa) max(fa)])
   end
   
   end
   if i1==1 title('MS parameter');legend(leg(1:k,:));end
end

