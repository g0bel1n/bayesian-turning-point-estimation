dMS=size(ZMS,2);
index_col=[1:size(ZG,2)];
index_perm2=zeros(dMS,1);
for iMS=1:dMS
   incol=index_col(all(all((ZG-ZMS(:,iMS*ones(size(ZG,2),1),:))==0,3)==1,1));
   if isempty(incol) 
      ['warning: column ' num2str(iMS) ' of ZMS not contained in ZG']
   else   
      index_perm2(iMS)=incol;
   end  

end   