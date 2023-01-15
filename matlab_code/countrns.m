function ct = countrns(sst,nst)
%

sstd = diff(sst);
ct = zeros(nst);

for tr = [sst(1:end-1);sst(2:end)]
    ct(tr(1),tr(2)) = ct(tr(1),tr(2)) + 1;
end

