function dataperiod = BinPeriod(dataperiod,thr_ipi)
% function dataperiod = BinPeriod(dataperiod,ipi)
% Combine the periods when the inter period interval is smaller than thr_ipi 
if size(dataperiod,1)<2
    return
end
tmp_end = find(dataperiod(1:end-1,2)>=dataperiod(2:end,1));
dataperiod(tmp_end,2) = dataperiod(tmp_end+1,1)-1;
be_tmp = sortrows([dataperiod(:),reshape(repmat(1:2,size(dataperiod,1),1),[],1)]);
dtmp = find((be_tmp(2:end,1)-be_tmp(1:end-1,1))<=thr_ipi & be_tmp(1:end-1,2)==2);
be_tmp([dtmp;dtmp+1],:)=[];
dataperiod = reshape(be_tmp(:,1),2,[])';