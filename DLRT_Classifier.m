function [pd pf] = DLRT_Classifier(H0, H1, Ht0, Ht1,k)

% H0 = H(find(target == 0),:);
% H1 = H(find(target == 1),:);

n0 = length(H0)/(length(H0) + length(H1));
n1 = length(H1)/(length(H0) + length(H1));

D = size(H0,2);

xTest = [Ht0;Ht1];

for len = 1:length(k)
for iHt = 1: length(xTest)
    
    Dist0 = Calculate_Norm(xTest(iHt,:), H0);
    Dist0 = sort(Dist0);
    
    Dist1 = Calculate_Norm(xTest(iHt,:), H1);
    Dist1 = sort(Dist1);
    DS(iHt) = log(n0/n1) + D*(log(Dist0(k(len))) - log(Dist1(k(len))));
end

DS0 = DS(1:length(Ht0));
DS1 = DS((length(Ht0)+1):length(DS));

DS0 = sort(DS0);
DS1 = sort(DS1);

pf(:,len) = 0:1/length(DS):1;

for count = 1:(length(DS)+1)
     Beta_Index = ceil(length(DS0) - pf(count)*length(DS0));
     if Beta_Index == 0
         pd(count,len) = 1;
     else
         Beta = DS0(Beta_Index);
         pd(count,len) = sum(find(DS1>Beta) > 0)/length(DS1);
     end
end



end
% V = max(DS) - min(DS);
% 
% for count = 1:length(DS)
%      k = min(DS) + count*V/length(DS)
%      pf(count) = sum(find(DS>k)<=size(Ht0,1))/size(Ht0,1);
%      pd(count) = sum(find(DS>k)>size(Ht0,1))/size(Ht1,1);
% end
end