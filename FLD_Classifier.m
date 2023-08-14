function [pd pf] = FLD_Classifier(H0, H1, Ht0, Ht1)

% H0 = H(find(target == 0),:);
% H1 = H(find(target == 1),:);

u0 = mean(H0);
u1 = mean(H1);

Cov0 = cov(H0);
Cov1 = cov(H1);

CovW = Cov0 + Cov1;

xTest = [Ht0;Ht1];

Lamd = inv(CovW)*(u1-u0)';

for i = 1:size(xTest,1)
    DS(i) = Lamd'*xTest(i,:)';
end

DS0 = DS(1:length(Ht0));
DS1 = DS((length(Ht0)+1):length(DS));

DS0 = sort(DS0);
DS1 = sort(DS1);

pf = 0:1/length(DS):1;

for count = 1:(length(DS)+1)
     Beta_Index = ceil(length(DS0) - pf(count)*length(DS0));
     if Beta_Index == 0
         pd(count) = 1;
     else
         Beta = DS0(Beta_Index);
         pd(count) = sum(find(DS1>Beta) > 0)/length(DS1);
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