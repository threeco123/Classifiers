function [pd pf] = Bayes_Classifier(H0, H1, Ht0, Ht1)

% H0 = H(find(target == 0),:);
% H1 = H(find(target == 1),:);

u0 = mean(H0);
u1 = mean(H1);

Cov0 = cov(H0);
Cov1 = cov(H1);

P0 = size(H0,1)/(size(H0,1) + size(H1,1));
P1 = size(H1,1)/(size(H0,1) + size(H1,1));

xTest = [Ht0;Ht1];

for i = 1:size(xTest,1)
    DS(i) = (-0.5*(xTest(i,:) - u1)*inv(Cov1)*(xTest(i,:) - u1)' - 0.5*log(det(Cov1)) + log(P1))...
        - (-0.5*(xTest(i,:) - u0)*inv(Cov0)*(xTest(i,:) - u0)' - 0.5*log(det(Cov0)) + log(P0));
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

end