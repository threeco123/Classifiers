function [Pd Pf] = LD_Classifier(H0,H1,Ht0,Ht1)

u0 = mean(H0);
u1 = mean(H1);

Cov0 = cov(H0);
Cov1 = cov(H1);

CovW = Cov0 + Cov1;

xTest = [Ht0;Ht1];

Lamd = inv(CovW)*(u1-u0)';

sigmoidFcn = @(x)1./(1+exp(-x));
xH0 = zeros(length(H0), (size(H0,2)+1));
xH0(:,1) = 1;
xH0(:,2:(size(H0,2)+1)) = H0;

xH1 = zeros(length(H1), (size(H1,2)+1));
xH1(:,1) = 1;
xH1(:,2:(size(H1,2)+1)) = H1;

xHt0 = zeros(length(Ht0), (size(Ht0,2)+1));
xHt0(:,1) = 1;
xHt0(:,2:(size(Ht0,2)+1)) = Ht0;

xHt1 = zeros(length(Ht1), (size(Ht1,2)+1));
xHt1(:,1) = 1;
xHt1(:,2:(size(Ht1,2)+1)) = Ht1;

xTest = [xHt0; xHt1];
xTrain = [xH0; xH1];
LamdN = [0 Lamd'];

Beta_old = zeros(1,size(H0,2)+1);
Beta_new = LamdN;

target = zeros(length(xTrain),1);
target((length(target)/2+1):length(target)) = 1;

coun = 0;
while(norm(Beta_new - Beta_old)> 0.01)

Beta_old = Beta_new;

P = xTrain*Beta_old';
L = 1./(1+exp(-P));

% for i = 1:size(xTrain,1)
%       P(i) = Beta_old*xTrain(i,:)';
% end

w = zeros(length(xTrain),length(xTrain));

for i = 1:size(xTrain,1)
      w(i,i) = L(i)*(1-L(i));
end

Beta_new = Beta_old + (inv(xTrain'*w*xTrain)*xTrain'*(target - L))';
coun = coun + 1;
end

for i = 1:size(xTest,1)
      Num = Beta_new*xTest(i,:)';
      DS(i) = sigmoidFcn(Num);
end

DS0 = DS(1:length(Ht0));
DS1 = DS((length(Ht0)+1):length(DS));

DS0 = sort(DS0);
DS1 = sort(DS1);

Pf = 0:1/length(DS):1;

for count = 1:(length(DS)+1)
     Beta_Index = ceil(length(DS0) - Pf(count)*length(DS0));
     if Beta_Index == 0
         Pd(count) = 1;
     else
         Beta = DS0(Beta_Index);
         Pd(count) = sum(find(DS1>Beta) > 0)/length(DS1);
     end
end











