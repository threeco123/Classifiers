function Beta_new = LD_DSurface(H, target)

% xH = zeros(length(H), size(H,2)+1);
% xH(:,1) = 1;
% xH(:,2:(size(H,2)+1)) = H;

H0 = H(find(target == 0),:);
H1 = H(find(target == 1),:);

u0 = mean(H0);
u1 = mean(H1);

Cov0 = cov(H0);
Cov1 = cov(H1);

CovW = Cov0 + Cov1;

sigmoidFcn = @(x)1./(1+exp(-x));

Lamd = inv(CovW)*(u1-u0)';

xTest = zeros(length(H), (size(H,2)+1));
xTest(:,1) = 1;
xTest(:,2:(size(H,2)+1)) = H;
LamdN = [0 Lamd'];


Beta_old = zeros(1,size(H,2)+1);
Beta_new = LamdN;
count = 0;

while(norm(Beta_new - Beta_old)> 0.01)

Beta_old = Beta_new;

P = xTest*Beta_old';
L = 1./(1+exp(-P));
% for i = 1:size(xTest,1)
%       P(i) = Beta_old*xTest(i,:)';
% end
w = zeros(length(xTest),length(xTest));

for i = 1:size(xTest,1)
      w(i,i) = L(i)*(1-(L(i)));
end

Beta_new = Beta_old + (inv(xTest'*w*xTest)*xTest'*(target - L))';
% test_num = norm(Beta_new - Beta_old)
count = count + 1;
end


% x1 = linspace(min(H(:,1))*1.2,max(H(:,1))*1.2,251);
% x2 = linspace(min(H(:,2))*1.2,max(H(:,2))*1.2,251);
% [xTest1,xTest2] = meshgrid(x1,x2);
% xTest= [xTest1(:) xTest2(:)];
% 
% xTestN = zeros(length(xTest), size(xTest,2)+1);
% xTestN(:,1) = 1;
% xTestN(:,2:(size(xTest,2)+1)) = xTest;
% 
% for i = 1:size(xTestN,1)
%       DS(i) = Lamd*xTestN(i,:)';
% end
% 
% DS= reshape(DS,length(x2),length(x1));
% 
% imagesc(x1([1 end]),x2([1 end]),DS);
% 
%  hold on 
%  plot(H0(:,1),H0(:,2),'b*'); % H0
% 
%  plot(H1(:,1),H1(:,2),'ro') % H1

% Beta_new = [-0.2183    0.8266    0.6405];
 
 x1 = linspace(min(H(:,1))*1.2,max(H(:,1))*1.2,251);
x2 = linspace(min(H(:,2))*1.2,max(H(:,2))*1.2,251);
[xTest1,xTest2] = meshgrid(x1,x2);
xTest= [xTest1(:) xTest2(:)];

xTestN = zeros(length(xTest), size(xTest,2)+1);
xTestN(:,1) = 1;
xTestN(:,2:(size(xTest,2)+1)) = xTest;

for i = 1:size(xTestN,1)
      Num = Beta_new*xTestN(i,:)';
      DS(i) = sigmoidFcn(Num);
end

DS= reshape(DS,length(x2),length(x1));

imagesc(x1([1 end]),x2([1 end]),DS);
title('LD Decision Surface');
% set(gca,'CLim',[-1 1])

 hold on 
 plot(H0(:,1),H0(:,2),'b*'); % H0

 plot(H1(:,1),H1(:,2),'ro') % H1
 legend('H0', 'H1')
