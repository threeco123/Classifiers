function DS = Bayes_DSurface(H, target)

H0 = H(find(target == 0),:);
H1 = H(find(target == 1),:);

u0 = mean(H0);
u1 = mean(H1);

Cov0 = cov(H0);
Cov1 = cov(H1);

P0 = size(H0,1)/(size(H0,1));
P1 = size(H1,1)/(size(H,1));


x1 = linspace(min(H(:,1))*1.2,max(H(:,1))*1.2,251);
x2 = linspace(min(H(:,2))*1.2,max(H(:,2))*1.2,251);
[xTest1,xTest2] = meshgrid(x1,x2);
xTest= [xTest1(:) xTest2(:)];

for i = 1:size(xTest,1)
    DS(i) = (-0.5*(xTest(i,:) - u1)*inv(Cov1)*(xTest(i,:) - u1)' - 0.5*log(det(Cov1)) + log(P1))...
        - (-0.5*(xTest(i,:) - u0)*inv(Cov0)*(xTest(i,:) - u0)' - 0.5*log(det(Cov0)) + log(P0));
end

DS= reshape(DS,length(x2),length(x1));

imagesc(x1([1 end]),x2([1 end]),DS);

 hold on 
 plot(H0(:,1),H0(:,2),'b*'); % H0

 plot(H1(:,1),H1(:,2),'ro') % H1
 
end