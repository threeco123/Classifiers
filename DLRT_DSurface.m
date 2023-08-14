function DS = DLRT_DSurface(H, target,k)

H0 = H(find(target == 0),:);
H1 = H(find(target == 1),:);

n0 = length(H0)/(length(H0) + length(H1));
n1 = length(H1)/(length(H0) + length(H1));

D = size(H,2);

x1 = linspace(min(H(:,1))*1.2,max(H(:,1))*1.2,251);
x2 = linspace(min(H(:,2))*1.2,max(H(:,2))*1.2,251);
[xTest1,xTest2] = meshgrid(x1,x2);
xTest= [xTest1(:) xTest2(:)];

% for len = 1:length(k)
    
for iHt = 1: length(xTest) 

    Dist0 = Calculate_Norm(xTest(iHt,:), H0);
   
    Dist0 = sort(Dist0);
    
    Dist1 = Calculate_Norm(xTest(iHt,:), H1);
   
   Dist1 = sort(Dist1);
   DS(iHt) = log(n0/n1) + D*(log(Dist0(k)) - log(Dist1(k)));
end

DS= reshape(DS,length(x2),length(x1));

imagesc(x1([1 end]),x2([1 end]),DS);

 hold on 
 plot(H0(:,1),H0(:,2),'b*'); % H0

 plot(H1(:,1),H1(:,2),'ro') % H1
 
end