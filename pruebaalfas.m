%%%betas
alfas=[];
celprome=[];

for jjjj=5:5:65
     jjjj
    
load([num2str(['alfa' num2str(jjjj)])]);

alfas=[alfas; jjjj max(monimeris(:,5))];
celprome=[celprome; jjjj sum(xoric(2:end))/(prr)];
end

for i=1:length(VXT(1,1:end))
plot(sqrt(VXT(1:nc,i).^2+VYT(1:nc,i).^2))
axis([0 length(VXT(1,1:end)) 0 .3]);
axis square
pause(0.1)
end