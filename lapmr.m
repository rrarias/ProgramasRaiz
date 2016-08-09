function[lapm]=lapm(x,y,nvert,tria,neib)


Adjc=zeros(nvert);
Degree=Adjc;

for j=1:nvert   
% neibr=[];
% for q=1:length(tria)
%     for jj=1:3
%         if x(j)==x(tria(q,jj))
%             neibr=[neibr;tria(q,jj+1)];
%         end
%     end
% end
% neibr=unique(neibr);
neibr=neib{j}';
neibr=neibr(1:find(neibr<=nvert, 1, 'last' ));

for nn=1:length(neibr)

Adjc(j,neibr(nn)) = 1/distance([x(j),y(j)],[x(neibr(nn)),y(neibr(nn))]);
Adjc(neibr(nn),j)=Adjc(j,neibr(nn));   
end
end


for diag=1:nvert
    Adjc(diag,diag)=0;
    Degree(diag,diag)=sum(Adjc(diag,:));
end

lapm=Adjc-Degree;

for diag=1:nvert
    
        
    lapm(diag,:)=(1/Degree(diag,diag))*lapm(diag,:);
    
end
    
    