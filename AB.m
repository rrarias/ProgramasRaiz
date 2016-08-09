%==========Programa para encontrar los sitios de frontera
celint=[];

%Escojo los que tienen contacto con la frontera
for q=1:length(tri)
    if any(tri(q,:)>nc)
        AA=find(tri(q,:)<nc+1);
        celint=[celint;tri(q,AA)'];
    end
end

JJ=unique(celint);
xb=x(JJ);
yb=y(JJ);

% [vv,cc]=voronoin([x,y]); 
% 
% for i=1:length(JJ)
% patch(vv(cc{JJ(i)},1),vv(cc{JJ(i)},2),i); % use color i.
% end
% 
