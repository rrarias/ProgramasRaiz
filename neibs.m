function[neib]=neibs(x,y,nvert,tria)

neib=[];
for j=1:nvert   
neibs=[];
for q=1:length(tria)
    for jj=1:3
        if x(j)==x(tria(q,jj))
            neibs=[neibs;tria(q,jj+1)];
        end
    end
end
neibs=unique(neibs);
%neibs=neibs(1:find(neibs<=nvert, 1, 'last' ));
neib=[neib;{[neibs']}];
end