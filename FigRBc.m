
%%  Figura para las dinamicas vistas desde la punta de la raiz
figure(1)
clf
[aa bb]=sort(y(1:xori+length(xv2)));
em=max(Ep(bb));
hold on
plot(y(bb),(Ep(bb))/em,'b');
plot(y(bb),c(bb)/max(c(bb)),'r')
plot(y(bb),r(bb).*sqrt(c(bb))/gamma,'k')
pause(.01)
plot(y(1),0,'*r')
hold off

 [vv,cc] = voronoin([x,y]);
    xv3=x5;
    yv3=y5;
    [IN ON] = inpolygon(x,y,xv3,yv3);
for i =1:length(cc)
    if ON(i)==0;
        nuevovoro{i}=cc{i};
        nuevovoro{i}=horzcat(nuevovoro{i},nuevovoro{i}(1));
    end
end

%% figura para el potencial
figure(2)
clf
Cau=Ep;%/sum(Ep);
fEc=(Cau);

Z=nuevovoro;
for i=1:length(nuevovoro)
    tt=Z{i};
    tt(:)=1;
    Z{i}=fEc(i)*tt;
end

hold on
%voronoi(x,y,'k')
axis xy ,colormap 'jet'
axis equal, axis off, view(0,90);

for i=1:length(nuevovoro)
    line(vv(nuevovoro{i},1),vv(nuevovoro{i},2),'Color','k','LineWidth',1.1)
    patch(vv(nuevovoro{i},1),vv(nuevovoro{i},2),Z{i}',Z{i}')
    
    A=vv(nuevovoro{i},1);
    B=vv(nuevovoro{i},2);
    E=Z{i}';
    Xi=[A,B,E];
    Xf=[A,B,0*E];
    hold on
    for j=1:length(A)
        line('xdata',[Xi(j,1),Xf(j,1)],'ydata',[Xi(j,2),Xf(j,2)],'zdata',[Xi(j,3),Xf(j,3)],'linestyle',':')
    end
end
% 
%  plot3(x(1),y(1),fEc(1),'.w','markersize',15,'LineWidth',1.1)
%  plot3(x(2),y(2),fEc(2),'.w','markersize',15,'LineWidth',1)
%  plot3(x(1),y(1),fEc(1),'ok','markersize',5,'LineWidth',1.1)
%  plot3(x(2),y(2),fEc(2),'ok','markersize',5,'LineWidth',1.1)
%  plot3(x(1),y(1),fEc(1),'dw','markersize',6,'LineWidth',2.1)
%  plot3(x(2),y(2),fEc(2),'dw','markersize',6,'LineWidth',2.1)
%    plot3(x(JJ),y(JJ),fEc(JJ),'.w','markersize',18,'LineWidth',1.1)
%    plot3(x(JJ),y(JJ),fEc(JJ),'ok','markersize',7,'LineWidth',1.1)

plot3(x(JJ),y(JJ),fEc(JJ)+max(Ep),'o','markersize',7,'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1.1)
axis([-5+min(x) max(x)+5 min(y)-.2 max(y)+.2]),
%axis off,
%colorbar
hold off
Q=getframe(gca);
QN(prr)=Q;

%%   figura para la auxina
figure(3)
clf
c(c==0)=0.0001;
Cau=c;%log(c*0.01);

%Cau=(c);
fEc=Cau;

Z=nuevovoro;
for i=1:length(nuevovoro)
    tt=Z{i};
    tt(:)=1;
    Z{i}=fEc(i)*tt;
    
end

hold on
axis xy ,colormap 'jet'
axis equal, axis off, view(0,90);

for i=1:length(nuevovoro)
    line(vv(nuevovoro{i},1),vv(nuevovoro{i},2),'Color','k','LineWidth',1.1)
    patch(vv(nuevovoro{i},1),vv(nuevovoro{i},2),Z{i}',Z{i}')
    
    A=vv(nuevovoro{i},1);
    B=vv(nuevovoro{i},2);
    E=Z{i}';
    Xi=[A,B,E];
    Xf=[A,B,0*E];
    hold on
    for j=1:length(A)
        line('xdata',[Xi(j,1),Xf(j,1)],'ydata',[Xi(j,2),Xf(j,2)],'zdata',[Xi(j,3),Xf(j,3)],'linestyle',':')
    end
end

%  plot3(x(1),y(1),fEc(1),'.w','markersize',15,'LineWidth',1.1)
%  plot3(x(2),y(2),fEc(2),'.w','markersize',15,'LineWidth',1)
%  plot3(x(1),y(1),fEc(1),'ok','markersize',5,'LineWidth',1.1)
%  plot3(x(2),y(2),fEc(2),'ok','markersize',5,'LineWidth',1.1)
%  plot3(x(1),y(1),fEc(1),'dw','markersize',6,'LineWidth',2.1)
%  plot3(x(2),y(2),fEc(2),'dw','markersize',6,'LineWidth',2.1)
plot3(x(JJ),y(JJ),fEc(JJ)+1,'o','markersize',7,'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1.1)
%colorbar
hold off
Q1=getframe(gca);
QN1(prr)=Q1;

%%  figura para el flujo

  [vv,cc] = voronoin([x,y]);
    xv3=x5;
    yv3=y5;
    [IN ON] = inpolygon(x,y,xv3,yv3);
    for i =1:length(cc)
        if ON(i)==0;
            nuevovoro{i}=cc{i};
            nuevovoro{i}=horzcat(nuevovoro{i},nuevovoro{i}(1));
        end
    end
    
    fluxtot=fluxt;    % se hace promedio temporal para el flujo, depende de mo
    normt=sqrt(fluxtot(:,3).^2+fluxtot(:,4).^2);
	normt1=normt;
	normt(:)=5*normt(:).^(-1);


    figure(4)  
clf
    
anpol1=0*fluxtot(:,3);
anpol1(:)=sqrt(fluxtot(:,3).^2+((fluxtot(:,4)).^(2)));
anpol1(anpol1==0)=.001;

hold on
Z=nuevovoro;

    for i=1:length(nuevovoro)
		
		anpol(i)=(anpol1(i)/(max(anpol1)));
        tt=Z{i};
        tt(:)=1;
        Z{i}=anpol1(i)*tt;
		patch(vv(nuevovoro{i},1),vv(nuevovoro{i},2),[1 1 .85],'LineWidth',1.1)
    end

 for   j=1:xori+length(xv2)

  if fluxtot(j,4)<0

gils=quiver(x(j),y(j),normt(j)*(fluxtot(j,3)),normt(j)*(fluxtot(j,4)),.06,'Color','r','LineWidth',1.01);
addquiverarrowhead(gils,7);   
  else 

gils=quiver(x(j),y(j),normt(j)*(fluxtot(j,3)),normt(j)*(fluxtot(j,4)),.06,'Color','b','LineWidth',1.01);
addquiverarrowhead(gils,7);
  end
  end
%.1*min(anpol(j)/max(anpol))+.001   %scale for flux
 axis equal, axis off
 view(0,90)
 axis([-5+min(x) max(x)+5 min(y)-.2 max(y)+.2])
 set(gca,'xtick',[],'ytick',[])
hold off
Q2=getframe(gca);
QN2(prr)=Q2;
   
%% figura para la polarizacion
 figure(5)
  clf
  hold on
for i=1:length(nuevovoro)    
    
    line(vv(nuevovoro{i},1),vv(nuevovoro{i},2),'Color','k','LineWidth',1.1)
    tri= delaunay(x,y);
    unique (tri);
        tria=[tri tri(:,1)];
        neib=neibs(x,y,length(x),tria);      
		patch(vv(nuevovoro{i},1),vv(nuevovoro{i},2),[1 1 .85],'LineWidth',1.1)
    neibr=neib{i}';
    
for nn=1:length(neibr)
            if neibr(nn)<=length(nuevovoro)
               
                pin(i,neibr(nn))=heaviside((Ep(neibr(nn))-Ep(i))*(c(neibr(nn))-c(i)));
                if pin(i,neibr(nn))>0.2 && Ep(i)>=sum(Ep(1:xori))/(xori)% && fluxt(i)fluxt(nn)
                    line([x(i) (x(neibr(nn))+x(i))/2],[y(i) (y(neibr(nn))+y(i))/2],'Color','b','LineWidth',1.5,'LineStyle','-')
                end
                if pin(i,neibr(nn))>=0.2 && Ep(i)<sum(Ep(1:xori))/(xori) %&& fluxt(i)<fluxt(nn)
                    line([x(i) (x(neibr(nn))+x(i))/2],[y(i) (y(neibr(nn))+y(i))/2],'Color','r','LineWidth',1.5,'LineStyle','-')
                end
            end
        end
end
axis equal, axis off
axis([-5+min(x) max(x)+5 min(y)-.2 max(y)+.2])
set(gca,'xtick',[],'ytick',[])
box on
hold off
Q3=getframe(gca);
QN3(prr)=Q3;