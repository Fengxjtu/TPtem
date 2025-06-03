
%clear
%load FASM.mat;
L = 230;
alpha = 0.01;
for i = 1:length(cj)-L+1
    [b,bint] = regress(cj(i:i+L-1),[Year(i:i+L-1) ones(size(Year(i:i+L-1)))],alpha);
    tr1(i) = b(1);
    tr2(i) = b(2);
    tr1int(1,i) = bint(1,1);
    tr1int(2,i) = bint(1,2);
end
I=find(tr1==min(tr1));
minY = Year(I+L-1);
mintr = tr1(I);
fe = Year(I:I+L-1)*tr1(I)+tr2(I);

figure('position',[100 100 560 300])
bsvl = min(tr1int(1,:))-.01;
ha(2) = area(Year(L:end),tr1int(2,:),'linestyle','none','basevalue',bsvl,'ShowBaseLine','off','FaceColor',[.5 .5 .5]);
hold on
area(Year(L:end),tr1int(1,:),'linestyle','none','basevalue',bsvl,'ShowBaseLine','off','FaceColor',[1 1 1]);
hold on
ha(1) = plot(Year(L:end),tr1,'linewidth',1.5);
hold on
legend(ha,{'slope evaluation';'99% confidence'})
xl = [min(Year(L:end))-1 max(Year(L:end))+1];
%yl = [-1.4 1.8];
yl = [-0.004 0.006];
%text(xl(1)+(xl(2)-xl(1))*.035,yl(2)-(yl(2)-yl(1))*.08,['min slope: ',num2str(roundn(mintr,-3))]);
%text(xl(1)+(xl(2)-xl(1))*.035,yl(2)-(yl(2)-yl(1))*.16,['min slope year: ',num2str(minY-L+1),'-',num2str(minY)]);
set(gca,'xlim',xl,'ylim',yl,'Fontsize',30)
title(['L=',num2str(L)],'Fontsize',30);
xlabel('Year','Fontsize',30);
ylabel('Slope','Fontsize',30);
set(gcf, 'PaperPositionMode', 'auto');
print(gcf,'-djpeg',['Ç÷ÊÆ',num2str(L),'Äê'],'-r600');