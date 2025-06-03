clear
recon=xlsread('d:\mat131\FenDuanjianyan\Book2.xls','sheet2');
  cal1=1961;cal2=1990;
  ver1=1991;ver2=2010;
  %cal1=1981;cal2=2010;
  %ver1=1961;ver2=1980;
  d11=find(recon(:,1)==cal1);d12=find(recon(:,1)==cal2);
  d13=find(recon(:,1)==ver1);d14=find(recon(:,1)==ver2);
  d21=d11;d22=d12;d23=d13;d24=d14;
  
  x1=[ones(d12-d11+1,1) recon(d11:d12,4)];
  y1=recon(d21:d22,2);
  a=regress(y1,x1);  
  y=[ones(d14-d13+1,1) recon(d13:d14,4)]*a;
  y11=[ones(d12-d11+1,1) recon(d11:d12,4)]*a;
  %y12=[ones(1996-1960+1,1) recon(d13:d12,4)]*a;
  y2=recon(d23:d24,2);
  yc=mean(y1);
  yv=mean(y2);
  
%signtest
%%[p3,h3_1,r_signtest]=signtest(y,y2,0.01);
y111=y1-yc;
y112=y11-yc;
rst=0;
for i=1:cal2-cal1+1
    if (y111(i)>=0)&(y112(i)>=0)
        rst=rst+1;
    else if (y111(i)<0)&(y112(i)<0)
            rst=rst+1;
        end
    end
end
    y221=y2-yv;
y222=y-yv;
vst=0;
for i=1:ver2-ver1+1
    if (y221(i)>=0)&(y222(i)>=0)
        vst=vst+1;
    else if (y221(i)<0)&(y222(i)<0)
            vst=vst+1;
        end
    end
    end
%r
%[r_all,p_all]=corrcoef(y12,recon(d23:d22,2));
[r_ver,p_ver]=corrcoef(y,y2);
[r_cal,p_cal]=corrcoef(y1,y11);
%R^2
SS_R=0;SS_T=0;
for i=1:cal2-cal1+1
    SS_R=SS_R+(y11(i)-y1(i))^2;
    SS_T=SS_T+(y1(i)-yc)^2;
end
R2=1-SS_R/SS_T;
%RE
SSR=0;SSM=0;
for i=1:ver2-ver1+1
    SSR1=(y(i,1)-y2(i,1))^2;
    SSM1=(y2(i,1)-yc)^2;
    SSR=SSR+SSR1;
    SSM=SSM+SSM1;
end
RE=1-SSR/SSM;
%CE
SSR=0;SSv=0;
for i=1:ver2-ver1+1
    SSR1=(y(i,1)-y2(i,1))^2;
    SSv1=(y2(i,1)-yv)^2;
    SSR=SSR+SSR1;
    SSv=SSv+SSv1;
end
CE=1-SSR/SSv; 

calibration(1)=r_cal(1,2);
calibration(2)=R2;
calibration(3)=rst;
ver(1)=r_ver(1,2);
ver(2)=RE;
ver(3)=CE;
ver(4)=vst;
