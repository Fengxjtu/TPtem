function [ r_et,r_max,p_rmax,t_rmax,obs_rmax,yr_sp ] = ...
    corr_analysis( chr,chr_yr,obs,obs_yr,time_scale,time_span,var_name,IPLT )
%[ r_et,r_max,p_rmax,t_rmax,obs_rmax,yr_sp ]=corr_analysis( chr,chr_yr,obs,obs_yr,time_scale,time_span,var_name,IPLT )
%   chr: the chronology with form of column vector
%   chr_yr: the year span of chronology
%   obs: the observation data for correlation analysis (i.e. temperature,
%        precipitation, PDSI and so on) (NOTICE: (1) If you want to plot
%        the more than one variables on one figure, you should pooled the 
%        obs matrices into three dimensional like {obs_new(:,:,1)=obs_temp;
%        obs_new(:,:,2)=obs_prec;} before inputting data, which third di-
%        mension denotes the different variables. (2) Before inputting 
%        data, you must average the two-day's data of 28-29th Feb. into 
%        one day in leap years.)
%   obs_yr: the year span of observation data
%   time_scale: time scale, a string only can be set as 'm', 'd', 'h', 'x'.
%   time_span: the time span, usually from the last year to current year
%        (The input must be a string and formed as P(C)XXXC(P)XXX, here 
%        the XXX is the number of time at each scale and the P denotes the 
%        previous year, C for current. ADVISE: Set the time_span as long 
%        as possible in the interval from previous years to current years.)
%   var_name: variable names (Usually:
%           tmean   平均温度
%           tmax    最大温度
%           tmin    最低温度
%           pr      降水
%           rhmean  平均相对湿度
%           PDSI    干旱指数 )
%   IPLT: whether plot(1 for do, 0 for not)
%   r_et: correlation coefficient(r) of each time combination
%   r_max: the maximum r for all the possible time intervals
%   p_rmax: the p-value of maximum r
%   t_rmax: the time intervals of maximum r formed as 'P(C)X(XX)-P(C)X(XX)'
%   obs_rmax: the combinational sequence of observed data with maximum r
%   yr_sp: the actual year span for correlation analysis
% Some constraints of input datum:
%    TIME_      SIZE OF OB-     FORM OF TIME SPAN        EXPLANATION
%    SCALE      SERVED DATA        (EXAMPLES)
%  m(month月)    year×12       PXXCXX  (P8C10)    上年8月至当年10月
%                year×365                   idem (同上)
%  d(day日)      year×365      PXXXCXXX(P315C80)  上年第315日至当年第80日
%  h(penkad侯)   year×73       PXXCXX  (P32C50)    上年第32侯至当年第50侯
%                year×365                   idem (同上)
%  x(dekad旬)    year×36       PXXLXX  (P10L30)    上年第10旬至上年第30旬
%                year×365                   idem (同上)
%
%  CODE WRITER: Fang Congxi
%  CONTACT: sfh_st_cn2@163.com
%  ADRESS: Tree-ring Laboratory, Insititute of Earth Environment CAS, 
%          Xi'an, China (postcode: 710061).

if nargin<8
    IPLT = 1;
end
if nargin<7
    for i = 1:size(obs,3)
        var_name{i} = ['VAR',num2str(i)];
    end
end
if nargin<6
    time_span = [];
end
if nargin<5
    time_scale = [];
end
if nargin<4
    error('No enought input');
end

if isempty(time_scale)
    time_scale = 'M';
end
time_scale = upper(time_scale);
if isempty(time_span)
    switch time_scale
        case 'M'
            time_span = 'P1C10';
        case 'D'
            time_span = 'P1C304';
        case 'H'
            time_span = 'P1C61';
        case 'X'
            time_span = 'P1C30';
    end
end
time_span = upper(time_span);

MN = ['Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'];
MD = [31; 28; 31; 30; 31; 30; 31; 31; 30; 31; 30; 31 ];

[mc,~] = size(chr);
if mc~=1
    chr = chr';
end
lc = length(chr);
[nyo,nto,V] = size(obs);
if ~(nto==12 || nto==365 || nto==73 || nto==36)
    error(['The column of observed data(obs) must be 12 or 365 or 73 ', ...
        'or 36, here it is: ',num2str(nto)])
end
if lc~=length(chr_yr)
    error(['The year length of chr_yr data must be same as that of ' ...
        'chronology, here L_chr-yr=',num2str(length(chr_yr)),', L_chr=', ...
        num2str(lc)])
end
if nyo~=length(obs_yr)
    error(['The year length of chr_yr data must be same as that of ' ...
        'chronology, here L_obs-yr=',num2str(length(obs_yr)),', NY_obs=', ...
        num2str(nyo)])
end

switch time_scale
    case 'M'
        if ~(nto==12 || nto==365)
            error(['The column of observed data(obs) must be 12 or ', ... 
                '365 for monthly time scale, here it is: ',num2str(nto)])
        end
        t_tot = 12;
        obs_adj = zeros(nyo,t_tot,V);
        if nto==365
            i = 0;
            for j=1:t_tot
                obs_adj(:,j,:) = mean(obs(:,i+1:i+MD(j),:),2,'omitnan');
                i=i+MD(j);
            end
        elseif nto==12
            obs_adj = obs;
        end
        
    case 'D'
        if ~(nto==365)
            error(['The column of observed data(obs) must be 365 ', ...
                'for daily time scale, here it is: ',num2str(nto)])
        end
        t_tot = 365;
        obs_adj = obs;
        
    case 'H'
        if ~(nto==365 || nto==73)
            error(['The column of observed data(obs) must be 73 ', ...
                'or 365 for pentad time scale, here it is: ',num2str(nto)])
        end
        t_tot = 73;
        obs_adj = zeros(nyo,t_tot,V);
        if nto==365
            for i=1:t_tot
                obs_adj(:,i,:) = mean(obs(:,(i-1)*5+1:i*5,:),2,'omitnan');
            end
        elseif nto==73
            obs_adj = obs;
        end
        
    case 'X'
        if ~(nto==365 || nto==36)
            error(['The column of observed data(obs) must be 36 or ' ...
                '365 for pentad time scale, here it is: ',num2str(nto)])
        end
        t_tot = 36;
        obs_adj = zeros(nyo,t_tot,V);
        if nto==365
            i = 0;
            for j=1:length(MD)
                obs_adj(:,(j-1)*3+1,:) = mean(obs(:,i+1:i+10,:),2,'omitnan');
                obs_adj(:,(j-1)*3+2,:) = mean(obs(:,i+11:i+20,:),2,'omitnan');
                obs_adj(:,(j-1)*3+3,:) = mean(obs(:,i+21:i+MD(j),:),2,'omitnan');
                i=i+MD(j);
            end
        elseif nto==36
            obs_adj = obs;
        end
        
    otherwise
        error(['This is not a valid time_scale: ',time_scale])
end

%  Get time span.
IP = find(time_span=='P');
IC = find(time_span=='C');
t_staf = 'P';
t_endf = 'C';
if isempty(IP)
    t_staf = 'C';
    t_sta = str2double(time_span(IC(1)+1:IC(2)-1));
    t_end = str2double(time_span(IC(2)+1:end));
elseif isempty(IC)
    t_endf = 'P';
    t_sta = str2double(time_span(IP(1)+1:IP(2)-1));
    t_end = str2double(time_span(IP(2)+1:end));
else
    t_sta = str2double(time_span(IP+1:IC-1));
    t_end = str2double(time_span(IC+1:end));
end

%  Adjust the observation data.
if strcmp(t_staf,'P') && strcmp(t_endf,'C')
    t_num(1:t_tot-t_sta+1) = t_sta+(0:t_tot-t_sta);
    t_flg(1:t_tot-t_sta+1) = {t_staf};
    t_sgn(1:t_tot-t_sta+1) = 0;
    t_num(t_tot-t_sta+2:t_tot-t_sta+1+t_end) = 1:t_end;
    t_flg(t_tot-t_sta+2:t_tot-t_sta+1+t_end) = {t_endf};
    t_sgn(t_tot-t_sta+2:t_tot-t_sta+1+t_end) = 1;
    obs_adj2(2:nyo,1:t_tot-t_sta+1,:) = obs_adj(1:nyo-1,t_sta+(0:t_tot-t_sta),:);
    obs_adj2(:,t_tot-t_sta+2:t_tot-t_sta+1+t_end,:) = obs_adj(1:nyo,1:t_end,:);
    obs_adj2(1,1:t_tot-t_sta+1,:) = nan;
    obs_yr_adj = obs_yr;
elseif strcmp(t_staf,t_endf) && (strcmp(t_endf,'P')|| strcmp(t_endf,'C'))
    t_num(1:t_end-t_sta+1) = t_sta+(0:t_end-t_sta);
    t_flg(1:t_end-t_sta+1) = {t_staf};
    if strcmp(t_endf,'P')
        t_sgn(1:t_end-t_sta+1) = 0;
        obs_adj2(:,1:t_end-t_sta+1,:) = obs_adj(1:nyo,t_sta+(0:t_end-t_sta),:);
        obs_yr_adj = obs_yr+1;
    else
        t_sgn(1:t_end-t_sta+1) = 1;
        obs_adj2(:,1:t_end-t_sta+1,:) = obs_adj(1:nyo,t_sta+(0:t_end-t_sta),:);
        obs_yr_adj = obs_yr;
    end
else
    error('That is irrational!')
end
iy0 = find(chr_yr==max(chr_yr(1),obs_yr_adj(1)));
iy1 = find(chr_yr==min(chr_yr(end),obs_yr_adj(end)));
chr_c = chr(iy0:iy1)';
iy0 = find(obs_yr_adj==max(chr_yr(1),obs_yr_adj(1)));
iy1 = find(obs_yr_adj==min(chr_yr(end),obs_yr_adj(end)));
obs_adj2 = obs_adj2(iy0:iy1,:,:);

%  Give the initial value for outputs.
r_et = zeros(length(t_num),length(t_num),V)*nan;
r_max = zeros(1,V)*nan;
p_rmax = zeros(1,V)*nan;
t_rmax(1,1:V) = {[]};
obs_rmax = [];

%  Deal for NANs and adjust the chronology and observation data into
%  uniform.
IQ = squeeze(min(~isnan(obs_adj2),[],2));
kv = 0;
IV = [];
for i = 1:size(IQ,2)
    if sum(IQ(:,i))>=20
        kv = kv+1;
        IV(kv) = i;
    end
end
if ~isempty(IV)
    yr_sp = obs_yr_adj(iy0:iy1);
    obs_adj3 = obs_adj2(:,:,IV);
else
    yr_sp = [];
    return
end
%IQC = min(min(~isnan(obs_adj3),[],2),[],3);
%obs_adj4 = obs_adj3(IQC,:,:);
%obs_yr_sp = obs_yr_adj(IQC);
%ky = 0;
%if length(obs_yr_sp)>=20
%    for i = 1:length(obs_yr_sp)
%        IY = find(chr_yr==obs_yr_sp(i));
%        if ~isempty(IY)
%            ky = ky+1;
%            yr_sp(ky) = chr_yr(IY);
%            chr_adj(ky) = chr(IY);
%            obs_adj5(ky,:,:) = obs_adj4(i,:,:);
%        end
%    end
%else
%    yr_sp = [];
%    return
%end
%obs_adj5 = squeeze(obs_adj5);

%  Calculate the critical R for correlative test on certain significant
%  level.
m = nyo-1;
n = length(t_num);
[ r_sg ] = rinv(.05,m-2);

%  The core part for combination correlation coefficient calculation.
for k=1:kv
    [r_et0(:,:,k),r_max0(k),p_rmax0(k),loc_rmax,obs_rmax0(:,k)] = ...
        SFR( chr_c,squeeze(obs_adj3(:,:,k)),time_scale );
    if ~isempty(loc_rmax)
        if loc_rmax(2)==loc_rmax(1)
            t_rmax0(k) = {[char(t_flg(loc_rmax(1))),num2str(t_num(loc_rmax(1)))]};
        else
            t_rmax0(k) = {[char(t_flg(loc_rmax(1))),num2str(t_num(loc_rmax(1))),'-',...
                char(t_flg(loc_rmax(2))),num2str(t_num(loc_rmax(2)))]};
        end
    else
        t_rmax0(k) = {'none'};
    end
    if ~(strcmpi(var_name{k},'PR')||strcmpi(var_name{k},'EVPL')||strcmpi(var_name{k},'EVPS'))
        obs_rmax0(:,k) = obs_rmax0(:,k)/(loc_rmax(2)-loc_rmax(1)+1);
    end
    r_eb(:,k) = diag(r_et0(:,:,k));
end
r_et(:,:,IV) = r_et0;
r_max(IV) = r_max0;
p_rmax(IV) = p_rmax0;
t_rmax(IV) = t_rmax0;
obs_rmax = zeros(size(obs_rmax0,1),V)*nan;
obs_rmax(:,IV) = obs_rmax0;
V = kv;

if IPLT==0
    return
end

switch time_scale
    case 'M'
        mltv = 1/2;
        mltp = 1;
        ts = 'month';
    case 'D'
        mltv = 5;
        mltp = 50;
        ts = 'day';
    case 'H'
        mltv = 3/2;
        mltp = 10;
        ts = 'pentad';
    case 'X'
        mltv = 1;
        mltp = 5;
        ts = 'dekad';
end
X = 1:n;
XV = n+mltv+mltp/2;
XC = [X,XV];
xlmin = min(XC)-mltv;
xlmax = max(XC)+mltp/2;
r_c = [r_eb;r_max0];
for k = 1:V
    bscolor(k,:) = [(k-1)/V (k-1)/V 1-(k-1)/V];
end
bscolor(1,:) = [0 0 0];
vrn_lgd = tool_legend_varname(var_name(IV));
if strcmp(time_scale,'M') && V<=5
    hb = bar(XC,r_c,'DisplayName','X','barwidth',.68);
    for k=1:V
        hb(k).FaceColor = bscolor(k,:);
    end
    hl = legend(hb,vrn_lgd,'location','northwest','Position',[.1 .83 .69 .05],'Orientation','Horizontal','box','off');
    hold on
else
    for k=1:V
        hp(k) = plot(X,r_eb(:,k),'Color',bscolor(k,:));
        hold on
    end
    hl = legend(hp,vrn_lgd,'location','northwest','Orientation','Horizontal','box','off');
    r_v = zeros(2,V);
    r_v(1,:) = r_max0;
    hb = bar([XV,XV+mltp],r_v,'DisplayName','X','barwidth',.68);
    for k=1:V
        hb(k).FaceColor = bscolor(k,:);
    end
    hold on
    plot([xlmin xlmax],[0 0],'Color','k')
    hold on
end
plot([xlmin xlmax],[r_sg r_sg],'k--');
hold on
plot([xlmin xlmax],[-r_sg -r_sg],'k--');
hold on
yl = get(gca,'ylim');

% Set ticks
if strcmp(time_scale,'M')
    xt = XC;
    for i=1:n
        xtl(i) = {MN(t_num(i),1)};
    end
    xtl(n+1) = {[]};
else
    if strcmp(time_scale,'D')
        nump = 15:50:365;
        numc = 50:50:365;
    elseif strcmp(time_scale,'H')
        nump = 3:10:73;
        numc = 10:10:73;
    elseif strcmp(time_scale,'X')
        nump = 1:5:36;
        numc = 5:5:36;
    end
    i = 0;
    for num = nump
        I = find((t_num==num)+(~t_sgn)==2);
        if ~isempty(I)
            i = i+1;
            xt(i) = XC(I);
            xtl(i) = {num2str(num)};
        end
    end
    for num = numc
        I = find((t_num==num)+(t_sgn)==2);
        if ~isempty(I)
            i = i+1;
            xt(i) = XC(I);
            xtl(i) = {num2str(num)};
        end
    end
    xt(i+1) = XV;
    xtl(i+1) = {[]};
    set(gca,'XMinorTick','on');
end
P = [.09 .2 .835 .7];
axis([xlmin xlmax yl(1) yl(2)+.1]);
set(gca,'position',P,'xtick',xt,'xticklabel',xtl,'Fontname','Times New Roman');
yl = get(gca,'ylim');

% Add some necessary elements.
for k=1:V
    text(double(XV-mltp/5*2),double(yl(1)-(yl(2)-yl(1))*.04*k),t_rmax0(k),'HorizontalAlignment','left','fontsize',8,'Fontname','Times New Roman');
end
if strcmp(t_staf,'P') && strcmp(t_endf,'C')
    I_seg = (find(~t_sgn, 1, 'last' )+find(t_sgn, 1, 'first' ))/2;
    I_end = XC(n)+mltv;
    seg_pos = P(1)+P(3)/(xlmax-xlmin)*(I_seg-xlmin);
    end_pos = P(1)+P(3)/(xlmax-xlmin)*(I_end-xlmin);
    annotation('doublearrow',[P(1)+.005 seg_pos-.005],[P(2)-.07 P(2)-.07],'HeadSize',6);
    annotation('doublearrow',[seg_pos+.005 end_pos-.005],[P(2)-.07 P(2)-.07],'HeadSize',6);
    annotation('line',[P(1) P(1)],[P(2)-.09 P(2)-.05]);
    annotation('line',[seg_pos seg_pos],[P(2)-.09 P(2)-.05]);
    annotation('line',[end_pos end_pos],[P(2)-.09 P(2)-.05]);
    text((xlmin+I_seg)/2,min(yl)-(max(yl)-min(yl))*.13,'Previous year','HorizontalAlignment','center','Fontname','Times New Roman');
    text((I_seg+I_end)/2,min(yl)-(max(yl)-min(yl))*.13,'Current year','HorizontalAlignment','center','Fontname','Times New Roman');
elseif strcmp(t_staf,'P') && strcmp(t_endf,'P')
    I_end = XC(n)+mltv;
    end_pos = P(1)+P(3)/(xlmax-xlmin)*(I_end-xlmin);
    annotation('doublearrow',[P(1)+.005 end_pos-.005],[P(2)-.07 P(2)-.07],'HeadSize',6);
    annotation('line',[P(1) P(1)],[P(2)-.09 P(2)-.05]);
    annotation('line',[end_pos end_pos],[P(2)-.09 P(2)-.05]);
    text((xlmin+I_end)/2,min(yl)-(max(yl)-min(yl))*.13,'Previous year','HorizontalAlignment','center','Fontname','Times New Roman');
elseif strcmp(t_staf,'C') && strcmp(t_endf,'C')
    I_end = XC(n)+mltv;
    end_pos = P(1)+P(3)/(xlmax-xlmin)*(I_end-xlmin);
    annotation('doublearrow',[P(1)+.005 end_pos-.005],[P(2)-.07 P(2)-.07],'HeadSize',6);
    annotation('line',[P(1) P(1)],[P(2)-.09 P(2)-.05]);
    annotation('line',[end_pos end_pos],[P(2)-.09 P(2)-.05]);
    text((xlmin+I_end)/2,min(yl)-(max(yl)-min(yl))*.13,'Current year','HorizontalAlignment','center','Fontname','Times New Roman');
end
xlabel(ts,'position',[(xlmin+I_end)/2 min(yl)-(max(yl)-min(yl))*.19]);
ylabel('correlation coefficient');
text(xlmax+.01*(xlmax-xlmin),r_sg,'p=0.05','Fontname','Times New Roman')
text(xlmax+.01*(xlmax-xlmin),-r_sg,'p=0.05','Fontname','Times New Roman')
plot([XC(n)+mltv XC(n)+mltv],[min(yl) max(yl)],'color','k')
end

function [r_et,r_max,p_rmax,loc_rmax,obs_rmax] = SFR( CHR,OBS,TS )
% The core subfunction of corr_analysis.

switch TS
    case 'M'
        maxj = 12;
    case 'D'
        maxj = 365;
    case 'H'
        maxj = 73;
    case 'X'
        maxj = 36;
end

[M,N] = size(OBS);
r_max = 0;
p_rmax = 1;
r_et = zeros(N,N)*nan;
loc_rmax = [];
obs_rmax = zeros(M,1)*nan;
for i=1:N
    o_to = zeros(M,1);
    for j = i:N
        o_to = o_to+OBS(:,j);
        chr_c = CHR( ~isnan(CHR) & ~isnan(o_to) );
        oto_c = o_to( ~isnan(CHR) & ~isnan(o_to) );
        [r_co,p_co] = corrcoef(chr_c,oto_c);
        r_et(j,i) = r_co(2);
        if abs(r_co(2))>abs(r_max)
            r_max = r_co(2);
            p_rmax = p_co(2);
            loc_rmax = [i j];
            obs_rmax = o_to;
        end
        if j-i+1>=maxj
            break
        end
    end
end
end