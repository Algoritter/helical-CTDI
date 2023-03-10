%% WRITTEN BY ALEXANDER KRINGS (2023)

clear, clc, close all hidden
load("Messwerte.mat");
format shortEng
tic
% Samples = 1000;
% VALUE PAIRS [time,value]
mz = size(Messwerte,2)/2;
% Ausnahme, Messwerte umrechnen zu mGy
% Messwerte einmalig für 2. Spalte umrechnen
Messwerte{:,2} = Messwerte{:,2}/1000;

%% ====================================================
% TABLE TO 3D-ARRAY
X = reshape(table2array(Messwerte),[],2,mz);
% CONVERSION INTO mGy
X(:,2,:) = X(:,2,1:mz)/1000;


%% PLOT AND EVAL DATA FOR A SINGLE MEASUREMENT

for m = 8 %m = 1:mz
    PlotDose(X,m); hold on
    [P.pks,P.locs,P.widths,P.proms] = FindThePeaks(X,m,10);
    DrawPeakMarkers(X,P,m)
    for n = 1:numel(P.locs)
        DrawArea(X,P,m,n)
        P.CTDI(n,1,m) = CalcDose(X,P,m,n);
        WriteDose(X,P,m,n)
    end
    hold off
    SetXLimit(X,P,m)
end


%% Find Distance of Peaks
% locs = FindPeakLocs(X,10,5);
% Dist = diff(locs);


%% PLOT THE MOVING SUMS OF DOSES
% PlotDoseSums(X,M1,M2,T)

% L = {'Pos.A';'Pos.B';'Pos.C';'Pos.D';'Pos.E';...
%     'Pos.A';'Pos.B';'Pos.C';'Pos.D';'Pos.E'};
% 
% PlotDoseSums(X,1,5,"Abdomen")
% legend(L{1:5},'Location', 'northeastoutside')

% PlotDoseSums(X,6,10,"Kopf")
% legend(L{6:10},'Location', 'northeastoutside')


%% POLYNOMIAL FIT TO SUM-CURVES
% PlotPolyCurves(X,6,10)
% PlotPolyCurves(X,3,5)
% PlotPolyCurves(X,2,2)
% PlotPolyCurves(X,1,1)

%% WRITE DELTAS INTO EXCEL
% T = DeltaTable(X,10,4000);

%% MOVING SUM
% MoveSum(X,M)
% W = MoveSum(X,1,5);

%% ExcelDoseTable(X,M,N)
% CTDI100 = DoseTable(X,10,10);
% Safe2Excel(CTDI100);

%% EXCEL CHECK
% format longEng
% CheckSumAgainstExcel(X,m,zeitidx,samples)
% CheckSumAgainstExcel(X,6,4700,1000)
% [T] = CalcTotalDose(X);

%% CLEANUP
clear mz m n 
toc

%% ===================== FUNCTION DEFINITIONS =============================

function [pks,locs,widths,proms] = FindThePeaks(X,m,N)
[~,~,widths,~] = findpeaks(X(:,2,m));
wlimit = mean(widths)*0.7;
[pks,locs,widths,proms] = findpeaks(X(:,2,m),...
    'SortStr','descend','MinPeakWidth',wlimit,...
    'MinPeakDistance',950,'NPeaks',N);
end


function [locs] = FindPeakLocs(X,N,m)
[~,~,widths,~] = findpeaks(X(:,2,m));
wlimit = mean(widths)*0.7;
[~,locs] = findpeaks(X(:,2,m),...
    'SortStr','none','MinPeakWidth',wlimit,...
    'MinPeakDistance',950,'NPeaks',N);
end


function DrawPeakMarkers(X,P,m,~)
P.raise = mean(P.pks)*.03;
plot(X(P.locs,1,m), P.pks + P.raise,...
    'MarkerFaceColor',[0 0.4470 0.7411],...
    'MarkerEdgeColor','none',...
    'Marker','v',...
    'MarkerSize',8,...
    'LineStyle','none');
text(X(P.locs,1,m)+100, P.pks + P.raise,...
    num2str((1:numel(P.pks))'))
legend(['Messung: ',num2str(m)],'Peaks')
end


function DrawArea(X,P,m,n)
x1 = P.locs(n)-499;
x2 = P.locs(n)+500;
area(X(x1:x2,1,m),X(x1:x2,2,m),EdgeColor="none" );
end


function [S] = CalcDose(X,P,m,n)
x1 = P.locs(n)-499;
x2 = P.locs(n)+500;
S = sum(X(x1:x2,2,m));
end


function WriteDose(X,P,m,n)
val1 = P.CTDI(n,1,m);
val2 = round(val1,4,"decimals");
str = [num2str(val2),' mGy'];
text(X(P.locs(n),1,m), P.pks(n)/2,str,...
    'Rotation',-40,'HorizontalAlignment','center',...
    'FontSize',14,'FontWeight','bold')
end


function SetXLimit(X,P,m)
x1 = P.locs(end-1:end);
x2 = sort(x1,1,'ascend');
x3 = X(x2,1,m);
x4 = round(x3,-2,"decimals");
x5 = [x4(1)-1000,x4(2)+1000];
xlim(x5)
end


function [W] = MoveSum(X,M1,M2)
var1 = size(X);
W = nan(var1(1),var1(3));
for m = M1:M2
    var2 = sum(not(isnan(X(:,1,m))));
    x1 = 499;
    x2 = 500;
    for x = x2:var2-x2
        W(x,m) = sum(X(x-x1:x+x2,2,m));
    end
end
end


function PlotDose(X,m)
figure("Name","Neu",Position=[100,100,1300,600])
plot(X(:,1,m), X(:,2,m),LineWidth=1,Color='#000')
title(['CTDI Messung Nr. ',num2str(m)])
xlabel('Zeit (ms)')
ylabel('Dosisrate (mGy/s)')
legend Location northeastoutside
end


function PlotDoseSums(X,M1,M2,T)
W = MoveSum(X,M1,M2);
figure("Name",T, Position=[100,100,1300,400])
title(['CTDI (flow) Messungen: ',num2str(M1),' bis ',num2str(M2)])
xlabel('Zeit (ms)')
ylabel('CTDI_{100} (mGy)')
plot(X(:,1),W(:,M1:M2),LineWidth=1)
% CROP PLOT TO ROI
var1 = max(W,[],1,"omitnan");
var2 = sort(var1,'ascend');
var3 = var2(1)*0.90;
var4 = max(var1)*1.02;
ylim([var3,var4])
end


function PlotPolyCurves(X,M1,M2)
figure(Position=[100,100,1300,600])
title(['CTDI , Messung: ',num2str(M1),' bis ',num2str(M2)])
subtitle('Polynomfit mit 95% Konfidenzintervall und relativem Fehler')
xlabel('Zeit (ms)')
ylabel('CTDI_{100} (mGy)')
W = MoveSum(X,M1,M2);
PolyDraw(W,M1,M2,4000)
end


function PolyDraw(W,M1,M2,Z)
[M,I] = max(W,[],1,"omitnan");
hold on
for m = M1:M2
    x1 = I(m) -Z/2;
    x2 = I(m) +Z/2;
    x = x1:x2;
    [P,s] = polyfit(x,W(x,m),2);
    [fit,delta] = polyval(P,x,s);
    S.x = -P(2)/2/P(1);
    S.y = P(3)-P(2)*P(2)/4/P(1);
        set(gca,'ColorOrderIndex',m)
    plot(x,W(x,m),LineWidth=1)
        set(gca,'ColorOrderIndex',m)
    plot(x,fit,LineWidth=1)
    plot(x,fit+2*delta,'m--',x,fit-2*delta,'m--')
        set(gca,'ColorOrderIndex',m)
    plot(I(m),M(m),'o')
        set(gca,'ColorOrderIndex',m)
    plot(S.x,S.y,'|',MarkerSize=20)
    % CALCULATE RELATIVE ERROR FOR 95% CONFIDENCE
    val1 = mean(delta);
    val2 = 2*val1/S.y*100;
    val3 = round(val2,4,"decimals");
    str = ['± ',num2str(val3),' %'];
    text(S.x,S.y-S.y*0.007,str,FontSize=14,...
        HorizontalAlignment='center',FontWeight='bold')
end
hold off
xlim([x1,x2])
end


function [val1,val2,val3] = CalcDelta(X,m,Z)
W = MoveSum(X,m,m);
[~,I] = max(W);
x1 = I(m) -Z/2;
x2 = I(m) +Z/2;
x = x1:x2;
[P,s] = polyfit(x,W(x,m),2);
[~,delta] = polyval(P,x,s);
S.y = P(3)-P(2)*P(2)/4/P(1);
val1 = mean(delta);
val2 = S.y;
val3 = 2*val1/S.y;
end


function [T] = DeltaTable(X,M2,Z)
T = nan(M2,5);
for m = 1:M2
    [val1,val2,val3] = CalcDelta(X,m,Z);
    T(m,1) = val1;      %Delta
    T(m,2) = val1*2;
    T(m,3) = val2;      %Vertex
    T(m,4) = val3;      %Relative Error
    T(m,5) = val3*100;
end
header = {'Delta','2 Delta','Vertex',...
    'Relative Error','Relative Error in %'};
T = array2table(T,"VariableNames",header);
filename = 'MATLAB.Deltas.xlsx';
writetable(T,filename,'Sheet',1)
end


function [CTDI] = DoseTable(X,M,N)
CTDI = nan(N,M);
for m = 1:M
    [~,locs,~,~] = FindThePeaks(X,m,N);
    for n = 1:numel(locs)
        x1 = locs(n)-499;
        x2 = locs(n)+500;
        CTDI(n,m) = sum(X(x1:x2,2,m));
    end
end
CTDI = sort(CTDI,1,"descend","MissingPlacement","last");
CTDI = array2table(CTDI);
end


function Safe2Excel(CTDI100)
filename = 'MATLAB.Results.xlsx';
writetable(CTDI100,filename,'Sheet',1)
end


function [S] = CheckSumAgainstExcel(X,m,zeitidx,samples)
x1 = find(X(:,1,m)==zeitidx);
S = sum(X(x1:x1+samples-1,2,m));
end


function [T] = CalcTotalDose(X)
T = squeeze(sum(X(1:end,2,:),1,"omitnan"));
filename = 'MATLAB.Results.xlsx';
writematrix(T,filename,'Sheet',4)
end
