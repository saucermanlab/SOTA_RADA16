function f = plotCorrelation_Tenaya(dataStruct,d1_conversion,annotate)
%% FUNCTION DESCRIPTION
% Written by Matthew Sutcliffe doi: 10.1038/s41598-018-19539-z
% Adapted by Alice Luanpaisanon
% This function plots the correlation vs. pixel offset distance, as shown
% in Figure 3B and 3D.
%
% INPUTS
% dataStruct : metric structure for a single cell with correlation data
% d1_conversion : pixels to microns conversion for a single dimension
% OUTPUTS
% f : function handle

%%
data = dataStruct.CorrelationData;
f = figure('units','pixels','position',[20 -10 600 400]);
hold on;
cmap = hsv(size(data,1));
%x = 0:0.1:10;
x = 0:0.1:40; %from master plotCorr
x = x.*d1_conversion;
for i = 1:size(data,1)
    plot(x,data(i,:),'Color',cmap(i,:));
end

if annotate
    plot(x,data(dataStruct.SarcomereOrientation,:),'Color',[0 0 0],'LineWidth',2);
end

if annotate
    [~,minLocs] = findpeaks(-data(dataStruct.SarcomereOrientation+1,:));
    [~,maxLocs] = findpeaks(data(dataStruct.SarcomereOrientation+1,:));
    
    plot([maxLocs(1).*d1_conversion.*1.15./10 minLocs(1).*d1_conversion./10],data(dataStruct.SarcomereOrientation+1,minLocs(1)).*[1 1],'k:','LineWidth',2);
    plot([maxLocs(1).*d1_conversion.*1.15./10 maxLocs(1).*d1_conversion./10],data(dataStruct.SarcomereOrientation+1,maxLocs(1)).*[1 1],'k:','LineWidth',2);
    
    plot([maxLocs(1).*d1_conversion.*1.15./10 maxLocs(1).*d1_conversion.*1.25./10],data(dataStruct.SarcomereOrientation+1,minLocs(1)).*[1 1],'k','LineWidth',2);
    plot([maxLocs(1).*d1_conversion.*1.15./10 maxLocs(1).*d1_conversion.*1.25./10],data(dataStruct.SarcomereOrientation+1,maxLocs(1)).*[1 1],'k','LineWidth',2);
    plot(maxLocs(1).*d1_conversion./10.*[1 1],[data(dataStruct.SarcomereOrientation+1,maxLocs(1)).*1.05 data(dataStruct.SarcomereOrientation+1,maxLocs(1)).*1.2],'k','LineWidth',2);
end
set(gca,'FontName','Arial');
set(gca,'FontSize',12);
set(gca,'Linewidth',1);
xlabel('Offset distance');
ylabel('Haralick correlation');
set(gca,'TickDir','Out');
set(gca,'XColor','k','YColor','k');
set(gca,'Color','None');
%xlim([0 15]);
ylim([-0.2 1]);
set(gca,'TickLength',[0.0125 0.025]);
if annotate
    cmap(dataStruct.SarcomereOrientation-1,:) = [0 0 0];
    cmap(dataStruct.SarcomereOrientation,:) = [0 0 0];
    cmap(dataStruct.SarcomereOrientation+1,:) = [0 0 0];
end
colormap(cmap);
b = colorbar('Ticks',[0 0.25 0.5 0.75 1],'TickLabels',{'0°','45°','90°','135°','180°'});
b.LineWidth = 1;
b.TickLength = 0.015;
b.Color = 'k';
end