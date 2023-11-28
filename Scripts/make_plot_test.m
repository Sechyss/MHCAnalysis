clear
clc
clf

load('datasets_1.mat')

mal_rec_max=max(mal_rec,[],2);
non_mal_rec_max=max(non_mal_rec,[],2);


ToPlot_total=[];
ToPlot_noMal=[];
ToPlot_withMal=[];

for i=1:1:114

ToPlot_noMal=[ToPlot_noMal;[ones(non_mal_rec_max(i),1).*i,[1:1:non_mal_rec_max(i)]',ones(non_mal_rec_max(i),1).*sim_hum_gen(i)]];

ToPlot_withMal=[ToPlot_withMal;[ones(mal_rec_max(i),1).*i,[1:1:mal_rec_max(i)]',ones(mal_rec_max(i),1).*sim_hum_gen(i)]];

ToPlot_total=[ToPlot_total;[ones(totalvariants(i),1).*i,[1:1:totalvariants(i)]',ones(totalvariants(i),1).*sim_hum_gen(i)]];


end

subplot(2,1,1)

scatter(ToPlot_total(:,1),ToPlot_total(:,2),[],ToPlot_total(:,3));

hold on

scatter(ToPlot_withMal(:,1),ToPlot_withMal(:,2),[],ToPlot_withMal(:,3),'filled');

xticks(1:10:114)
xticklabels(274+5:10:274+114+5)
cb=colorbar();
ylabel(cb,'% Similarity to human peptidome','Rotation',90)

title ('With historical malaria exposure')
ylabel({'Number of variant 11mers';'(filled = recognised'; 'by high frequency HLAs)'},'Rotation',0)

subplot(2,1,2)

scatter(ToPlot_total(:,1),ToPlot_total(:,2),[],ToPlot_total(:,3));

hold on

scatter(ToPlot_noMal(:,1),ToPlot_noMal(:,2),[],ToPlot_noMal(:,3),'filled');


xticks(1:10:114)
xticklabels(274+5:10:274+114+5)

xlabel ('Amino acid position of centre of 11mer (3D7)')
ylabel({'Number of variant 11mers';'(filled = recognised'; 'by high frequency HLAs)'},'Rotation',0)

title ('No historical malaria exposure')
cb=colorbar();
ylabel(cb,'% Similarity to human peptidome','Rotation',90)