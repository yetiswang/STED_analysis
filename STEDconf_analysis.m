%% These codes were developed by Yuyang Wang to analyze the data presented in the ACS paper Multicolor Super-resolution microscopy of protein corona on single nanoparticles.
%  The following codes involve full image analysis of STED and confocal microscopy data collected from Abberior Expertline STED microscope and Imspector software. 
%  A simple outline of the script is as below:
%  ----- Image import with Bioformat Toolbox for Matlab, extraction of useful metadata for automated plotting and calculation
%  ----- Multicolor image plotting 
%  ----- Single particle recognition based on circle detection algorithm and near-neighbour distance filtering. A visualization handle is used for ease of evaluation
%  ----- Cross-talk correction based on experimentally determined parameters and generation of corrected results.
% For any questions, please contact Yuyang at y.wang8@tue.nl. 


%Load image from file and get metadata
data = bfopen();
%dir = '';
%mk dir
num_stacks = size(data,1);


%%
for i = 1 : num_stacks
    metadata{i} = data{i,4};
    
    name{i} = char(metadata{i}.getImageName(i - 1));
    % Print all file names and select images for analysis
    disp(['Frame ', num2str(i), ' is ',name{i}])
    
    stackSizeX(i) = metadata{i}.getPixelsSizeX(i - 1).getValue(); % image width, pixels
    stackSizeY(i) = metadata{i}.getPixelsSizeY(i - 1).getValue(); % image height, pixels
    stackSizeZ(i) = metadata{i}.getPixelsSizeZ(i - 1).getValue(); % number of Z slices
    
    stacks(i) = data{i,1}(1,1);
    
    voxelSizeXdefaultValue(i) = metadata{i}.getPixelsPhysicalSizeX(i - 1).value(ome.units.UNITS.NANOMETER);           % returns value in default unit
    %voxelSizeXdefaultUnit{i} = char(metadata{i}.getPixelsPhysicalSizeX(i - 1).unit().getSymbol()); % returns the default unit type. Note: .st files use mm as unit, instead of anstrom.
    voxelSizeXdouble(i) = double(voxelSizeXdefaultValue(i));
    
    
end 
%% 
% Visualize 3 overlaid images

stks2overlay = [2, 4, 6]; 

% show images separately with scales

figure
subplot(1,3,1)
imshow(uint8((stacks{stks2overlay(1)})))
colorbar
subplot(1,3,2)
imshow(uint8((stacks{stks2overlay(2)})))
colorbar
subplot(1,3,3)
imshow(uint8((stacks{stks2overlay(3)})))
colorbar
saveas(gcf,'stacks.fig')
saveas(gcf,'stacks.png')

% show images overlay
rgbImage = cat(3, uint8(stacks{stks2overlay(1)}),uint8(stacks{stks2overlay(2)}),uint8(stacks{stks2overlay(3)}));
figure
imshow(rgbImage)
axis on
saveas(gcf, 'fused RGB.fig')
saveas(gcf, 'fused RGB.png')

%% Find single particles and aggregates

% Find all particles using imfindcircles
[centers, radii, metric] = imfindcircles(stacks{stks2overlay(2)},[3 8], 'method','PhaseCode', 'Sensitivity', 1, 'EdgeThreshold', 0.05);
%centersStrong5 = centers(1:end,:); 
%radiiStrong5 = radii(1:end);
%metricStrong5 = metric(1:end);
figure
imshow(rgbImage)
axis on
viscircles(centers, radii,'EdgeColor','b');
hold on 
save centers centers
save radii radii
save metric metric

%% Categorize single and aggregates based on near-neighbor distances among % centers. 
figure
scatter(centers(:,1), centers(:,2)) 

for i = 1 : length(centers)
    % label particles
    text(centers(i,1) + 10 ,centers(i,2), num2str(i))
    fullpool = centers(:,:); 
    fullpool (i,:) = [0,0]; 
    querypoints(i,:) = centers(i,:);
    Idx(i) = knnsearch(fullpool, querypoints(i,:));
    hold on 
    scatter(centers(Idx(i),1),centers(Idx(i),2),'r.')
    text(centers(Idx(i),1)+ 20, centers(Idx(i),2) + 10, ['np ', num2str(i)])
    nndist(i) = sqrt((centers(i,1) - centers(Idx(i),1))^2 + (centers(i,2) - centers(Idx(i),2))^2 );
    
end 

saveas(gcf,'nn plots.fig')
saveas(gcf,'nn plots.png')
%% show near neighbor distances
figure
histogram(nndist.*voxelSizeXdouble(stks2overlay(1)), 50 )
xlabel('NN distances (nm)')
ylabel('Occurrence')
saveas(gcf,'nn histogram.fig' )
saveas(gcf,'nn histogram.png')
%% 
% filter out clusters and particles close to edges
ROI_size = 7 ; 
edge = ROI_size/2; 
singles = find(nndist >= 20);
% From singles set filter out particles to close to edge
xb = [ edge stackSizeX(stks2overlay(1)) - edge stackSizeX(stks2overlay(1)) - edge edge ];
yb = [ stackSizeY(stks2overlay(1)) - edge stackSizeY(stks2overlay(1)) - edge edge edge ]; 
not_on_edge_singles = inpolygon( centers(singles,1),centers(singles,2), xb, yb );

singles = singles(not_on_edge_singles); 

figure
imshow(rgbImage)
axis on
hold on
viscircles(centers(singles,:),radii(singles), 'EdgeColor','b')
for i = 1 : length(singles)
    text(centers(singles(i),1)+ 10, centers(singles(i),2), num2str(singles(i)),'Color','w')
end 
hold off
saveas(gcf,'singles_not_on_edge.fig')
saveas(gcf,'singles_not_on_edge.png')
% Plot sizes of single particles
figure
histogram(radii(singles).*2.*voxelSizeXdouble(stks2overlay(1)),10 )
xlabel('Diameter (nm)')
ylabel('Occurrence')
saveas(gcf,'sizes.fig')
saveas(gcf,'sizes.png')

%% Image partitioning based on single particle localizations

for j = stks2overlay
    for i = 1 : length(singles)
        pcle{1,j}(:,:,i) = define_ROI(stacks{1,j},centers(singles(i),:),ROI_size);
    end
end

save pcle_ROIs pcle

%% Plot 3 particles in every channel 
for i = 1 : 3
rand_array(1,i) = randi(length(singles));
end 
pcles2plot = rand_array;
%pcles2plot = [4 8 13];
% Scaled images plot
figure
ax1 = subplot( 3,3,1);
imagesc(pcle{1,stks2overlay(1)}(:,:,pcles2plot(1)))
colorbar
title(['Pcle Num.: ', num2str(singles(pcles2plot(1)))])

ax2 = subplot(3,3,2); 
imagesc(pcle{1,stks2overlay(2)}(:,:,pcles2plot(1)))
colorbar

ax3 = subplot(3,3,3);
imagesc(pcle{1,stks2overlay(3)}(:,:,pcles2plot(1)))
colorbar


ax4 = subplot( 3,3,4);
imagesc(pcle{1,stks2overlay(1)}(:,:,pcles2plot(2)))
colorbar
title(['Pcle Num.: ', num2str(singles(pcles2plot(2)))])

ax5 = subplot(3,3,5);
imagesc(pcle{1,stks2overlay(2)}(:,:,pcles2plot(2)))
colorbar

ax6 = subplot(3,3,6);
imagesc(pcle{1,stks2overlay(3)}(:,:,pcles2plot(2)))
colorbar


ax7 = subplot( 3,3,7); 
imagesc(pcle{1,stks2overlay(1)}(:,:,pcles2plot(3)))
colorbar
title(['Pcle Num.: ', num2str(singles(pcles2plot(3)))])

ax8 = subplot(3,3,8); 
imagesc(pcle{1,stks2overlay(2)}(:,:,pcles2plot(3)))
colorbar

colorbar
ax9 = subplot(3,3,9);
imagesc(pcle{1,stks2overlay(3)}(:,:,pcles2plot(3)))
colorbar

colormap(ax1,'hot')
colormap(ax4,'hot')
colormap(ax7,'hot')

colormap(ax2,'summer')
colormap(ax5,'summer')
colormap(ax8,'summer')

colormap(ax3,'cool')
colormap(ax6,'cool')
colormap(ax9,'cool')

saveas(gcf,'pcles_ROIs.fig')
saveas(gcf,'pcles_ROIs.png')
%% Scaled images plot
figure
ax1 = subplot( 3,3,1);
surf(pcle{1,stks2overlay(1)}(:,:,pcles2plot(1)),'FaceColor','Interp')
colorbar
title(['Pcle Num.: ', num2str(singles(pcles2plot(1)))])

ax2 = subplot(3,3,2); 
surf(pcle{1,stks2overlay(2)}(:,:,pcles2plot(1)),'FaceColor','Interp')
colorbar

ax3 = subplot(3,3,3);
surf(pcle{1,stks2overlay(3)}(:,:,pcles2plot(1)),'FaceColor','Interp')
colorbar


ax4 = subplot( 3,3,4);
surf(pcle{1,stks2overlay(1)}(:,:,pcles2plot(2)),'FaceColor','Interp')
colorbar
title(['Pcle Num.: ', num2str(singles(pcles2plot(2)))])

ax5 = subplot(3,3,5);
surf(pcle{1,stks2overlay(2)}(:,:,pcles2plot(2)),'FaceColor','Interp')
colorbar

ax6 = subplot(3,3,6);
surf(pcle{1,stks2overlay(3)}(:,:,pcles2plot(2)),'FaceColor','Interp')
colorbar


ax7 = subplot( 3,3,7); 
surf(pcle{1,stks2overlay(1)}(:,:,pcles2plot(3)),'FaceColor','Interp')
colorbar
title(['Pcle Num.: ', num2str(singles(pcles2plot(3)))])

ax8 = subplot(3,3,8); 
surf(pcle{1,stks2overlay(2)}(:,:,pcles2plot(3)),'FaceColor','Interp')
colorbar

colorbar
ax9 = subplot(3,3,9);
surf(pcle{1,stks2overlay(3)}(:,:,pcles2plot(3)),'FaceColor','Interp')
colorbar

colormap(ax1,'hot')
colormap(ax4,'hot')
colormap(ax7,'hot')

colormap(ax2,'summer')
colormap(ax5,'summer')
colormap(ax8,'summer')

colormap(ax3,'cool')
colormap(ax6,'cool')
colormap(ax9,'cool')

saveas(gcf,'pcles_ROIs_surf.fig')
saveas(gcf,'pcles_ROIs_surf.png')

%% Convert Cartesian images to polar images
% % cartesian meshing 
% [col, row] = size(stacks{stks2overlay(1)}); 
% x = - col/2 : col/2; 
% y = - row/2 : row/2; 
% [X, Y] = meshgrid(x,y);
% figure
% imagesc(x,y,pcle{1,stks2overlay(1)}(:,:,pcles2plot(1)))
% [theta,rho] = cart2pol(x,y);
% mesh(x,y,pcle{1,stks2overlay(1)}(:,:,pcles2plot(1)))

%% Intensity analysis 

% Red channel mean intensities
reds = mean(pcle{1,stks2overlay(1)}(:,:,:),[1,2]);

oranges = mean(pcle{1,stks2overlay(2)}(:,:,:),[1,2]);

blues = mean(pcle{1,stks2overlay(3)}(:,:,:),[1,2]);

figure
histogram(reds,'FaceColor','red', 'BinWidth', 5)
hold on
histogram(oranges, 'FaceColor','yello','BinWidth', 5)
histogram(blues, 'FaceColor','blue','BinWidth', 5)
legend('BSA(Star Red)','IgG(Star Orange)','Tf(Chromeo 494)')
xlabel('Average Photon counts')
ylabel('Occurrence')

saveas(gcf,'intensity_dist.fig')
saveas(gcf,'intensity_dist.png')
figure
subplot(1,3,1)
scatter(reds,oranges)
xlabel('BSA(star red) intensity')
ylabel('IgG(star orange) intensity')
%xlabel('IgG(star red) intensity')
%ylabel('Tf(star orange) intensity')
subplot(1,3,2)
scatter(reds,blues)
xlabel('IgG(star red) intensity')
ylabel('Tf(Chromeo 494) intensity')
%xlabel('BSA(star red) intensity')
%ylabel('Tf(Chromeo 494) intensity')
subplot(1,3,3)
scatter(blues,oranges)
xlabel('Tf(Chromeo 494) intensity')
ylabel('IgG(star orange) intensity')
%xlabel('BSA(Chromeo 494) intensity')
%ylabel('Tf(star orange) intensity')
saveas(gcf,'correlation.fig')
saveas(gcf,'correlation.png')
% plot normalized intensity plot
reds_norm = reds./max(reds);
oranges_norm = oranges./max(oranges);
blues_norm = blues./max(blues);

figure
histogram(reds_norm,'FaceColor','red', 'BinWidth', 0.01)
hold on
histogram(oranges_norm, 'FaceColor','yello','BinWidth', 0.01)
histogram(blues_norm, 'FaceColor','blue','BinWidth', 0.01)
legend('BSA(Star Red)','IgG(Star Orange)','Tf(Chromeo 494)')
%Tlegend('IgG(Star Red)','Tf(Star Orange)','BSA(Chromeo 494)')

xlabel('Normalized Photon counts')
ylabel('Occurrence')
saveas(gcf,'norm_distribution.fig')
saveas(gcf,'norm_distribution.png')
% plot normalized correlation
figure
subplot(1,3,1)
scatter(reds_norm,oranges_norm)
xlabel('BSA(star red) intensity')
ylabel('IgG(star orange) intensity')
%xlabel('IgG(star red) intensity')
%ylabel('Tf(star orange) intensity')
subplot(1,3,2)
scatter(reds_norm,blues_norm)
%xlabel('IgG(star red) intensity')
%ylabel('BSA(Chromeo 494) intensity')
xlabel('BSA(star red) intensity')
ylabel('Tf(Chromeo 494) intensity')
subplot(1,3,3)
scatter(blues_norm,oranges_norm)
%xlabel('BSA(Chromeo 494) intensity')
%ylabel('IgG(star red) intensity')
xlabel('BSA(star red) intensity')
ylabel('Tf(Chromeo 494) intensity')

saveas(gcf,'norm_correlation.fig')
saveas(gcf,'norm_correlation.png')

save reds reds
save oranges oranges
save blues blues
save stacks stacks

%% Calculate cross talk corrected intensities
load('C:\Users\Yuyang\Dropbox\Postdoc\Publications\STED nanoparticles\data\20210302 1C crosstalks\CFmatrix.mat');
intens_matrix = [squeeze(reds) squeeze(oranges) squeeze(blues)]';
intenscor_matrix = inv(CF)*intens_matrix; 

figure
sgtitle('Corrected intensities')
subplot(1,3,1)
scatter(intenscor_matrix(1,:),intenscor_matrix(2,:))
ylim([0 max(intenscor_matrix(1,:))])
xlabel('Red')
ylabel('Orange')
subplot(1,3,2)
scatter(intenscor_matrix(1,:),intenscor_matrix(3,:))
ylim([-5 max(intenscor_matrix(1,:))])
xlabel('Red')
ylabel('Blue')
subplot(1,3,3)
scatter(intenscor_matrix(2,:),intenscor_matrix(3,:))
ylim([-5 max(intenscor_matrix(1,:))])
ylim([-5 max(intenscor_matrix(1,:))])
xlabel('Orange')
ylabel('Blue')
set(gcf,'position',[100, 100, 2000, 500])
saveas(gcf,'corrected intensity correlations.fig')
saveas(gcf,'corrected intensity correlations.png')

figure
histogram(intenscor_matrix(1,:),'FaceColor','red','BinWidth', 2)
hold on
histogram(intenscor_matrix(2,:), 'FaceColor','yello','BinWidth', 2)
histogram(intenscor_matrix(3,:), 'FaceColor','blue','BinWidth', 2)
legend('BSA(Star Red)','IgG(Star Orange)','Tf(Chromeo 494)')
%legend('IgG(Star Red)','Tf(Star Orange)','BSA(Chromeo 494)')
xlabel('Corrected photon counts')
ylabel('Occurrence')
saveas(gcf,'corintens_distribution.fig')
saveas(gcf,'corintens_distribution.png')

save cor_intensity_matrix intenscor_matrix

%% Correct further with extinction coefficients of dyes, DOL and quantum yield

% !!! FOR labeling strategy I !!!
ext_red = 120000; % M-1cm-1
ext_orange = 95000; 
ext_chromeo = 55000; 
phi_red = 0.55;
phi_orange = 0.55;
phi_chromeo = 0.15;
DOL_red = 1.1;
DOL_orange = 2.3;
DOL_chromeo = 1;
%DOL_red = 0.67;
%DOL_orange = 2.5;
%DOL_chromeo = 0.99;
intenscor_matrix(3,intenscor_matrix(3,:)<=0) = NaN;
cormat = [ intenscor_matrix(1,:)./(ext_red*phi_red*DOL_red); ...
    intenscor_matrix(2,:)./(ext_orange*phi_orange*DOL_orange); ...
    intenscor_matrix(3,:)./(ext_chromeo*phi_chromeo*DOL_chromeo)];

% figure 
% subplot(1,3,1)
% scatter(intenscor_matrix(1,:),intenscor_matrix(2,:))
% subplot(1,3,2)
% scatter(intenscor_matrix(2,:),intenscor_matrix(3,:))
% subplot(1,3,3)
% scatter(intenscor_matrix(1,:),intenscor_matrix(3,:))
%including correction for laser power
%cormat = [ intenscor_matrix(1,:)./(ext_red*phi_red*DOL_red); ...
%    intenscor_matrix(2,:)./(ext_orange*phi_orange*DOL_orange).*0.93/0.41; ...
%    intenscor_matrix(3,:)./(ext_chromeo*phi_chromeo*DOL_chromeo).*0.93/0.61];

figure
histogram(cormat(1,:))
hold on 
histogram(cormat(2,:))
histogram(cormat(3,:))

figure
subplot(1,3,1)
scatter(cormat(1,:),cormat(2,:))
subplot(1,3,2)
scatter(cormat(2,:),cormat(3,:))
subplot(1,3,3)
scatter(cormat(1,:),cormat(3,:))

Sum = cormat(1,:) + cormat(2,:) + cormat(3,:);
Red_perc = cormat(1,:)./Sum;
Orange_perc = cormat(2,:)./Sum; 
Blue_perc = cormat(3,:)./Sum ; 

figure
histogram(Red_perc,'FaceColor','red', 'Binwidth', 0.02)
hold on 
histogram(Orange_perc,'FaceColor','yellow', 'Binwidth', 0.02)
histogram(Blue_perc,'FaceColor','blue', 'Binwidth', 0.02)
legend('BSA(Star Red)','IgG(Star Orange)','Tf(Chromeo 494)')
xlabel('Protein content')
ylabel('Occurrence')
pbaspect([1.2 1 1])
saveas(gcf,'finalcornorm_distribution.fig')
saveas(gcf,'finalcornorm_distribution.png')


figure
subplot(1,3,1)
scatter(Red_perc,Orange_perc)
subplot(1,3,2)
scatter(Orange_perc,Blue_perc)
subplot(1,3,3)
scatter(Red_perc,Blue_perc)

save cormat cormat

%save cormat_norm cormat_norm

%% !!! FOR labeling strategy II !!!
ext_red = 120000; % M-1cm-1
ext_orange = 95000; 
ext_chromeo = 55000; 
phi_red = 0.55;
phi_orange = 0.55;
phi_chromeo = 0.15;
% DOL_red = 2.86;
% DOL_orange = 1.88;
% DOL_chromeo = 2.98;
DOL_red = 2.2;
DOL_orange = 1.6;
DOL_chromeo = 3;
%intenscor_matrix(3,intenscor_matrix(3,:)<=0) = NaN;
 cormat = [ intenscor_matrix(1,:)./(ext_red*phi_red*DOL_red); ...
     intenscor_matrix(2,:)./(ext_orange*phi_orange*DOL_orange); ...
     intenscor_matrix(3,:)./(ext_chromeo*phi_chromeo*DOL_chromeo)];


figure
histogram(cormat(1,:), 'Binwidth', 0.0001)
hold on 
histogram(cormat(2,:), 'Binwidth', 0.0001)
histogram(cormat(3,:), 'Binwidth', 0.0001)

figure
histogram(cormat(1,:)./max(max(cormat)),'FaceColor','red', 'Binwidth', 0.02)
hold on 
histogram(cormat(2,:)./max(max(cormat)),'FaceColor','yellow', 'Binwidth', 0.02)
histogram(cormat(3,:)./max(max(cormat)),'FaceColor','blue', 'Binwidth', 0.02)
legend('IgG(Star Red)','Tf(Star Orange)','BSA(Chromeo 494)')
xlabel('Protein content')
ylabel('Occurrence')
pbaspect([1.2 1 1])
saveas(gcf,'finalcornorm_distribution.fig')
saveas(gcf,'finalcornorm_distribution.png')


figure
subplot(1,3,1)
scatter(cormat(1,:),cormat(2,:))
subplot(1,3,2)
scatter(cormat(2,:),cormat(3,:))
subplot(1,3,3)
scatter(cormat(1,:),cormat(3,:))

save cormat cormat


%% Plot cytofluorograms for all three channels

h = figure; 
subplot(1,3,1)
scatter(stacks{1,4}(:),stacks{1,5}(:))
ylim([0,max(stacks{1,4}(:))])
xlabel('Blue')
ylabel('Red')
title(['Ratio of average: ', num2str(mean(mean(stacks{1,5}))/mean(mean(stacks{1,4})))])
subplot(1,3,2)
scatter(stacks{1,4}(:),stacks{1,6}(:))
ylim([0,max(stacks{1,4}(:))])
xlabel('Blue')
ylabel('Orange')
title(['Ratio of average: ', num2str(mean(mean(stacks{1,6}))/mean(mean(stacks{1,4})))])
subplot(1,3,3)
scatter(stacks{1,5}(:),stacks{1,6}(:))
xlim([0,max(stacks{1,4}(:))])
ylim([0,max(stacks{1,4}(:))])
xlabel('Red')
ylabel('Orange')
sgtitle('Cytofluorograms Chromeo 494 only')
set(gcf,'position',[100, 100, 2000, 500])
saveas(gcf,'cytofluorograms.fig')
saveas(gcf,'cytofluorograms.png')

%% Plot cytofluorograms for single particles


h = figure; 
subplot(1,3,1)
reds1 = squeeze(reds1);
reds2 = squeeze(reds2);
reds3 = squeeze(reds3);
reds4 = squeeze(reds4);
oranges1 = squeeze(oranges1);
oranges2 = squeeze(oranges2);
oranges3 = squeeze(oranges3);
oranges4 = squeeze(oranges4);
blues1 = squeeze(blues1);
blues2 = squeeze(blues2);
blues3 = squeeze(blues3);
blues4 = squeeze(blues4);
reds = [reds1; reds2 ; reds3; reds4 ];
oranges = [oranges1 ; oranges2 ; oranges3; oranges4];
blues = [blues1 ; blues2; blues3; blues4]; 
scatter(reds,oranges)
ylim([0,max(reds)])
xlabel('Red')
ylabel('Orange')
title(['Ratio of average: ', num2str(mean(oranges)/mean(reds))])
subplot(1,3,2)
scatter(reds,blues)
ylim([0,max(reds)])
xlabel('Red')
ylabel('Blue')
title(['Ratio of average: ', num2str(mean(blues)/mean(reds))])
subplot(1,3,3)
scatter(oranges,blues)
xlim([0,max(reds)])
ylim([0,max(reds)])
xlabel('Orange')
ylabel('Blue')
sgtitle('Cytofluorograms Star Red only')
set(gcf,'position',[100, 100, 2000, 500])
saveas(gcf,'cytofluorograms.fig')
saveas(gcf,'cytofluorograms.png')
