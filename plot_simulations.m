%% Clean workspace, and open files


clear all
close all
clc


% Select simulation to load
k = 2;

disp([num2str(k)])

% Check how many files are inside simulation folder
a = dir(['Sim' num2str(k) '/Gen_space*.txt']);
Size = numel(a);

% Set number of voxels per dimension
N = 80;
% Set carrying capacity
K = 2e5;
% And threshold number of cells to consider a voxel occupied
threshold = 0.2*K;

% Create cell to store data at all time steps
A = cell(1,Size);

% Change to 1 if you want to save videos
videos = 0;


for i = 1:Size
    
    disp(num2str(i))
    f = ['Sim' num2str(k) '/Gen_space_' num2str(i*20) '.txt'];
    delimiterIn = ' ';
    headerlinesIn = 0;
    A{i} = importdata(f,delimiterIn,headerlinesIn);
    
end


% Decode data into 3d structures
poptot = cell(1,Size);
pops = cell(1,Size);
nec = cell(1,Size);
act = cell(1,Size);

for i = 1:Size
    
    disp(num2str(i))
    poptot{i} = zeros(N,N,N);
    pops{i} = zeros(N,N,N,8);
    nec{i} = zeros(N,N,N);
    act{i} = zeros(N,N,N);

   
    
    for e = 1:8
        % Coordinates of occupied voxels
        if A{i}(find(A{i}(:,e+3)),e+3)
            occ = A{i}(find(A{i}(:,e+3)),[1:3]); 
            index = length(occ(:,1));
            
            x = occ(:,1);
            y = occ(:,2);
            z = occ(:,3); 
            
            popgen = A{i}(find(A{i}(:,e+3)),e+3);         
            for j = 1:index
                
                % Save info about cell number for each cell population in each voxel
                pops{i}(x(j),y(j),z(j),e) = popgen(j);
                % And for total cell number at each voxel
                poptot{i}(x(j),y(j),z(j)) = poptot{i}(x(j),y(j),z(j)) + pops{i}(x(j),y(j),z(j),e);
                
            end
            
        end

    end

    % Save info about cell activity
    occ = A{i}(find(A{i}(:,12)),[1:3]);  
    x = occ(:,1);
    y = occ(:,2);
    z = occ(:,3);
    actvox = A{i}(find(A{i}(:,12)),12);
    for j = 1:length(occ(:,1))
        act{i}(x(j),y(j),z(j)) = actvox(j);
    end
    
    % And also about necrotic cells
    occ = A{i}(find(A{i}(:,12)),[1:3]);  
    x = occ(:,1);
    y = occ(:,2);
    z = occ(:,3);
    necvox = A{i}(find(A{i}(:,12)),12);
    for j = 1:length(occ(:,1))
        nec{i}(x(j),y(j),z(j)) = necvox(j);
    end
    

    
end

% Get volume (number of occupied voxels)
VOL2 = zeros(1,Size);
% And total activity at each time step
totact = zeros(1,Size);

for i = 1:Size
    mascP{i} = poptot{i}+nec{i} > threshold;
    VOL2(i) = sum(sum(sum(mascP{i})));
    totact(i) = sum(sum(sum(act{i})));
end

% Uncomment code below if you want to save these variables into files
% fileID = fopen(['Sim' num2str(k) '/Vol_matlab.txt'],'w');
% fprintf(fileID,'%8.2f\n',VOL2);
% fclose(fileID);
% 
% fileID = fopen(['Sim' num2str(k) '/Act_matlab.txt'],'w');
% fprintf(fileID,'%8.2f\n',totact);
% fclose(fileID);


%%

% HEATMAP - TABLE

popg = zeros(8,Size);

for i = 1:Size
    popg(:,i) = sum(sum(sum(pops{i})));
end




%% ACTIVITY SLICE

figure()
subplot(221)
imagesc(act{Size-20}(:,:,40))
colormap gray
axis equal
axis([0 80 0 80])
caxis([0 5e3])
% title('Activity','Fontsize',20)
colorbar

subplot(222)
imagesc(act{Size-15}(:,:,40))
colormap gray
axis equal
axis([0 80 0 80])
caxis([0 5e3])
colorbar


subplot(223)
imagesc(act{Size-10}(:,:,40))
colormap gray
axis equal
axis([0 80 0 80])
caxis([0 5e3])
colorbar


subplot(224)
imagesc(act{Size}(:,:,40))
colormap gray
axis equal
axis([0 80 0 80])
caxis([0 5e3])
colorbar

% print(gcf,['Sim' num2str(k) '/Activity_lobules.png'],'-dpng','-r600');



%% TIME OF APPEARANCE OF EACH MUTATION

popT = zeros(Size,8);

for i = 1:Size
    popG = squeeze(sum(pops{i},1));
    popG = squeeze(sum(popG,1));
    popG = squeeze(sum(popG,1));
    popT(i,:) = popG;
end

Tp53 = find(popT(:,2),1);
Trb = find(popT(:,3),1);
Tp53rb = find(popT(:,4),1);
Trtk = find(popT(:,5),1);
Trtkp53 = find(popT(:,6),1);
Trtkrb = find(popT(:,7),1);
Tall = find(popT(:,8),1);

if isempty(Tall)
    Tall = 200;
end

if isempty(Tp53rb)
    Tp53rb = 200;
end


%% HETEROGENEITY

shannon = zeros(1,Size);
simpson = zeros(1,Size);

% popG = squeeze(sum(pops{Size},1));
% popG = squeeze(sum(popG,1));
% popG = squeeze(sum(popG,1));

for i = 1:Size
    popG = squeeze(sum(pops{i},1));
    popG = squeeze(sum(popG,1));
    popG = squeeze(sum(popG,1));
    
    for e = 1:8
        if popG(e) > 0
            shannon(i) = shannon(i) - popG(e)/sum(popG)*log(popG(e)/sum(popG));
            simpson(i) = simpson(i) + (popG(e)/sum(popG))^2; 
        end
    end
end

time = linspace(0,100,length(shannon));

figure()
hold on
plot(time,shannon,'Linewidth',2)
plot(time,simpson,'Linewidth',2)
plot([Tp53 Tp53], [0 1.2],'k--');
text(Tp53-2,0.9,'P53')
plot([Trb Trb], [0 1.2],'k--');
text(Trb-2,0.85,'Rb')
plot([Tp53rb Tp53rb], [0 1.2],'k--');
text(Tp53rb-4,0.9,'P53+Rb')
plot([Trtk Trtk], [0 1.2],'k--');
text(Trtk-2,0.8,'RTK')
plot([Trtkp53 Trtkp53], [0 1.2],'k--');
text(Trtkp53-4,0.8,'P53+RTK')
plot([Trtkrb Trtkrb], [0 1.2],'k--');
text(Trtkrb-4,0.85,'Rb+RTK')
plot([Tall Tall], [0 1.2],'k--');
text(Tall-2,0.8,'ALL')
xlabel('Time [%]','FontSize',20)
hold off
set(gcf,'color','w');
set(gca,'Fontsize',18);
% print(gcf,['Sim' num2str(k) '/Heterogeneity.png'],'-dpng','-r600');


%% MORPHOLOGY MEASUREMENTS

% K = 1e3/0.2;
K = 2e5;
threshold = 0.2*K;
mascP = cell(1,Size);
stats = cell(1,Size);
SPHER = zeros(1,Size);
SOLID = zeros(1,Size);
VOL = zeros(1,Size);
CONVOL = zeros(1,Size);
VOL_SPHERE = zeros(1,Size);
Centroid = zeros(3,Size);
Ncontr = zeros(1,Size);
SUVmax = zeros(1,Size);
iSUVmax = zeros(1,Size);
i3dSUVmax = zeros(3,Size);
distCentSUV = zeros(1,Size);
distCentSUVref = zeros(1,Size);
RME = zeros(1,Size);
SURFACE = zeros(1,Size);
SURFACE_REG = zeros(1,Size);




PX = 1; % mm
PY = 1; % mm
PZ = 1; % mm

m = 80;
n = 80;
p = 80;

[mallaX,mallaY,mallaZ] = meshgrid(PX/2 : PX : (PX*m-PX/2), ...
                                  PY/2 : PY : (PY*n-PY/2), ...
                                  PZ/2 : PZ : (PZ*p-PZ/2));



% figure()
% plot(VOL2)

minVol = 1e2;
ind = find(VOL2 > minVol);
limit = ind(1);


for i = limit:Size
    mascP{i} = poptot{i}+nec{i} > threshold;
    
    Ncontr(i) = sum(sum(sum(mascP{i})));
    
    Centroid(:,i) = [sum(sum(sum(mallaX .* mascP{i}))), ...
             sum(sum(sum(mallaY .* mascP{i}))), ...
             sum(sum(sum(mallaZ .* mascP{i})))] / Ncontr(i);
         
    
    [SUVmax(i) iSUVmax(i)] = max(act{i}(:));
    i3dSUVmax(:,i) = [mallaX(iSUVmax(i)),mallaY(iSUVmax(i)),mallaZ(iSUVmax(i))];
    
    distCentSUV(i) = sqrt((Centroid(1,i)-i3dSUVmax(1,i))^2 + (Centroid(2,i)-i3dSUVmax(2,i))^2 + (Centroid(3,i)-i3dSUVmax(3,i))^2 );
    
    RME(i) = (3/(4*pi)*Ncontr(i))^(1/3);
    distCentSUVref(i) = distCentSUV(i)/(RME(i));         
       

    stats{i} = regionprops3(mascP{i},'all');
    
    SOLID(i) = stats{i}.Solidity;
    VOL(i) = stats{i}.Volume;
    CONVOL(i) = stats{i}.ConvexVolume;
    VOL_SPHER(i) = 4/3*pi*(stats{i}.SurfaceArea/(4*pi))^(3/2);
    SPHER(i) = VOL(i)/VOL_SPHER(i);
    SURFACE(i) = stats{i}.SurfaceArea;
    SURFACE_REG(i) = 6*sqrt(pi)*(VOL(i)/(stats{i}.SurfaceArea^(3/2)));
end

t = limit:Size;
t2 = 100*t/Size;

figure()
plot(t2,SPHER(t),'ro-','Linewidth',1)
xlabel('Time [%]','Fontsize',18)
ylabel('Sphericity','Fontsize',18)
set(gcf,'color','w');
set(gca,'FontSize',16);
% print(gcf,['Sim' num2str(k) '/Sphericity.png'],'-dpng','-r600');

figure()
plot(t2,VOL2(t),'Linewidth',2)
xlabel('Time [%]','Fontsize',18)
ylabel('Volume [voxels]','Fontsize',18)
set(gcf,'color','w');
set(gca,'FontSize',16);
% print(gcf,['Sim' num2str(k) '/Volumes.png'],'-dpng','-r600');


% distCentSUV = distCentSUV(~isnan(distCentSUV));
figure()
plot(t2,distCentSUV(t),'rx-')
xlabel('Time [%]','FontSize',18)
ylabel('Distance [mm]','FontSize',18)
set(gcf,'color','w');
set(gca,'FontSize',16);
% print(gcf,['Sim' num2str(k) '/distSUVmax_Centroid.png'],'-dpng','-r600');

% distCentSUVref = distCentSUVref(~isnan(distCentSUVref));                  
figure()
plot(t2,distCentSUVref(t),'bx-')
xlabel('Time [%]','FontSize',18)
ylabel('Distance/Radius','FontSize',18)
set(gcf,'color','w');
set(gca,'FontSize',16);
% print(gcf,['Sim' num2str(k) '/distSUVmax_Centroid_ref.png'],'-dpng','-r600');

% RME = RME(find(RME));
figure()
plot(t2,RME(t))
xlabel('Time [%]','FontSize',18)
ylabel('MSR [mm]','FontSize',18)
set(gcf,'color','w');
set(gca,'FontSize',16);
% print(gcf,['Sim' num2str(k) '/MSRadius.png'],'-dpng','-r600');


figure()
plot(t2,SURFACE_REG(t))
xlabel('Time [%]','FontSize',18)
ylabel('Surface Regularity','FontSize',18)
set(gcf,'color','w');
set(gca,'FontSize',16);
% print(gcf,['Sim' num2str(k) '/SurfReg.png'],'-dpng','-r600');


figure()
plot(t2,SOLID(t))
title('SOLIDITY','FontSize',16)

% save(['Sim' num2str(k) '/distSUVmax_Centroid.mat'],'distCentSUV');
% save(['Sim' num2str(k) '/distSUVmax_Centroid_ref.mat'],'distCentSUVref');
% save(['Sim' num2str(k) '/Sphericity.mat'],'SPHER');
% save(['Sim' num2str(k) '/SurfReg.mat'],'SURFACE_REG');


   
%% CELL NUMBER SLICE

if videos == 1

    figure()
    colormap gray
    for i = 1:Size

        h1 = slice(poptot{i}(:,:,:,1),[],40,40);
        set(h1,'FaceColor','interp');
        set(h1,'EdgeColor','none');
        axis equal
        axis tight
        axis([0 80 0 80])
        axis off
        caxis([0 2e5]);
        %text(5,75,40,[num2str(i*10),' days'],'Color','white','FontSize',14)
        title('Cell number','FontSize',20)
        xlabel('x [mm]','FontSize',16)
        ylabel('y [mm]','FontSize',16)
        colorbar
        set(gcf,'color','w');

        drawnow
        G(i) = getframe(gcf);

        %saveas(gcf,'PopTot_slice.png')
    end
    % create the video writer with 1 fps
    %writerObj = VideoWriter('Rho_ActAlpha_Evolution_BoneColor.avi');
    writerObj = VideoWriter(['Sim' num2str(k) '/Tumor_slice.avi']);
    writerObj.FrameRate = 10;
    % set the seconds per image
    % open the video writer
    open(writerObj);
    % write the frames to the video
    for i=1:length(G)
    % convert the image to a frame
    frame = G(i) ;    
    writeVideo(writerObj, frame);
    end
    % close the writer object
    close(writerObj);
    
end


%% CELL NUMBER SLICE (LOG)

if videos == 1

    figure()
    colormap gray
    for i = 1:Size
    %     imagesc(poptot{i}(:,:,40))
    %     axis equal
    %     axis([0 80 0 80])
    %     caxis([0 2e5])
    %     title('Cell number','Fontsize',18)
    %     colorbar

        h1 = slice(log(poptot{i}(:,:,:)),[],40,40);
%         set(h1,'FaceColor','interp');
        set(h1,'EdgeColor','none');
        axis equal
        axis tight
        axis([0 80 0 80])
        axis off
        caxis([0 log(2e5)]);
        %text(5,75,40,[num2str(i*10),' days'],'Color','white','FontSize',14)
        title('Cell number','FontSize',20)
        xlabel('x [mm]','FontSize',16)
        ylabel('y [mm]','FontSize',16)
        colorbar
        set(gcf,'color','w');

        drawnow
        G(i) = getframe(gcf);

        %saveas(gcf,'PopTot_slice.png')
    end
    % create the video writer with 1 fps
    %writerObj = VideoWriter('Rho_ActAlpha_Evolution_BoneColor.avi');
    writerObj = VideoWriter(['Sim' num2str(k) '/Tumor_slice_log.avi']);
    writerObj.FrameRate = 10;
    % set the seconds per image
    % open the video writer
    open(writerObj);
    % write the frames to the video
    for i=1:length(G)
    % convert the image to a frame
    frame = G(i) ;    
    writeVideo(writerObj, frame);
    end
    % close the writer object
    close(writerObj);
    
end


%% CELL NUMBER SLICE ALTERATIONS

if videos == 1

    figpos = get(groot, 'ScreenSize');
    ventana = figure('Name','Evolution of mutations','Units','pixels','OuterPosition',[0 0 figpos(4)/sqrt(2) figpos(4)]);



    colormap parula
    for i = 1:Size

        subplot(3,3,1)
        h1 = slice(pops{i}(:,:,:,1),[],40,40);
        set(h1,'FaceColor','interp');
        set(h1,'EdgeColor','none');
        axis equal
        axis([0 80 0 80])
        axis off
        caxis([0 max(max(pops{i}(:,:,40,1)))]);
        title('Wild-type','FontSize',18)
        xlabel('x [mm]','FontSize',16)
        ylabel('y [mm]','FontSize',16)

        subplot(3,3,2)
        h2 = slice(pops{i}(:,:,:,2),[],40,40);
        set(h2,'FaceColor','interp');
        set(h2,'EdgeColor','none');
        axis equal
        axis([0 80 0 80])
        axis off
        caxis([0 max(max(pops{i}(:,:,40,2)))]);
        title('P53','FontSize',18)

        subplot(3,3,3)
        h3 = slice(pops{i}(:,:,:,3),[],40,40);
        set(h3,'FaceColor','interp');
        set(h3,'EdgeColor','none');
        axis equal
        axis([0 80 0 80])
        axis off
        caxis([0 max(max(pops{i}(:,:,40,3)))]);
        title('Rb','FontSize',18)

        subplot(3,3,4)
        h4 = slice(pops{i}(:,:,:,4),[],40,40);
        set(h4,'FaceColor','interp');
        set(h4,'EdgeColor','none');
        axis equal
        axis([0 80 0 80])
        axis off
        caxis([0 max(max(pops{i}(:,:,40,4)))]);
        title('P53+Rb','FontSize',18)

        subplot(3,3,5)
        h5 = slice(pops{i}(:,:,:,5),[],40,40);
        set(h5,'FaceColor','interp');
        set(h5,'EdgeColor','none');
        axis equal
        axis([0 80 0 80])
        axis off
        caxis([0 max(max(pops{i}(:,:,40,5)))]);
        title('RTK','FontSize',18)

        subplot(3,3,6)
        h6 = slice(pops{i}(:,:,:,6),[],40,40);
        set(h6,'FaceColor','interp');
        set(h6,'EdgeColor','none');
        axis equal
        axis([0 80 0 80])
        axis off
        caxis([0 max(max(pops{i}(:,:,40,6)))]);
        title('P53+RTK','FontSize',18)

        subplot(3,3,7)
        h7 = slice(pops{i}(:,:,:,7),[],40,40);
        set(h7,'FaceColor','interp');
        set(h7,'EdgeColor','none');
        axis equal
        axis([0 80 0 80])
        axis off
        caxis([0 max(max(pops{i}(:,:,40,7)))]);
        title('Rb+RTK','FontSize',18)

        subplot(3,3,8)
        h8 = slice(pops{i}(:,:,:,8),[],40,40);
        set(h8,'FaceColor','interp');
        set(h8,'EdgeColor','none');
        axis equal
        axis([0 80 0 80])
        axis off
        caxis([0 max(max(pops{i}(:,:,40,8)))]);
        title('P53+Rb+RTK','FontSize',18)

        set(gcf,'color','w');

        drawnow
        G(i) = getframe(gcf);

        %saveas(gcf,'PopTot_slice.png')
    end
    % create the video writer with 1 fps
    %writerObj = VideoWriter('Rho_ActAlpha_Evolution_BoneColor.avi');
    writerObj = VideoWriter(['Sim' num2str(k) '/TumorAlterations_slice.avi']);
    writerObj.FrameRate = 10;
    % set the seconds per image
    % open the video writer
    open(writerObj);
    % write the frames to the video
    for i=1:length(G)
    % convert the image to a frame
    frame = G(i) ;    
    writeVideo(writerObj, frame);
    end
    % close the writer object
    close(writerObj);
    
end

%% CELL NUMBER 3D

if videos == 1

    % Generate color palette
    colP53 = [1,0,0];
    colRb = [0,0,1];
    colRTK = [0,1,0];

    figure()
    colormap parula
    for i = 1:Size
        hold on

        pt = patch(isosurface(smooth3(poptot{i}),1e4));
        set(pt,'FaceAlpha',0.8,'EdgeAlpha',0.8,'FaceColor','cyan','EdgeColor','none');

        %vol3d('CData',poptot{i},'Texture','3d')
        view(3)
        axis equal
        axis([0 80 0 80 0 80])
        axis off
        camlight
        lighting phong

        hold off
        drawnow
        G(i) = getframe(gcf);

        %saveas(gcf,'PopTot_slice.png')
    end
    % create the video writer with 1 fps
    % writerObj = VideoWriter('Rho_ActAlpha_Evolution_BoneColor.avi');
    writerObj = VideoWriter(['Sim' num2str(k) '/Tumor_3d.avi'],'Uncompressed AVI');
    % writerObj.Quality = 100;
    writerObj.FrameRate = 10;
    % set the seconds per image
    % open the video writer
    open(writerObj);
    % write the frames to the video
    for i=1:length(G)
    % convert the image to a frame
    frame = G(i) ;    
    writeVideo(writerObj, frame);
    end
    % close the writer object
    close(writerObj);

    % %% CELL NUMBER 3D ALL GENS
    % 
    % % Generate color palette
    % colP53 = [0.9,0.1,0.1];
    % colRb = [0.1,0.1,0.9];
    % colRTK = [0.1,0.9,0.1];
    % 
    % figure()
    % colormap parula
    % for i = 1:Size
    %     hold on
    % %     pt = patch(isosurface(smooth3(poptot{i}),1e4));
    % %     set(pt,'FaceAlpha',0.15,'EdgeAlpha',0.15,'FaceColor','cyan','EdgeColor','none');
    % %     
    %     p1 = patch(isosurface(smooth3(pops{i}(:,:,:,1)),1.1e4));
    %     set(p1,'FaceAlpha',0.9,'EdgeAlpha',0.9,'FaceColor','cyan','EdgeColor','none');
    %     
    %     p2 = patch(isosurface(smooth3(pops{i}(:,:,:,2)),1.1e4));
    %     set(p2,'FaceAlpha',0.9,'EdgeAlpha',0.9,'FaceColor',colP53,'EdgeColor','none');
    %     
    %     p3 = patch(isosurface(smooth3(pops{i}(:,:,:,3)),1.1e4));
    %     set(p3,'FaceAlpha',0.9,'EdgeAlpha',0.9,'FaceColor',colRb,'EdgeColor','none');
    %     
    %     p4 = patch(isosurface(smooth3(pops{i}(:,:,:,4)),1.1e4));
    %     set(p4,'FaceAlpha',0.9,'EdgeAlpha',0.9,'FaceColor',(colP53+colRb)./2,'EdgeColor','none');
    %     
    %     p5 = patch(isosurface(smooth3(pops{i}(:,:,:,5)),1.1e4));
    %     set(p5,'FaceAlpha',0.9,'EdgeAlpha',0.9,'FaceColor',colRTK,'EdgeColor','none');
    %     
    %     p6 = patch(isosurface(smooth3(pops{i}(:,:,:,6)),1.1e4));
    %     set(p6,'FaceAlpha',0.9,'EdgeAlpha',0.9,'FaceColor',(colP53+colRTK)./2,'EdgeColor','none');
    %     
    %     p7 = patch(isosurface(smooth3(pops{i}(:,:,:,7)),1.1e4));
    %     set(p7,'FaceAlpha',0.9,'EdgeAlpha',0.9,'FaceColor',(colRTK+colRb)./2,'EdgeColor','none');
    %     
    %     p8 = patch(isosurface(smooth3(pops{i}(:,:,:,8)),1.1e4));
    %     set(p8,'FaceAlpha',0.9,'EdgeAlpha',0.9,'FaceColor',(colP53+colRTK+colRb)./3,'EdgeColor','none');
    % %     h = slice(poptot{i},N/2,N/2,N/2);
    % %     set(h,'EdgeColor','none','FaceAlpha',0.2);
    % 
    %     %vol3d('CData',poptot{i},'Texture','3d')
    %     view(3)
    %     axis equal
    %     axis([0 80 0 80 0 80])
    %     axis off
    %     camlight
    %     lighting phong
    %     
    %     hold off
    %     drawnow
    %     G(i) = getframe(gcf);
    % 
    %     %saveas(gcf,'PopTot_slice.png')
    % end
    % % create the video writer with 1 fps
    % % writerObj = VideoWriter('Rho_ActAlpha_Evolution_BoneColor.avi');
    % writerObj = VideoWriter('Tumor_3d_AllGens.avi','Archival');
    % % writerObj.Quality = 100;
    % writerObj.FrameRate = 10;
    % % set the seconds per image
    % % open the video writer
    % open(writerObj);
    % % write the frames to the video
    % for i=1:length(G)
    % % convert the image to a frame
    % frame = G(i) ;    
    % writeVideo(writerObj, frame);
    % end
    % % close the writer object
    % close(writerObj);
    
end


%% ACTIVITY SLICE

if videos == 1

    figure()
    colormap gray
    for i = 1:Size
        imagesc(act{i}(:,:,40))
        axis equal
        axis([0 80 0 80])
        caxis([0 1e4])
        title('Activity','Fontsize',18)
        colorbar
        drawnow
        G(i) = getframe(gcf);

        %saveas(gcf,'PopTot_slice.png')
    end
    % create the video writer with 1 fps
    %writerObj = VideoWriter('Rho_ActAlpha_Evolution_BoneColor.avi');
    writerObj = VideoWriter(['Sim' num2str(k) '/Activity_slice.avi']);
    writerObj.FrameRate = 10;
    % set the seconds per image
    % open the video writer
    open(writerObj);
    % write the frames to the video
    for i=1:length(G)
    % convert the image to a frame
    frame = G(i) ;    
    writeVideo(writerObj, frame);
    end
    % close the writer object
    close(writerObj);

end





%% SINGLE SLICE 

% Custom maps
bluemap = [0 0 0 % For P53
    0 0 .1
    0 0 .2
    0 0 .3
    0 0 .4
    0 0 .5
    0 0 .6
    0 0 .7
    0 0 .8
    0 0 .9
    0 0 1];
% bluemap(:,1:2) = bluemap(:,1:2)+1;

redmap = [0 0 0 % For Rb
   .1 0 0 
   .2 0 0 
   .3 0 0 
   .4 0 0 
   .5 0 0 
   .6 0 0 
   .7 0 0 
   .8 0 0 
   .9 0 0 
   1 0 0 ];
% redmap(:,2:3) = redmap(:,2:3)+1;

greenmap = [0 0 0 % For RTK
   0 .1 0 
   0 .2 0 
   0 .3 0 
   0 .4 0 
   0 .5 0 
   0 .6 0 
   0 .7 0 
   0 .8 0 
   0 .9 0 
   0 1 0 ];
% greenmap(:,1) = greenmap(:,1)+1;
% greenmap(:,3) = greenmap(:,3)+1;

cyanmap = [0 0 0 % Mix of green and blue -> P53+RTK
   0 .1 .1 
   0 .2 .2 
   0 .3 .3 
   0 .4 .4 
   0 .5 .5 
   0 .6 .6 
   0 .7 .7 
   0 .8 .8 
   0 .9 .9 
   0 1 1 ];
% cyanmap(:,1) = cyanmap(:,1)+1;

yellowmap = [0 0 0 % Mix of green and red -> Rb+RTK
   .1 .1 0 
   .2 .2 0 
   .3 .3 0 
   .4 .4 0 
   .5 .5 0 
   .6 .6 0 
   .7 .7 0 
   .8 .8 0 
   .9 .9 0 
   1 1 0 ];
% yellowmap(:,3) = yellowmap(:,3)+1;

magentamap = [0 0 0 % Mix of blue and red -> Rb+P53
   .1 0 .1 
   .2 0 .2 
   .3 0 .3 
   .4 0 .4 
   .5 0 .5 
   .6 0 .6 
   .7 0 .7 
   .8 0 .8 
   .9 0 .9 
   1 0 1 ];
% magentamap(:,2) = magentamap(:,2)+1;



i = Size;

figure()
h = slice(poptot{i}(:,:,:),[],40,40);
colormap gray
set(h,'FaceColor','interp');
set(h,'EdgeColor','none');
axis equal
axis([0 80 0 80])
axis tight
axis off
caxis([0 2e5]);
title('Tumor slice','FontSize',18)
xlabel('x [mm]','FontSize',16)
ylabel('y [mm]','FontSize',16)
colorbar
set(gcf,'color','w');
% print(gcf,['Sim' num2str(k) '/SingleSlice.png'],'-dpng','-r600');


figpos = get(groot, 'ScreenSize');
ventana = figure('Name','Evolution of mutations','Units','pixels','OuterPosition',[0 0 figpos(4)/sqrt(2) figpos(4)]);
% colormap parula
subplot(3,3,1)
h1 = slice(pops{i}(:,:,:,1),[],40,40);
set(h1,'FaceColor','interp');
set(h1,'EdgeColor','none');
axis equal
axis([0 80 0 80])
caxis([0 max(max(pops{i}(:,:,40,1)))]);
title('Wild-type','FontSize',18)
xlabel('x [mm]','FontSize',16)
ylabel('y [mm]','FontSize',16)
colormap(gca,flipud('gray'))

subplot(3,3,2)
h2 = slice(pops{i}(:,:,:,2),[],40,40);
set(h2,'FaceColor','interp');
set(h2,'EdgeColor','none');
axis equal
axis([0 80 0 80])
caxis([0 max(max(pops{i}(:,:,40,2)))]);
title('P53','FontSize',18)
xlabel('x [mm]','FontSize',16)
ylabel('y [mm]','FontSize',16)
colormap(gca,bluemap)

subplot(3,3,3)
h3 = slice(pops{i}(:,:,:,3),[],40,40);
set(h3,'FaceColor','interp');
set(h3,'EdgeColor','none');
axis equal
axis([0 80 0 80])
caxis([0 max(max(pops{i}(:,:,40,3)))]);
title('Rb','FontSize',18)
xlabel('x [mm]','FontSize',16)
ylabel('y [mm]','FontSize',16)
colormap(gca,redmap)

subplot(3,3,4)
h4 = slice(pops{i}(:,:,:,4),[],40,40);
set(h4,'FaceColor','interp');
set(h4,'EdgeColor','none');
axis equal
axis([0 80 0 80])
caxis([0 max(max(pops{i}(:,:,40,4)))]);
title('P53+Rb','FontSize',18)
xlabel('x [mm]','FontSize',16)
ylabel('y [mm]','FontSize',16)
colormap(gca,magentamap)

subplot(3,3,5)
h5 = slice(pops{i}(:,:,:,5),[],40,40);
set(h5,'FaceColor','interp');
set(h5,'EdgeColor','none');
axis equal
axis([0 80 0 80])
caxis([0 max(max(pops{i}(:,:,40,5)))]);
title('RTK','FontSize',18)
xlabel('x [mm]','FontSize',16)
ylabel('y [mm]','FontSize',16)
colormap(gca,greenmap)

subplot(3,3,6)
h6 = slice(pops{i}(:,:,:,6),[],40,40);
set(h6,'FaceColor','interp');
set(h6,'EdgeColor','none');
axis equal
axis([0 80 0 80])
caxis([0 max(max(pops{i}(:,:,40,6)))]);
title('P53+RTK','FontSize',18)
xlabel('x [mm]','FontSize',16)
ylabel('y [mm]','FontSize',16)
colormap(gca,cyanmap)

subplot(3,3,7)
h7 = slice(pops{i}(:,:,:,7),[],40,40);
set(h7,'FaceColor','interp');
set(h7,'EdgeColor','none');
axis equal
axis([0 80 0 80])
caxis([0 max(max(pops{i}(:,:,40,7)))]);
xlabel('x [mm]','FontSize',16)
ylabel('y [mm]','FontSize',16)
title('Rb+RTK','FontSize',18)
colormap(gca,yellowmap)

subplot(3,3,8)
h8 = slice(pops{i}(:,:,:,8),[],40,40);
set(h8,'FaceColor','interp');
set(h8,'EdgeColor','none');
axis equal
axis([0 80 0 80])
caxis([0 max(max(pops{i}(:,:,40,8)))]);
xlabel('x [mm]','FontSize',16)
ylabel('y [mm]','FontSize',16)
title('P53+Rb+RTK','FontSize',18)
colormap(gca,flipud('gray'))

set(gcf,'color','w');
% print(gcf,['Sim' num2str(k) '/SingleSlice_AllMuts.png'],'-dpng','-r600');





%% RIM DISTANCE

K = 2e5;
threshold1 = 0.2*K;
threshold2 = 100;
mascP1 = cell(1,Size);
stats1 = cell(1,Size);
mascP2 = cell(1,Size);
stats2 = cell(1,Size);
rim = zeros(1,Size);


PX = 1; % mm
PY = 1; % mm
PZ = 1; % mm

m = 80;
n = 80;
p = 80;

[mallaX,mallaY,mallaZ] = meshgrid(PX/2 : PX : (PX*m-PX/2), ...
                                  PY/2 : PY : (PY*n-PY/2), ...
                                  PZ/2 : PZ : (PZ*p-PZ/2));

for i = 1:Size


    mascP1{i} = poptot{i}+nec{i} > threshold1;
    mascP2{i} = nec{i} > threshold1;

    VOLrim1(i) = sum(sum(sum(mascP1{i})));
    VOLrim2(i) = sum(sum(sum(mascP2{i})));

    rim(i) = 0.62*( (VOLrim1(i))^(1/3) - (VOLrim2(i))^(1/3) );

end


figure()
subplot(121)
imagesc(mascP1{Size}(:,:,40)-mascP2{Size}(:,:,40))
title('Contrast Enhancing Volume','Fontsize',18)
subplot(122)
imagesc(mascP2{Size}(:,:,40))
title('Necrotic Volume','Fontsize',18)
colormap gray

figure()
hold on
plot(VOLrim1)
plot(VOLrim2)
hold off

t = linspace(0,100,length(rim));
figure()
plot(t,rim,'LineWidth',2)
xlabel('Time [%]','FontSize',18)
ylabel('Rim width [mm]','FontSize',18)
set(gca,'FontSize',16)
set(gcf,'color','w');
% print(gcf,['Sim' num2str(k) '/RimWidth.png'],'-dpng','-r600');


% save(['Sim' num2str(k) '/RimDistance.mat'],'rim');

