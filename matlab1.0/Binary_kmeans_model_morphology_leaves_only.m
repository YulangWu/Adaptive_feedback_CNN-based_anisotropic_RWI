clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%
% Strategy: maximize the number of clusters at the first and second steps
% Randomly choose the final number of clusters to output model
%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0. set experiment parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%0.1 Set parameters and input starting model
%Marmousi
nz=256;
nx=256*3;
%
nz=192;
nx=192*3;
% water depth of the model
water_depth = 6; 

%0.2 k-means parameters
max_iteration = 300;
kmean_num = 2; %binary segmentation
% num_label = 100; %(abandon) this number should be large enough for segmentation 
max_num_group = 2; %Given the num_label, keep the first k biggest groups (HERE IT SHOULD BE 2s)
tree_depth = 14; %maximum depth of tree
threshold = floor(0.01*nz*nx); % the children will be found when sum(sign(map)) > threshold
dir = 'Marmousi'; %'right_marmousi';

filename =  '0th_mig_Marmousi';
% filename =  '1th_mig_Marmousi';
% filename =  '2th_mig_Marmousi';
% filename =  '3th_mig_Marmousi';

disp(filename);

vp_mig=dlmread([dir '/' filename '.dat']);

matfilename = [dir '_' filename '_binary_tree' '_threshold' num2str(threshold) '.dat'];
velocity_dir = 'velocity';
system(['mkdir ' velocity_dir])

filename2 = '0th_true_Marmousi'; %'mig_10Marmousi256x768'; %'true_Marmousi256x768';
vp_true=dlmread([dir '/' filename2 '.dat']);

%0.3 smoothing parameters
smooth_num = 200; %smooth number for smoothing clustered velocity
smooth_filter_size = 3;

%0.4 create new model
std_factor = 2; %v = v_mean _ v_std*factor
num_new_model = 32; %number of new velocity
%min and max velocity in each new generated model
min_val = 1.4;
max_val = 5.0;
%min and max r2 of new model with respect to starting model (vp_mig)
min_r2 = 0.5;
max_r2 = 1.0;

disp(['tree_depth = ' num2str(tree_depth)])
disp(['threshold = ' num2str(threshold)])
disp(['num_new_model = ' num2str(num_new_model)])
disp(['std_factor = ' num2str(std_factor)])
disp(['min_val = ' num2str(min_val)])
disp(['max_val = ' num2str(max_val)])
disp(['min_r2 = ' num2str(min_r2)])
disp(['max_r2 = ' num2str(max_r2)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. binary-classification of model using k-means
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist(matfilename, 'file') == 0 %if file exists do not repeatedly create
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % vp_mig = dlmread('new_FWI_model199.dat');
    % nz=128;
    % nx=256;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %1.1 Input starting model used for LSRTM
    vp_mig = reshape(vp_mig,nz,nx);
    vp_mig(1:water_depth,:) = 1.5;

    disp('Input a 2-D model output a 2-D model only')
    [nz nx] = size(vp_mig);
    vp_mig = reshape(vp_mig,nz*nx,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %1.2 Define tree structure
    tree = struct('lchild',0,'rchild',0,'value',0);

    root = tree;
    root.value = vp_mig;
    lchild = tree;
    rchild = tree;
    %1.3 Get binary segmentation for root, as no blank area, so the null_num = 0
    [lchild.value, rchild.value]= get_seperated_two_models_for_tree_structure_only(root.value,kmean_num, 0, nz, nx, max_num_group);
    root.lchild = lchild;
    root.rchild = rchild;
    %1.4 Get binary segmentation for children, as blank area exists, so the null_num = 1
    root.lchild = recurr_kmean_decomposition_for_tree_structure_only(root.lchild, kmean_num, 1, nz, nx, max_num_group, tree_depth, threshold);
    root.rchild = recurr_kmean_decomposition_for_tree_structure_only(root.rchild, kmean_num, 1, nz, nx, max_num_group, tree_depth, threshold);
    
    %1.5 save the leaves as array to disk
    segmented_layer = get_model_from_tree_at_leaves(root, nz, nx, tree_depth);
    %1.6 Get feature maps with unique feature (Not 0's)
    num_feature_layers = 0;
    for i = 1:length(segmented_layer(1,:))
        if sum(segmented_layer(:,i)) ~= 0
            num_feature_layers = num_feature_layers + 1;
%         else
%             disp(i)
%             break;
        end
    end
    %1.7 Crop feature layers by deleting 0's feature map
    segmented_layer = segmented_layer(:,1:num_feature_layers);
    
    fid=fopen([matfilename(1:length(matfilename)-4) '.dat'],'wt');
    fprintf(fid,'%15.12f',segmented_layer);
    fclose(fid);
    
%     try
%         disp('Try to save mat file.....');
%         save([matfilename(1:length(matfilename)-4) '.mat']);
%     catch ME
%         disp('mat file is too large to be saved to disk! Save segmented_layer instead!');
%         fid=fopen(matfilename,'wt');
%         fprintf(fid,'%15.12f',segmented_layer);
%         fclose(fid);
%     end
end
%1.6 load the saved array to memory

% try
%     disp('Try to load mat file.....');
%     load([matfilename(1:length(matfilename)-4) '.mat']);
% catch ME
%     disp('mat file does not exist! Load segmented_layer instead!');
%     segmented_layer = dlmread([matfilename(1:length(matfilename)-4) '.dat']);
%     segmented_layer = reshape(segmented_layer,nz*nx,floor(length(segmented_layer)/(nz*nx)));
% end 
segmented_layer = dlmread([matfilename(1:length(matfilename)-4) '.dat']);
segmented_layer = reshape(segmented_layer,nz*nx,floor(length(segmented_layer)/(nz*nx)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. load segmented models and plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. check whether  segmentation is complete (reconstruction is correct)
cmap = combine_kmeans_layers_boundary( segmented_layer);

% cmap = vp_mig; %only if the difference can be visually ignored!!!
disp(['error=' num2str(sum(cmap-vp_mig))])
if sum(cmap-vp_mig) < 30
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %3. plot k-means result
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %3.1 plot the result
% % %     close all;
% % %     min_c = 0;
% % %     max_c = 4.5;
% % %     show2dlayer(cmap, 100, 1, 1, nz, nx, min_c, max_c);
% % % %     show2dlayer( segmented_layer(:,19:82), 1, 8, 8, nz, nx, min_c, max_c);
% % % %     show2dlayer( segmented_layer(:,83:146), 2, 8, 8, nz, nx, min_c, max_c);
% % %     
% % %     %show puzzled feature map:
% % % %     %Marmousi with 0.1 threshold
% % %     puzzle_layer1 = get_puzzle_layer(segmented_layer(:,155-64-64+1:155-64),6,6,nz,nx, 0, 5.0,123);
% % %     puzzle_layer2 = get_puzzle_layer(segmented_layer(:,155-64+1:155),6,6,nz,nx, 0, 5.0,321);
% % %     puzzle_layer3 = get_puzzle_layer(segmented_layer(:,1:28),7,4,nz,nx, 0, 5.0,321);

%     puzzle_layer2 = get_puzzle_layer(segmented_layer(:,1:num_feature_layer),7,5,nz,nx, 0, 5,123);

%     %BP with 0.1 threshold
%     puzzle_layer1 = get_puzzle_layer(segmented_layer(:,3:74),9,8,nz,nx, -1, 1,123);
%     puzzle_layer2 = get_puzzle_layer(segmented_layer(:,74:146),9,8,nz,nx, -1, 1,321);
    
%     puzzle_layer2 = get_puzzle_layer(segmented_layer(:,83:146),8,8,nz,nx, -1, 1,123);

% 
    %3.2 (abandon)plot original model with boundaries 
    % for combination, do not use the functions
    % (extract_individual_kmeans_layers_boundary.m)
    % and (combine_kmeans_layers.m) 
    % Instead, use function (extract_whole_model_boundary.m)
    % segmented_layer_boundary = extract_individual_kmeans_layers_boundary(sign(segmented_layer), nz, nx, threshold/2);
    % segmented_layer_combine_boundary = combine_kmeans_layers(segmented_layer_boundary)*1000;

% % %     %3.3 (adopted) plot original model with boundaries 
% % %     [boundary_map,color_map] = extract_whole_model_boundary(segmented_layer, nz, nx);
% % %     min_c = 1.5;
% % %     max_c = 5;
% % %     boundary_map = reshape(boundary_map,nz*nx,1);
% % %     show2dlayer(cmap+boundary_map*max_c, 101, 1, 1, nz, nx, 1.5, max_c);
% % %     set(gca,'XAxisLocation','top');
% % %     set(gca, 'XTick', [1 nx/3 nx/3*2 nx])            
% % %     set(gca,'XTickLabel',{'0.00','2.88','5.76','8.64'}) 
% % %     set(gca, 'YTick', [1 nx/6 nx/3])          
% % %     set(gca,'YTickLabel',{'0.00','1.44','2.88'})  
% % %     axis equal; 
% % %     xlabel('Distance (km)'); 
% % %     xlim([1 nx]);ylim([1 nz]);
% % %     ylabel('Depth (km)');
% % %     c = colorbar();
% % %     colorTitleHandle = get(c,'YLabel');
% % %     titleString = 'Velocity (km/s)';
% % %     set(colorTitleHandle ,'String',titleString);
% % %     set(c,'YTick',[1.5,2.5,3.5,4.5])
% % %     set(c,'YTickLabels',{'1.5','2.5','3.5','4.5'})
% % %     set(c,'Position',[0.841765873015872, 0.277777777777778, ...
% % %     0.0272817460317471, 0.469246031746029]) %get from the code of figure
% % % 
% % %     show2dlayer(color_map, 102, 1, 1, nz, nx, -5, length(segmented_layer(1,:))+10);
% % %     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %4. Get each group's mean and std
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %4.3 get mean and std for each feature map 
    % Notice that matlab-built-in mean and variance is not proper since
    % feature map has blank area which should be excluded in calculation
    hash_map_value = get_mean_std( segmented_layer);
%     figure(1000);plot(hash_map_value(:,1),hash_map_value(:,2),'o');
%     xlabel('Mean velocity (km/s)');ylabel('std velocity (km/s)');
%     set(gca,'XAxisLocation','top');
%     set(gca, 'XTick', [1.5 2.5 3.5 4.5])            
%     set(gca,'XTickLabel',{'1.5','2.5','3.5','4.5'}) 
    
    fid=fopen(['new_model_mean_std_.dat'],'wt');
    fprintf(fid,'%9.6f',hash_map_value);
    fclose(fid);
    %for read and plot a = ['new_model_mean_std_.dat']; a = reshape(a,34,2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %5. Create random Marmousi model
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %get random value for each group which satisfies normal distribution
    r2_value = zeros(num_new_model,2);
    figure(10000)
    
    vp_new_set = zeros(nz*nx,num_new_model);
    k = 1;
    while k <= num_new_model
        hash_map_value_randn = normal_dist(hash_map_value(:,1),hash_map_value(:,2)*std_factor,min_val,max_val);
        vp_new = zeros(nz*nx,1);
        for i = 1:length(segmented_layer(1,:))
                vp_new = vp_new + sign(segmented_layer(:,i))*hash_map_value_randn(i);
        end
        r2_value(k,1) = R2(vp_mig,vp_new);
        r2_value(k,2) = R2(vp_true,vp_new);
        if r2_value(k,1) <min_r2 || r2_value(k,1) > max_r2
            continue; %skip this generated model
        else
            fid=fopen([velocity_dir '/' 'new_model' num2str(k) 'r2_sm' num2str(floor(r2_value(k,1)*100)) 'r2_tr' num2str(floor(r2_value(k,2)*100)) '.dat'],'wt');
            fprintf(fid,'%9.6f',vp_new);
            fclose(fid);

            vp_new_set(:,k) = vp_new;
%             show2d(vp_new,nz,nx);caxis([1.5 4.5]);title([num2str(k) 'th model with R2 = ' num2str(r2_value(k,1))]);
            drawnow;
            pause(0.1);
            k = k + 1;
            
        end
    end
    
    %show puzzled feature map:
    puzzle_layer1 = get_puzzle_layer(vp_new_set,6,4,nz,nx, 1.5, 5,12345);
    
    %plot and save r2 value
    figure(10001)
    plot(r2_value);ylim([0.85 1]);
    xlabel('Number of training model');ylabel('Correlation coefficient')
    fid=fopen(['new_model_r2_.dat'],'wt');
    fprintf(fid,'%9.6f',r2_value);
    fclose(fid);
    
else  
    disp('Segmentation failed! Original image is not reconstructed!');
    %below is for debugging
    %show2dlayer( segmented_layer, 2000, 10, 10, nz, nx, min_c, max_c);
    figure(102)
    subplot(2,1,1)
    imagesc(reshape(color_map,nz,nx))
    caxis([-5  length(conventional_segmented_layer(1,:))+10])
    set(gca,'XAxisLocation','top');
    set(gca, 'XTick', [1 nx/3 nx/3*2 nx])            
    set(gca,'XTickLabel',{'0.00','2.88','5.76','8.64'}) 
    set(gca, 'YTick', [1 nx/6 nx/3])          
    set(gca,'YTickLabel',{'0.00','1.44','2.88'})  
    axis equal; 
    xlabel('Distance (km)'); 
    xlim([1 nx]);ylim([1 nz]);
    ylabel('Depth (km)');
    c = colorbar();
    colorTitleHandle = get(c,'XLabel');
    titleString = 'Index value';
    set(colorTitleHandle ,'String',titleString);
    set(c,'Position',[0.23363095238095 0.550958994708988 0.566964285714288 0.0300595238095237]) %get from the code of figure
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    layers_depth1 = root.value;
    layers_depth2 = zeros(nz*nx,2);
    layers_depth3 = zeros(nz*nx,4);
    layers_depth4 = zeros(nz*nx,8);
    layers_depth5 = zeros(nz*nx,16);
    layers_depth6 = zeros(nz*nx,32);
    layers_depth7 = zeros(nz*nx,64);

    layers_depth2 = get_model_from_tree_at_depth(root,2,layers_depth2,0);
    layers_depth3 = get_model_from_tree_at_depth(root,3,layers_depth3,0);
    layers_depth4 = get_model_from_tree_at_depth(root,4,layers_depth4,0);
    layers_depth5 = get_model_from_tree_at_depth(root,5,layers_depth5,0);
    layers_depth6 = get_model_from_tree_at_depth(root,6,layers_depth6,0);
    layers_depth7 = get_model_from_tree_at_depth(root,7,layers_depth7,0);
    min_c = 0;
    max_c = 4.5;
    show2d(layers_depth1, nz, nx);caxis([min_c max_c]);
    show2dlayer( layers_depth2, 2, 1, 2, nz, nx, min_c, max_c);
    show2dlayer( layers_depth3, 3, 2, 2, nz, nx, min_c, max_c);
    show2dlayer( layers_depth4, 4, 2, 4, nz, nx, min_c, max_c);
    show2dlayer( layers_depth5, 5, 4, 4, nz, nx, min_c, max_c);
    show2dlayer( layers_depth6, 6, 4, 8, nz, nx, min_c, max_c);
    show2dlayer( layers_depth7, 7, 8, 8, nz, nx, min_c, max_c);
    
    cmap1 = combine_kmeans_layers_boundary( layers_depth1);
    cmap2 = combine_kmeans_layers_boundary( layers_depth2);
    cmap3 = combine_kmeans_layers_boundary( layers_depth3);
    cmap4 = combine_kmeans_layers_boundary( layers_depth4);
    cmap5 = combine_kmeans_layers_boundary( layers_depth5);
    cmap6 = combine_kmeans_layers_boundary( layers_depth6);
    cmap7 = combine_kmeans_layers_boundary( layers_depth7);
    show2dlayer( cmap1 - vp_mig, 22, 1, 1, nz, nx, -max_c, max_c);
    show2dlayer( cmap2 - vp_mig, 23, 1, 1, nz, nx, -max_c, max_c);
    show2dlayer( cmap3 - vp_mig, 24, 1, 1, nz, nx, -max_c, max_c);
    show2dlayer( cmap4 - vp_mig, 25, 1, 1, nz, nx, -max_c, max_c);
    show2dlayer( cmap5 - vp_mig, 26, 1, 1, nz, nx, -max_c, max_c);
    show2dlayer( cmap6 - vp_mig, 27, 1, 1, nz, nx, -max_c, max_c);
    show2dlayer( cmap7 - vp_mig, 28, 1, 1, nz, nx, -max_c, max_c);
    
  
end