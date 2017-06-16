%% prepare MATLAB workspace
clc;
clear all;
close all;

%% reading images and scanning information
folder = 'I:\Iron Loading Assessment\data\01\1.5T\';
folder_info = dir(folder);
Img = [];
if folder(end-4:end-1) == '3.0T'
    for i = 3:length(folder_info)
        if folder_info(i).name(1:4) == 'Scan'
            for k = 1:8
                NIFTI = load_untouch_nii([folder folder_info(i).name '\' num2str(k) '.nii']);
                Img(:,:,k,str2num(folder_info(i).name(6))) = double(NIFTI.img);
                clear NIFTI
            end
        end
    end
elseif folder(end-4:end-1) == '1.5T'
    folder_info_new = folder_info(1:2);
    for i = 3:length(folder_info)
        if folder_info(i).name(1:4) == 'Scan'
            folder_info_new(length(folder_info_new)+1) = folder_info(i);
        end
    end
    if length(folder_info_new) == 3
        s=1;
        for k = 1:8
            NIFTI = load_untouch_nii([folder 'Scan 1\' num2str(k) '.nii']);
            Img(:,:,k) = double(NIFTI.img);
            clear NIFTI
        end
    else
        str = {folder_info_new(3:end).name};
        [s,v] = listdlg('PromptString','Select a file:',...
            'SelectionMode','single',...
            'ListString',str);
        clear folder_info_new
        for k = 1:8
            NIFTI = load_untouch_nii([folder str{s} '\' num2str(k) '.nii']);
            Img(:,:,k) = double(NIFTI.img);
            clear NIFTI
        end
    end
end 
    
clear folder_info
TE_TR = csvread([folder 'TE_TR.csv']);

%%  display images
for k = 1:size(Img,4)
    figure('name',['Image of scan ' num2str(k)],...
        'units','normalized','outerposition',[0 0 1 1]);
    Img_max = max(max(max(Img(:,:,:,k))));
    for i=1:8
        subplot(2,4,i);
        imagesc(Img(:,:,i,k)/Img_max);
        colormap(gray);
        axis off;
    end
    clear Image_max
    pause(1)
    close all
end

%% image registration
Image = [];
hWaitBar = waitbar(0,'Image registration in progress');
count = 0;
for k = 1:size(Img,4)
    Image(:,:,1,k) = Img(:,:,1,k);
    for i = 2:8
        [Max(i,k),X(i,k),Y(i,k),Image(:,:,i,k)] = Mcurve_My(Img(:,:,1,k),Img(:,:,i,k),5,30);
        count = count+1;
        waitbar(count/(size(Img,4)*7));
    end
end
delete(hWaitBar);
clear count

%%  display registered images
for k = 1:size(Image,4)
    figure('name',['Registered image of scan ' num2str(k)],...
        'units','normalized','outerposition',[0 0 1 1]);
    Image_max = max(max(max(Image(:,:,:,k))));
    for i=1:8
        subplot(2,4,i);
        imagesc(Image(:,:,i,k)/Image_max);
        colormap(gray);
        axis off;
    end
    clear Image_max
    pause(1);
    close all;
end

%% GAC process
for k = 1:size(Image,4);
    %% initial contour
    nrow = size(Image,1);
    ncol = size(Image,2);
    Image_max = max(max(Image(:,:,1,k)));
    %{
    f = figure(1);
    set(f,'Position',[67 500 560 420]);
    imagesc(Image(:,:,1,k)/Image_max);
    colormap(gray);
    axis off;
    close(f)
    %}
    g = figure;
    set(g,'Position',[644 500 560 420]);
    imagesc(Image(:,:,1,k)/Image_max);
    colormap(gray);
    axis off;
    
    %% obtain ventricle center
    [centers, radii, metric] = imfindcircles(Image(:,:,1,k),[15 30]);  % finding circles with radius between 15-30
    if size(centers,1) == 0
        title('Auto ventricle detection fail, please click on the center of blood pool');
        [x_0 y_0] = ginput(g);
        hold on
        scatter(x_0,y_0,'r*');
        title('Auto ventricle detection fail, please click on the epicardium');
        [x_end y_end] = ginput(g);
        radius_outer = sqrt((x_end-x_0).^2+(y_end-y_0).^2); %outer radius
        viscircles([x_0 y_0], radius_outer,'EdgeColor','r','LineWidth',1);
        title('Auto ventricle detection fail, please click on the endocardium');
        [x_epi y_epi] = ginput(g);
        radius_inner = sqrt((x_epi-x_0).^2+(y_epi-y_0).^2); %inner radius
        viscircles([x_0 y_0], radius_inner,'EdgeColor','r','LineWidth',1);
        pause(1.5)
        close all
    elseif size(centers,1) == 1
        % Construct a questdlg to check for correct detection
        title('Finish auto ventricle detection');
        x_0 = centers(1);
        y_0 = centers(2);
        radius_outer = radii;
        radius_inner = radii-7;
        hold on
        viscircles([x_0 y_0;x_0 y_0], [radius_outer;radius_inner],'EdgeColor','r','LineWidth',1);
        scatter(x_0,y_0,'r*');
        pause(1.5)
        choice = questdlg('Is the detection correct?', ...
            '', ...
            'Yes','No','Yes');
        % Handle response
        switch choice
            case 'Yes'
                pause(1.5)
                close all;
            case 'No'
                hold off
                imagesc(Image(:,:,1,k)/Image_max);
                axis off;
                hold on
                title('Auto ventricle detection fail, please click on the center of blood pool');
                [x_0 y_0] = ginput(g);
                scatter(x_0,y_0,'r*');
                title('Auto ventricle detection fail, please click on the epicardium');
                [x_end y_end] = ginput(g);
                radius_outer = sqrt((x_end-x_0).^2+(y_end-y_0).^2); %outer radius
                viscircles([x_0 y_0], radius_outer,'EdgeColor','r','LineWidth',1);
                title('Auto ventricle detection fail, please click on the endocardium');
                [x_epi y_epi] = ginput(g);
                radius_inner = sqrt((x_epi-x_0).^2+(y_epi-y_0).^2); %inner radius
                hold on
                viscircles([x_0 y_0], radius_inner,'EdgeColor','r','LineWidth',1);
                pause(1.5)
                close all
        end
    elseif size(centers,1) > 1
        % Construct a questdlg to check for correct detection
        title('Finish auto ventricle detection');
        x_0 = mean(centers(1:2,1));
        y_0 = mean(centers(1:2,2));
        radius_outer = radii(1);
        radius_inner = radii(1) - 7;
        hold on
        viscircles([x_0 y_0;x_0 y_0], [radius_outer;radius_inner],'EdgeColor','r','LineWidth',1);
        scatter(x_0,y_0,'r*');
        pause(1.5)
        choice = questdlg('Is the detection correct?', ...
            '', ...
            'Yes','No','Yes');
        % Handle response
        switch choice
            case 'Yes'
                pause(1.5)
                close all;
            case 'No'
                hold off
                imagesc(Image(:,:,1,k)/Image_max);
                axis off;
                hold on
                title('Auto ventricle detection fail, please click on the center of blood pool');
                [x_0 y_0] = ginput(g);
                scatter(x_0,y_0,'r*');
                title('Auto ventricle detection fail, please click on the epicardium');
                [x_end y_end] = ginput(g);
                radius_outer = sqrt((x_end-x_0).^2+(y_end-y_0).^2); %outer radius
                viscircles([x_0 y_0], radius_outer,'EdgeColor','r','LineWidth',1);
                title('Auto ventricle detection fail, please click on the endocardium');
                [x_epi y_epi] = ginput(g);
                radius_inner = sqrt((x_epi-x_0).^2+(y_epi-y_0).^2); %inner radius
                hold on
                viscircles([x_0 y_0], radius_inner,'EdgeColor','r','LineWidth',1);
                pause(1.5)
                close all
        end
    end
    initialLSF_outer = sdf2circle(nrow,ncol,y_0,x_0,radius_outer);
    initialLSF_inner = sdf2circle(nrow,ncol,y_0,x_0,radius_inner);
    
    %% plot outer and inner image
    u_outer = initialLSF_outer;
    u_inner = initialLSF_inner;
    figure('name',['Initial contour of scan ' num2str(k)],...
        'units','normalized','outerposition',[0 0 1 1]);
    for i=1:8
        subplot(2,4,i);
        imagesc(Image(:,:,i,k));
        colormap(gray);
        axis off;
        hold on;
        scatter(x_0,y_0,'r*');
        contour(u_outer,[0,0],'r');
        contour(u_inner,[0,0],'r');
    end
    pause(1)
    close all
    
    %% image smoothen
    sigma = 1.2;    % scale parameter in Gaussian kernel
    G = fspecial('gaussian',15,sigma);
    Imgsmooth = cell(1,8);
    f = cell(1,8);
    g = cell(1,8);
    for i=1:8
        Imgsmooth{i} = conv2(double(Image(:,:,i,k)),double(G),'same');  %image smoothen with Gaussian kernel
        [Ix Iy] = gradient(Imgsmooth{i}); 
        f{i} = Ix.^2+Iy.^2;
        g{i} = 1./(1+f{i});
    end
    
    clc
    [grow gcol] = size(g{1});
    gmax = g{1};
    gmin = g{1};
    fm = f{1};
    S = grow * gcol;
    for i=2:8
        for j = 1:S
            gmin(j) = min(gmin(j),g{i}(j));    %取最小值
            gmax(j) = max(gmax(j),g{i}(j));
        end
    end
    
    %% GAC iteration
    alf=0.5;  % coefficient of the pressure/balloon force 
    timestep=.1;  
    iterNum_evo = 10; % number of iterations in level set evolution before re-initialization
    iterNum_ri = 10; % number of iterations in re-initialization  
    figure('name',['GAC iteration of scan ' num2str(k)],...
        'units','normalized','outerposition',[0 0 1 1]);
    hWaitBar = waitbar(0,'GAC iteration in progress');
    count = 0;
    for n=1:100
        u_outer = VEVOLUTION_GAC(u_outer, gmax, alf, timestep, iterNum_evo);
        u_outer = reinit_SD(u_outer, 1, 1, .5, iterNum_ri);
        u_inner = VEVOLUTION_GAC(u_inner, gmin, alf, timestep, iterNum_evo);
        u_inner = reinit_SD(u_inner,1,1,.5,iterNum_ri);      %保证曲线满足符号函数
        pause(0.05);
        if mod(n,10)==0
            for i=1:8
                subplot(2,4,i);
                imagesc(Image(:,:,i,k));
                colormap(gray);
                axis off;
                hold on;
                [c2,h2] = contour(u_outer,[0 0],'r');
                [c3,h3] = contour(u_inner,[0 0],'r');
                scatter(x_0,y_0,'r*');
                iterNum=[num2str(n*iterNum_evo), ' iterations'];
                title(iterNum);
            end
            % iteration plotting
            figure('name',['GAC iteration ' num2str(n*iterNum_evo) ' of 1st echo time image']);
            imagesc(Image(:,:,1,k));
            colormap(gray);
            axis off;
            hold on;
            [c2,h2] = contour(u_outer,[0 0],'r','LineWidth',1.5);
            [c3,h3] = contour(u_inner,[0 0],'r','LineWidth',1.5);
            scatter(x_0,y_0,'r*');
            iterNum=[num2str(n*iterNum_evo), ' iterations'];
            title(iterNum);
            saveas(gcf,[folder 'Scan ' num2str(s) '\itr_result_' num2str(n*iterNum_evo) '.jpg']);
            close(gcf);
        end
        count = count+1;
        waitbar(count/100);
    end
    pause(1);
    close all
    delete(hWaitBar);
    clear count
    if folder(end-4:end-1) == '3.0T'
        save([folder 'Scan ' num2str(k) '\matfile']);
    else
        save([folder 'Scan ' num2str(s) '\matfile']);
    end
    
    %% display segmententation image
    for i = 1:8
        figure('name',['GAC iteration result of scan ' num2str(k)],...
            'position',[289 173 size(Image(:,:,i,k),1)*2 size(Image(:,:,i,k),2)*2]);
        imagesc(Image(:,:,i,k));
        colormap(gray);
        axis off;
        hold on;
        [c2,h2] = contour(u_outer,[0 0],'r');
        [c3,h3] = contour(u_inner,[0 0],'r');
        scatter(x_0,y_0,'r*');
        if folder(end-4:end-1) == '1.5T'
            saveas(gcf,[folder 'Scan ' num2str(s) '\seg_result_' num2str(i) '.jpg']);
        else
            saveas(gcf,[folder 'Scan ' num2str(k) '\seg_result_' num2str(i) '.jpg']);
        end
        close all;
    end
end