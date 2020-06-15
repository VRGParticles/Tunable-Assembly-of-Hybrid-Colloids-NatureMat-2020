%% Analysis of Brick (cmm) vs. Amorphous particle assemblies.
%     Copyright (C) 2020 Veronica Grebe, Weck Group at New York University.
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>

disp('Brick_Amorphous.m Copyright (C) 2020 Veronica Grebe, Weck Group at New York University')
disp('This program comes with ABSOLUTELY NO WARRANTY.')
disp('This is free software released under GNU GPL-3.0,')
disp('and you are welcome to redistribute it under certain conditions;')
disp('see <http://www.gnu.org/licenses/>')

%   Use: Label an optical micrograph into brick (cmm) symmetry, amorphous, shared and disorder particle assemblies.

%   Author Contact: Veronica Grebe vrg234@nyu.edu, New York University.

%   Function calls: This script calls the function 'bpass.pro' during image segmentation and was released under the fokkowibg lince:
%       This function, bpass.pro, is copyright:
%       John C. Crocker and David G. Grier.  It should be considered 'freeware'- and may be
%       distributed freely in its original form when properly attributed.
%       https://es.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/42573/versions/1/previews/bpass.first_neighbor/index.html

%%  Input parameters **EDIT HERE**
clearvars %clear any current variables in the workspace.

%file location and type
file_path = '/where/your/files/lives/'; %where your file(s) to be analyzed are located.
full_file_path = fullfile(file_path, '*.jpg'); % Change to whatever pattern your micrograph is.  Likely .jpg or .tif.  For one file, copy file name.  For multiple in the file_path, use a wildcard *.filetype.


%Image Preprocessing.
imgauss_param = 3; %gaussian filter parameter;
bwareaopen_param = 30; %how large must the area be to keep the blob.
bwareafilt_param = [60 150]; %lower and upper area bound of kept blobs.
bandpass_bpass = 1; %bandpass for bpass function
bandpass_strel = 2; %bandpass for oval/struture element

%distance cutoffs of centers of each particle.
% *** Note that in the script,
%constants are added to these numbers for the analysis of the images.
%These constants may need to be changed depending on the input images.***
distmin_dis = 0; %the min distance between particle centers for disordered particles
distmax_dis = 30; %the max distance between particle centers for disordered particles
distmin_parallel = 0; %the min distance between particle centers  for parallel particles
distmax_parallel = 25; %the max distance between particle centers for parallel particles

%angle cutoffs.
apara = 15; %parallel bounds, + or - about 0

%blob property cutoffs
malmin = 7; %minimum length a particle must be (in major axis)
areamin = 70; %minimum area a particle must have
perimin = 15; %minimum perimeter a particle must have

%other parameters

%CROPPING
crop = 0; %set 1 for yes, and to any other number for no.

%for cropping, sets pixel distance.
x1 = 0; %top left
y1 = 0; %bottom left
x2 = 200; %top right
y2 = 200; %bottom right

%% File Writing Setup **DO NOT EDI BELOW THIS LINE**

%Create directories to set files.
cd(file_path);
mkdir Image_Files; %folder for Images
mkdir Excel_Files; %folder for Excel files

% Get a list of all files in the folder with the inputfile name/pattern.
file_list = dir(full_file_path);

%% Begin Loop
for a = 1 : length(file_list)
    %read in each file in order it appears in the list.
    
    file_name_base = file_list(a).name;
    file_name_full = fullfile(file_path, file_name_base);
    fprintf(1, 'Now processing %s\n', file_name_base);
    
    %Read in the image to the workspace
    Im_original = imread(file_name_full);
    
    %crop image if user set crop to 1.
    if crop == 1
        Image_cropped = imcrop(Im_original, [x1 y1 x2 y2]);
        Image_processed = Image_cropped;
    else
        Image_processed = Im_original;
    end
    
    %rgb check, depends on image.
    rgb_check = size(Im_original,3);
    
    if rgb_check ==1
        rgbc = 0; %the image is greyscale
    else
        rgbc = 1; %the image is rgb
    end
    
    %% Image segmentation: binary conversion
    disp('Segmenting Image')
    
    if rgbc == 1
        Image_processed = rgb2gray(Image_processed); %convert to gray if the image is rgb
    end
    
    %apply Gaussian filter to blur image
    Image_processed = imgaussfilt(Image_processed,imgauss_param);
    
    %bandpass filter applied to image using function call bpass.m
    Image_processed = uint16(bpass(Image_processed,0,bandpass_bpass));
    
    %adjust image intensity to further segment out particles.
    Im_adjust = imadjust(Image_processed);
    
    greythresh_thresh = graythresh(Im_adjust); %find global image threshold
    
    imbw = Im_adjust > greythresh_thresh; %initial segmented, binary image
    
    %% Image segmentation: segment binary image
    
    se = strel('disk',bandpass_strel); %defines structure element to filter image based on shape.
    
    Image_opened_se = imopen(imbw>0,se); %morphological opening of image using structure element se.
    
    Image_opened_area = bwareaopen(Image_opened_se, bwareaopen_param); %remove initially small features
    
    Image_opened_area_filter = bwareafilt(Image_opened_area, bwareafilt_param); %remove small and large features outside of the accepted area range.
    
    %% Obtain blob properties
    
    [connected_binary,number_of_particles] = bwlabel(Image_opened_area_filter); %connect adjacent white pixels in segmented binary image
    
    stats = regionprops(connected_binary,'Area', 'Centroid', 'MajorAxisLength', 'Orientation', 'Perimeter'); %obtain defined properties for the connected regions of connected_binary
    %% Preallocation of arrays
    
    Amorphous_Total = zeros(number_of_particles,1); %Total count of amorphous particles
    area = zeros(number_of_particles,1); %area from regionprops
    Brick_Total = zeros(number_of_particles,1); %Total count of Brick particles
    cns = zeros(number_of_particles,2); %centroid (x,y) from regionprops
    Dif_Ang = zeros(number_of_particles); %orientation angle difference
    Disorder_Total = zeros(number_of_particles,1); %Total count of disordered particles
    mal = zeros(number_of_particles,1); %major axis length from regionprops
    Neighbors_Brick=  zeros(number_of_particles,1); %amount of brick neighbors
    ori = zeros(number_of_particles,1); %orientation from regionprops
    Parallel_Kept_0= zeros(number_of_particles,1); %count of particles neighboring para particles with 2 para neighbors.
    peri = zeros(number_of_particles,1); %perimeter from regionprops
    Perp_Kept = zeros(number_of_particles,1); %particles that are kept meeting perp. criteria
    Shared_Total = zeros(number_of_particles,1); %Total count of shared particles
    
    %% Index Region Stats
    %obtain centroids, orientation, major axis length, area and perimeter of
    %connected regions determined by 'stats' above
    for particle = 1:number_of_particles
        cns(particle,:) = stats(particle).Centroid;
        ori(particle,:) = stats(particle).Orientation;
        mal(particle,:) = stats(particle).MajorAxisLength;
        area(particle,:) = stats(particle).Area;
        peri(particle,:) = stats(particle).Perimeter;
    end
    
    %% Calculate distance and angle difference for all particles
    disp('Calculating Parameters')
    
    dist = squareform(pdist(cell2mat({stats(:).Centroid}'))); %computes center-to-center distance from all particles to all other particles
    
    [first_neighbor, second_neighbor] = find(dist < distmax_parallel + 30);   %determine which particles are neighboring, with some tolerance
    
    %Sweep through all neighboring particles and only calculate their
    %orientation angle difference if they meet any of the threshold
    %conditions for perimeter, major axis length, or area.
    
    for particle = 1:length(first_neighbor)
        particle1=first_neighbor(particle);
        particle2=second_neighbor(particle); 
        if mal(particle1) > malmin || area(particle1) > areamin || peri(particle1) > perimin %ensure the particle large enough in at least one of these properties
            
            Dif_Ang(particle1, particle2) = abs(stats(particle1).Orientation - stats(particle2).Orientation); %orientation angle difference
        end
        
    end
    
    
    %%  Parallel (Brick) particles
    
    AllowedDistance_Parallel = (distmin_parallel < dist & dist < distmax_parallel);  %neighboring particles that meet distance criteria
    
    AllowedAngles_Parallel =  (apara > Dif_Ang & Dif_Ang > -apara) ; %neighboring particles that meet angle criteria for parallel
    
    KeptParticlesParam_Parallel = (AllowedDistance_Parallel == AllowedAngles_Parallel) & (AllowedDistance_Parallel ==1); %determine if particles can be accepted by parallel if they meet both thresholds
    
    Parallel_Sum=sum(KeptParticlesParam_Parallel,2); % move to one column: how many neighbors each particle meets the above para criteria for.
    
    [first_neighbor, second_neighbor] = find(KeptParticlesParam_Parallel == 1); %find which particles meet the criteria.
    
    %Count for each particle1, which has at least 1 parallel neighbor,
    %count how many particle2 have at least 2 neighbors, if particle1 and
    %particle2 are parallel neighbors.
    
    for particle = 1:length(first_neighbor)
        particle1 = first_neighbor(particle);
        particle2 = second_neighbor(particle);
        
        if Parallel_Sum(particle1) >= 1 && Parallel_Sum(particle2) >= 2 && abs(Dif_Ang(particle1, particle2)) < apara
            Parallel_Kept_0(particle1) = Parallel_Kept_0(particle1) + 1; %count of how many particle pairs meet the criteria
        end
    end
    %If the particle meets the critera of having at least 2 para neighbors
    %that have at least 2 para neighbors, it is counted as a brick
    %particle (They are bundles of 3 para particles).
    Brick_Kept = Parallel_Kept_0 >= 2;
    
    
    %% Disorder particles
    disp('Finding Shapes')
    
    AllowedDistance_Disorder = (distmin_dis < dist & dist < distmax_dis);  %neighboring particles that meet distance criteria
    
    KeptParticlesParam_Disorder = (AllowedDistance_Disorder ==1); %determine if particles meet the distance threshold
    
    Disorder_Sum=sum(KeptParticlesParam_Disorder,2); % move to one column: how many neighbors each particle meets the above  disorder criteria.
    
    %Only mark particles as disordered if they have exactly 6 neighbors.
    Disorder_Kept = (Disorder_Sum ==6);
    
    %Mark all neighboring particles of particle1 (member of Disorder_Kept)
    %as disordered, also.  This functionally counts all particles that are
    %members of a set of 7 particles, and can overlap Brick_Kept.
    
    [first_neighbor, second_neighbor] = find(KeptParticlesParam_Disorder == 1);
    
    for particle = 1:length(first_neighbor)
        particle1 = first_neighbor(particle);
        particle2 = second_neighbor(particle);
        
        if Disorder_Sum(particle1) == 6
            Disorder_Kept(particle2) = 1;
        end
    end
    
    %% Shared Particles
    
    clear first_neighbor second_neighbor
    [first_neighbor,second_neighbor] = find(dist < distmax_parallel); %find all particles that meet distance threshold
    
    %for all neighbor pairs that meet distance threshold, find where both
    %particles are para AND accepted as brick.  Then, add one to particle1's
    %amount of brick neighbors
    for particle = 1:length(first_neighbor)
        particle1 = first_neighbor(particle);
        particle2 = second_neighbor(particle);
        
        if  Brick_Kept(particle1) == 1 && Brick_Kept(particle2) == 1 &&  Dif_Ang(particle1, particle2) < apara
            Neighbors_Brick(particle1) = Neighbors_Brick(particle1) + 1;
        end
    end
    
    %if a brick particle has at least 6 neighbors, but fewer than 3 or less brick neighbors,
    % mark the particle as shared (amorphous neighbors = [3, 6]). Thus,
    % Brick particles must have 5 or fewer total neighbors, OR have at
    % least 6 neighbors, with at least 4 being brick.
    Shared_Kept =  (Brick_Kept ==1 & Neighbors_Brick <=3 & Disorder_Sum >=6);
    
    %% Show results to figure
    disp('Creating figure')
    f = figure;
    if crop == 1
        imshow(Image_cropped);
    else
        imshow(Im_original);
    end
    
    hold on;
    %% Calculate lines to draw as colored labels on major axis of each particle
    %AMORPHOUS
    for particle = 1:number_of_particles
        if (Disorder_Kept(particle) == 1) && (Disorder_Kept(particle) ~= Brick_Kept(particle))  %only take disordered particles that are not brick particles (and thus not shared particles)
            %calculate line ot draw as red on major axis
            
            xMajor1 = stats(particle).Centroid(1) + (((stats(particle).MajorAxisLength)./2) * cosd(stats(particle).Orientation));
            xMajor2 = stats(particle).Centroid(1) - (((stats(particle).MajorAxisLength)./2) * cosd(stats(particle).Orientation));
            yMajor1 = stats(particle).Centroid(2) - (((stats(particle).MajorAxisLength)./2) * sind(stats(particle).Orientation));
            yMajor2 = stats(particle).Centroid(2) + (((stats(particle).MajorAxisLength)./2) * sind(stats(particle).Orientation));
            line([xMajor1 xMajor2],[yMajor1 yMajor2],'color','r', 'LineWidth',2)
        end
        
        %Brick
        if (Brick_Kept(particle) > 0) && (Shared_Kept(particle) == 0) %only brick particles
            %calculate line ot draw as blue on major axis
            
            xMajor1 = stats(particle).Centroid(1) + (((stats(particle).MajorAxisLength)./2) * cosd(stats(particle).Orientation));
            xMajor2 = stats(particle).Centroid(1) - (((stats(particle).MajorAxisLength)./2) * cosd(stats(particle).Orientation));
            yMajor1 = stats(particle).Centroid(2) - (((stats(particle).MajorAxisLength)./2) * sind(stats(particle).Orientation));
            yMajor2 = stats(particle).Centroid(2) + (((stats(particle).MajorAxisLength)./2) * sind(stats(particle).Orientation));
            line([xMajor1 xMajor2],[yMajor1 yMajor2],'color','blue', 'LineWidth',2);
            
        end
        
        %Shared
        if (Shared_Kept(particle) == 1) %Shared must be equal to 1
            %calculate line ot draw as green on major axis
            
            xMajor1 = stats(particle).Centroid(1) + (((stats(particle).MajorAxisLength)./2) * cosd(stats(particle).Orientation));
            xMajor2 = stats(particle).Centroid(1) - (((stats(particle).MajorAxisLength)./2) * cosd(stats(particle).Orientation));
            yMajor1 = stats(particle).Centroid(2) - (((stats(particle).MajorAxisLength)./2) * sind(stats(particle).Orientation));
            yMajor2 = stats(particle).Centroid(2) + (((stats(particle).MajorAxisLength)./2) * sind(stats(particle).Orientation));
            line([xMajor1 xMajor2],[yMajor1 yMajor2],'color','green', 'LineWidth', 2);
        end
        
        %NEITHER - amorphous particles
        if (Disorder_Kept(particle) == Brick_Kept(particle)) && (Brick_Kept(particle) == 0) %both perpendicular and parallel = 0
            plot(cns(particle,1),cns(particle,2),'blackX') %label as black X
        end
        
    end
    
    %legend
    hold on
    labels = {'Amorphous','Brick', 'Shared', 'Disorder'};
    category_type = categorical(labels); %Categories to be labeled in histogram
    category_type = reordercats(category_type,labels); %Preserves order of category_type.
    
    legend_dummy(1) = plot(nan, nan, 'red');
    legend_dummy(2) = plot(nan, nan, 'blue');
    legend_dummy(3) = plot(nan, nan, 'green');
    legend_dummy(4) = plot(nan, nan, 'blackX');
    
    legend(legend_dummy,labels); %create a legend.
    %Save labeled picture
    cd(file_path);
    cd Image_Files;
    saveas(f, [file_name_base, '_labeled.png']);
    
    %% Calculate quantitative data for each class
    
    
    for particle = 1:number_of_particles
        if (Disorder_Kept(particle) == 1) && (Brick_Kept(particle) == 0) % red particles
            Amorphous_Total(particle) = 1;
        else
            Amorphous_Total(particle) = 0;
        end
    end
    
    for particle = 1:number_of_particles
        if (Brick_Kept(particle) == 1) && (Shared_Kept(particle) == 0)
            Brick_Total(particle) = 1;
        else
            Brick_Total(particle) = 0; % Blue particles
        end
    end
    
    for particle = 1:number_of_particles
        if (Shared_Kept(particle) == 1)% Green particles
            Shared_Total(particle) = 1;
        else
            Shared_Total(particle) = 0;
        end
    end
    
    for particle = 1:number_of_particles
        if  (Brick_Kept(particle) == 0) && (Disorder_Kept(particle) == 0) %black particles
            Disorder_Total(particle) = 1;
        else
            Disorder_Total(particle) = 0;
        end
    end
    
    
    %% Output Quantitative Data
    
    %find percents
    Amorphous_Percent = ((sum(Amorphous_Total(:,1)/number_of_particles) * 100));
    Brick_Percent = ((sum(Brick_Total(:,1)/number_of_particles) * 100));
    Shared_Percent = ((sum(Shared_Total(:,1)/number_of_particles) * 100));
    Disorder_Percent = ((sum(Disorder_Total(:,1)/number_of_particles) * 100));
    
    %Write table
    colNames = {'Amorphous_Percent','Brick_Percent', 'Shared_Percent','Disorder_Percent', 'Amorphous_Total', ...
        'Brick_Total', 'Shared_Total', 'Disorder_Total',  'Particles_Total' };
    Output = table( Amorphous_Percent, Brick_Percent, Shared_Percent, Disorder_Percent,(sum(Amorphous_Total(:,1))),...
        (sum(Brick_Total(:,1))), (sum(Shared_Total(:,1))), (sum(Disorder_Total(:,1))), number_of_particles, 'VariableNames',colNames);
    
    
    cd(file_path);
    cd Excel_Files ;
    writetable(Output,[file_name_base, '_Percents', '.xlsx']); %writes analysis to excel file
    
    %write out parameters for file
    cd(file_path)
    cd Excel_Files;
    colNames2 = {'apara', 'areamin', 'bandpass_bpass', 'bandpass_strel', 'bwareafilt_param',...
        'bwareaopen_param', 'crop',  'distmax_dis', 'distmax_parallel', 'distmin_dis', ...
        'distmin_parallel', 'imgauss_param', 'malmin', 'perimin', ...
        'x1', 'x2', 'y1', 'y2'};
    
    Output2 = table(apara, areamin, bandpass_bpass, bandpass_strel, bwareafilt_param, ...
        bwareaopen_param, crop, distmax_dis, distmax_parallel, distmin_dis,...
        distmin_parallel,  imgauss_param, malmin, perimin, ...
        x1, x2, y1, y2, 'VariableNames',colNames2);
    writetable(Output2,[file_name_base, '_Parameters', '.xlsx']); %writes to excel
    
    disp('Done with file')
end

disp('Done with all files!')


