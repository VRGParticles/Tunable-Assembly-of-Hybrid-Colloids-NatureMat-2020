%% Analysis of Brick (cmm) vs. Square (p4m) particle assemblies
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

disp('Brick_Square.m Copyright (C) 2020 Veronica Grebe, Weck Group at New York University')
disp('This program comes with ABSOLUTELY NO WARRANTY.')
disp('This is free software released under GNU GPL-3.0,')
disp('and you are welcome to redistribute it under certain conditions;')
disp('see <http://www.gnu.org/licenses/>')

%   Use: Label an optical micrograph into cmm symmtery, pgg symmetry, boundaried
%       and disordered particle assemblies.

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
full_file_path = fullfile(file_path, '*.jpg'); % Change to whatever pattern your micrograph is.  Likely .jpg or .tif.  For one file, copy file name.  For multiple in the file_path, use a wildcard.

%Image Preprocessing.
imgauss_param = 3; %gaussian filter parameter;
bwareaopen_param = 30; %how large must the area be to keep the blob.
bwareafilt_param = [ 80 150 ]; %lower and upper area bound of kept blobs
bandpass_bpass = 1; %bandpass for bpass function
bandpass_strel = 2; %bandpass for oval/struture element

%distance cutoffs of centers of each particle.
% *** Note that in the script,
%constants are added to these numbers for the analysis of the images.
%These constants may need to be changed depending on the input images.***
distmin_parallel = 0; %the min distance between particle centers  for parallel particles, initially
distmax_parallel = 23; %the max distance between particle centers for parallel particles, initially
distmin_perp = 0; %the min distance between particle centers  for perp particles
distmax_perp = 27; %the max distance between particle centers for perp particles


%angle cutoffs.
apara = 15; %parallel bounds, + or - about 0
aperp1 = 70;   %perpendicular lower bounds set 1
aperp2 = 110;  %perpendicular upper bounds set 2

%blob property cutoffs
malmin = 7; %minimum length a particle must be (in major axis)
areamin = 70; %minimum area a particle must have
perimin = 15; %minimum perimeter a particle must have

%CROPPING
crop = 0; %set 1 for yes, and to any other number for no.
%for cropping, sets pixel distance.
x1 = 0; %top left
y1 = 200; %bottom left
x2 =  1400; %top right
y2 = 2000; %bottom right

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
    %bandpass filter applied to image using function call bpass.first_neighbor
    Image_processed = uint16(bpass(Image_processed,0,bandpass_bpass));
    
    %adjust image intensity to further segment out particles.
    Im_adjust = imadjust(Image_processed);
    
    greythresh_thresh = graythresh(Im_adjust); %find global image threshold
    
    imbw = Im_adjust > greythresh_thresh; %segmented, binary image
    
    %% Image segmentation: segment binary image
    
    se = strel('disk',bandpass_strel); %defines structure element to filter image based on shape.
    
    Image_opened_se = imopen(imbw>0,se); %morphological opening of image using structure element se.
    
    Image_opened_area = bwareaopen(Image_opened_se, bwareaopen_param); %remove small features
    
    Image_opened_area_filter = bwareafilt(Image_opened_area, bwareafilt_param); %remove small and large features with an accepted area range.
    
    %% Obtain blob properties
    
    [connected_binary,number_of_particles] = bwlabel(Image_opened_area_filter); %connect adjacent white pixels in segmented binary image
    
    stats = regionprops(connected_binary,'Area', 'Centroid', 'MajorAxisLength', 'Orientation', 'Perimeter'); %obtain defined properties for the connected regions of connected_binary
    
    
    %% Preallocation of arrays
    area = zeros(number_of_particles,1); %area from regionprops
    Brick_Total = zeros(number_of_particles,1); %Total number of Brick particles
    cns = zeros(number_of_particles,2); %centroid (x,y) from regionprops
    Dif_Ang = zeros(number_of_particles); %orientation angle difference
    Disorder_Total = zeros(number_of_particles,1); %Total number of Disorder particles
    Herringbone_Total = zeros(number_of_particles,1); %Total number of Herringbone particles
    mal = zeros(number_of_particles,1); %major axis length from regionprops
    Neighbors_Brick =  zeros(number_of_particles,1); %count of brick neighbors.
    Neighbors_Herringbone = zeros(number_of_particles,1); %count of knit neighbors.
    Neighbors_Parallel=  zeros(number_of_particles,1); %preallocate all the neighbors so the addition of tmps are reset to 0.
    Neighbors_Perp = zeros(number_of_particles,1); %count of perp neighbors.
    ori = zeros(number_of_particles,1); %orientation from regionprops
    Para_Check = zeros(number_of_particles,1); %Square particles that are parallel to each other (sides of a square).
    Brick_Kept = zeros(number_of_particles,1); %Particles marked as Brick
    Brick_Kept_0= zeros(number_of_particles,1); %particles that are kept meeting para. criteria, intitally
    Brick_Kept_not_neither =  zeros(number_of_particles,1); %particles that are classified as perp and are not neither.
    peri = zeros(number_of_particles,1); %perimeter from regionprops
    Square_Kept = zeros(number_of_particles,1);%Particles marked as Square
    Square_Kept_0= zeros(number_of_particles,1); %particles that are kept meeting perp. criteria, intitally
    Square_Kept_not_neither =  zeros(number_of_particles,1); %particles that are classified as perp and are not neither.
    Square_Kept_not_neither_2 =  zeros(number_of_particles,1); %particles that are classified as perp and are not neither, redefined.
    Shared_Total = zeros(number_of_particles,1); %Total number of Shared particles
    Square_Total = zeros(number_of_particles,1); %Total number of Square particles

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
    
    
    [first_neighbor, second_neighbor] = find(dist < (distmax_perp + 33));   %determine which particles are neighboring, with some tolerance
    
    %Sweep through all neighboring particles and only calculate their
    %orientation angle difference if they meet any of the threshold
    %conditions for perimeter, major axis length, or area.
    
    for particle = 1:length(first_neighbor)
        particle1=first_neighbor(particle); %read in the ith particle as a function
        particle2=second_neighbor(particle); %read in the jth particle as a function
        if mal(particle1) > malmin || area(particle1) > areamin || peri(particle1) > perimin %ensure the particle large enough in at least one of these properties
            
            Dif_Ang(particle1, particle2) = abs(stats(particle1).Orientation - stats(particle2).Orientation); %orientation angle difference
        end
        
    end
    
    %% Find Parallel edges of the perp shapes (solve for minimum of 3 neighbors)
    
    %First, find particles that are parallel to each other in a particular
    %distance threshold.
    [first_neighbor,second_neighbor] = find(dist > (distmax_perp -1) & dist < (distmax_perp + 3));
    for particle = 1:length(first_neighbor) %recalc neigbors
        particle1 = first_neighbor(particle);
        particle2 = second_neighbor(particle);
        
        %A particle pair is perpendicular to each other, and falls within
        %the extended range of being two parallel sides of a p4m configuration.
        
        if Dif_Ang(particle1, particle2) < apara && Dif_Ang(particle1, particle2) > 0 || Dif_Ang(particle1, particle2) < 180 && Dif_Ang(particle1, particle2) > (180 - apara) % only take parallel edges of the square.
            Para_Check(particle1) = 1;
            Para_Check(particle2) = 1;
        end
    end
    
    %% Plus Shape
    disp('Finding Shapes')
    
    AllowedDistance_Perp = (distmin_perp < dist & dist < distmax_perp); %neighboring particles that meet distance criteria
    
    AllowedAngles_Perp = (aperp1 < Dif_Ang & Dif_Ang < aperp2) | (-aperp1 > Dif_Ang & Dif_Ang > -aperp2); %neighboring particles that meet angle criteria for perp
    
    KeptParticlesParam_Perp = (AllowedDistance_Perp == AllowedAngles_Perp) & (AllowedDistance_Perp ==1);  %determine if particles can be accepted by perp if they meet both thresholds
    
    Perp_Sum=sum(KeptParticlesParam_Perp,2);% move to one column: how many neighbors each particle meets the above perp criteria.
    
    
    [first_neighbor, second_neighbor] = find(KeptParticlesParam_Perp == 1);  %find which particles meet the perp criteria.
    
    %Check for each pair, is it adajcent to a perp. particle that met the para
    %check? Meaning the formation of particle 1, particle2, and the previously
    %checked pair (particle1, and the other para particle) make up 3 sides of
    %a square structure (bundles of 3 as such: | _ | ).
    
    for particle=1:length(first_neighbor)
        
        particle1=first_neighbor(particle); 
        particle2=second_neighbor(particle); 
        
        %particle 1 has at least 2 perp. neighbors, and particle2 met the
        %Para_Check.
        if Perp_Sum(particle1) > 1 && Para_Check(particle2) == 1
            Square_Kept(particle1) = 1;
            Square_Kept(particle2) = 1;
        else
            Square_Kept(particle) = 0;
        end
    end
    
    %%  Parallel particles
    
    AllowedDistance_Parallel = (distmin_parallel < dist & dist < distmax_parallel);  %neighboring particles that meet distance criteria
    
    AllowedAngles_Parallel = (apara > Dif_Ang & Dif_Ang > -apara); %neighboring particles that meet angle criteria for parallel
    
    KeptParticlesParam_Parallel = (AllowedDistance_Parallel == AllowedAngles_Parallel) & (AllowedDistance_Parallel ==1); %determine if particles can be accepted by parallel if they meet both thresholds
    
    Parallel_Sum=sum(KeptParticlesParam_Parallel,2); % move to one column: how many neighbors each particle meets the above para criteria.
    
    [first_neighbor, second_neighbor] = find(KeptParticlesParam_Parallel == 1); %find which particles meet the para criteria.
    
    %Verify and count for each pair particle1, particle2, that they have
    %para neighbors (at least 1).
    for particle = 1:length(first_neighbor)
        particle1 = first_neighbor(particle);
        particle2 = second_neighbor(particle);
        if Parallel_Sum(particle1) > 0 && Parallel_Sum(particle2) > 0    %Neighbors must also be parallel.
            Brick_Kept_0(particle1) =  Brick_Kept_0(particle1) + 1;
        end
    end
    
    %If a particle has at least 2 parallel neighbors (that themselves have
    %para neighbors), they are marked as parallel particles.  This
    %establishes that parallel particles must be in bundles of 3.
    for particle= 1:length(Brick_Kept_0)
        if Brick_Kept_0(particle) > 1
            Brick_Kept(particle) = 1;  %NOTE: Brick_Kept is less stringent than Brick_Kept, will have more particles.
        end
    end
    
    
    
    %% Neighbor Requirement - parallel neigbhors of perp particles.
    clear first_neighbor second_neighbor
    [first_neighbor,second_neighbor] = find(dist < (distmax_parallel + 2));
    
    for particle = 1:length(first_neighbor)
        particle1 = first_neighbor(particle);
        particle2 = second_neighbor(particle);
        %If Particle 1 is perp, and particle2 is parallel, and the angle
        %difference between them is less than apara, particle 1 has a para
        %neighbor.
        if Square_Kept(particle1) == 1 && Brick_Kept(particle2) == 1 && Square_Kept(particle2) == 0 && (Dif_Ang(particle1, particle2) < apara)
            Neighbors_Parallel(particle1) = Neighbors_Parallel(particle1) + 1;
        end
    end
    
    %If a particle has at least 2 para neighbors, it is also para.
    for particle = 1:number_of_particles
        if Neighbors_Parallel(particle) > 1
            Brick_Kept(particle) = 1;
        end
    end
    
    %Repeat, accounting for the fact that there are now new para particles.
    for particle = 1:length(first_neighbor) %recalc neigbors
        particle1 = first_neighbor(particle);
        particle2 = second_neighbor(particle);
        if Square_Kept(particle1) == 1 && Brick_Kept(particle2) == 1 && Dif_Ang(particle1, particle2) < apara
            Neighbors_Parallel(particle1) = Neighbors_Parallel(particle1) + 1;
        end
    end
    
    %If there are new perp particles that have a new para neighbor, count
    %them as parallel.
    for particle= 1:number_of_particles
        if Neighbors_Parallel(particle) > 0
            Brick_Kept(particle) = 1;
        end
    end
    
    %%   Parallel Neighbors of Para particles to determine if non-structured particles should be classified as para.
    %reset neighbor counts.
    Neighbors_Parallel=  zeros(number_of_particles,1); %para neighbors of non-para particles.
    Neighbors_Parallel_para =  zeros(number_of_particles,1); %para neighbor of para particles.
    
    %If a para particle1 has a para neighbor, count it a para neighbor.
    for particle = 1:length(first_neighbor) %recalc neigbors
        particle1 = first_neighbor(particle);
        particle2 = second_neighbor(particle);
        if Brick_Kept(particle1) == 1 && Brick_Kept(particle2) == 1
            Neighbors_Parallel_para(particle1) = Neighbors_Parallel_para(particle1) + 1;
        end
    end
    
    %For non-perp and non-para particles, particle1, that have para neighbors, particle 2,
    %(who themselves have at
    %least one para neighbor), check if they are para.  If they are, that
    %currently non-structured particle1 has a para neighbor
    for particle = 1:length(first_neighbor) %recalc neigbors
        particle1 = first_neighbor(particle);
        particle2 = second_neighbor(particle);
        if Square_Kept(particle1) == 0 && Brick_Kept(particle1) == 0 && ...
                Brick_Kept(particle2) == 1 && Neighbors_Parallel_para(particle2) > 0 && Dif_Ang(particle1, particle2) < apara
            Neighbors_Parallel(particle1) =  Neighbors_Parallel(particle1) + 1;
        end
    end
    
    %If the particle meets the above criteria at least once and is
    %currently marked as non-para, mark it as para.
    for particle= 1:number_of_particles
        if Neighbors_Parallel(particle) > 0 && Brick_Kept(particle) == 0
            Brick_Kept(particle) = 1;
        end
    end
    
    %% Repeat the above process for the newly counted para particles.
    Neighbors_Parallel=  zeros(number_of_particles,1); %reset neighbors
    Neighbors_Parallel_para =  zeros(number_of_particles,1);
    
    %Now that select non-para and non-perp particles are marked as para, recalc para neighors of
    %para particles.
    for particle = 1:length(first_neighbor) %recalc neigbors
        particle1 = first_neighbor(particle);
        particle2 = second_neighbor(particle);
        if Brick_Kept(particle1) == 1 && Brick_Kept(particle2) == 1
            Neighbors_Parallel_para(particle1) = Neighbors_Parallel_para(particle1)+ 1;
        end
    end
    
    %Calculate for currently non-structured particles particle1, do they
    %have a para neighbor particle2 that has at least 1 para neighbor, AND
    %particle1 and 2 are parallel (repeat of above section, but with newly
    %counted particles).
    
    for particle = 1:length(first_neighbor)
        particle1 = first_neighbor(particle);
        particle2 = second_neighbor(particle);
        if Square_Kept(particle1) == 0 && Brick_Kept(particle1) == 0 && ...
                Brick_Kept(particle2) == 1 && Neighbors_Parallel_para(particle2) > 0 && Dif_Ang(particle1, particle2) < apara
            Neighbors_Parallel(particle1) =  Neighbors_Parallel(particle1) + 1;
        end
    end
    
    % If a parallel particle has only 0-1 para neighbor, remove it (bundles of 3+ for para neighbors).
    for particle= 1:number_of_particles
        if Neighbors_Parallel_para(particle) < 2 && Brick_Kept(particle) == 1
            Brick_Kept(particle) = 0;
        end
    end
    
    
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
    %Square
    for particle = 1:number_of_particles
        if (Square_Kept(particle) == 1) && (Brick_Kept(particle) == 0)  %only take perp particles
            %calculate line ot draw as red on major axis
            
            xMajor1 = stats(particle).Centroid(1) + (((stats(particle).MajorAxisLength)./2) * cosd(stats(particle).Orientation));
            xMajor2 = stats(particle).Centroid(1) - (((stats(particle).MajorAxisLength)./2) * cosd(stats(particle).Orientation));
            yMajor1 = stats(particle).Centroid(2) - (((stats(particle).MajorAxisLength)./2) * sind(stats(particle).Orientation));
            yMajor2 = stats(particle).Centroid(2) + (((stats(particle).MajorAxisLength)./2) * sind(stats(particle).Orientation));
            line([xMajor1 xMajor2],[yMajor1 yMajor2],'color','r', 'LineWidth',2)
        end
        
        %Brick
        if (Brick_Kept(particle) == 1) && (Square_Kept(particle) == 0) %only take para
            %calculate line ot draw as blue on major axis
            
            xMajor1 = stats(particle).Centroid(1) + (((stats(particle).MajorAxisLength)./2) * cosd(stats(particle).Orientation));
            xMajor2 = stats(particle).Centroid(1) - (((stats(particle).MajorAxisLength)./2) * cosd(stats(particle).Orientation));
            yMajor1 = stats(particle).Centroid(2) - (((stats(particle).MajorAxisLength)./2) * sind(stats(particle).Orientation));
            yMajor2 = stats(particle).Centroid(2) + (((stats(particle).MajorAxisLength)./2) * sind(stats(particle).Orientation));
            line([xMajor1 xMajor2],[yMajor1 yMajor2],'color','blue', 'LineWidth',2);
            
        end
        
        %BOTH
        if (Square_Kept(particle) == 1) && ( Brick_Kept(particle) == 1) % The particle is both para and perp
            %calculate line ot draw as green on major axis
            xMajor1 = stats(particle).Centroid(1) + (((stats(particle).MajorAxisLength)./2) * cosd(stats(particle).Orientation));
            xMajor2 = stats(particle).Centroid(1) - (((stats(particle).MajorAxisLength)./2) * cosd(stats(particle).Orientation));
            yMajor1 = stats(particle).Centroid(2) - (((stats(particle).MajorAxisLength)./2) * sind(stats(particle).Orientation));
            yMajor2 = stats(particle).Centroid(2) + (((stats(particle).MajorAxisLength)./2) * sind(stats(particle).Orientation));
            line([xMajor1 xMajor2],[yMajor1 yMajor2],'color','green', 'LineWidth', 2);
        end
        
        % Disorder
        if (Square_Kept(particle) == 0) && (Brick_Kept(particle) == 0) %both perp and parallel not met
            plot(cns(particle,1),cns(particle,2),'blackX') %label as black X
        end
        
    end
    
    %Save labeled picture
    cd(file_path);
    cd Image_Files;
    saveas(f, [file_name_base, '_labeled.png']);
    
    %% Calculate quantitative data for each class.
    for particle = 1:number_of_particles
        if Square_Kept(particle) == 1 && Brick_Kept(particle) == 1
            Shared_Total(particle) = 1;
        else
            Shared_Total(particle) = 0;
        end
    end
    
    for particle = 1:number_of_particles
        if Square_Kept(particle) == 1 && Shared_Total(particle) == 0
            Square_Total(particle) = 1;
        else
            Square_Total(particle) = 0;
        end
    end
    
    for particle = 1:number_of_particles
        if Brick_Kept(particle) == 1 && Shared_Total(particle) == 0
            Brick_Total(particle) = 1;
        else
            Brick_Total(particle) = 0;
        end
    end
    
    for particle = 1:number_of_particles
        if Square_Kept(particle) == 0 && Brick_Kept(particle) == 0
            Disorder_Total(particle) = 1;
        else
            Disorder_Total(particle) = 0;
        end
    end
    %% Output Quantitative Data
    %find percents
    Square_Percent = ((sum(Square_Total(:,1)/number_of_particles) * 100));
    Brick_Percent = ((sum(Brick_Total(:,1)/number_of_particles) * 100));
    Shared_Percent = ((sum(Shared_Total(:,1)/number_of_particles) * 100));
    Disorder_Percent = ((sum(Disorder_Total(:,1)/number_of_particles) * 100));
    
    %Write table
    colNames = {'Square_Percent','Brick_Percent', 'Shared_Percent', 'Disorder_Percent',  'Square_Total', ...
        'Brick_Total','Shared_Total', 'Disorder_Total',  'Particles_Total' };
    Output = table( Square_Percent, Brick_Percent, Shared_Percent, Disorder_Percent,(sum(Square_Total(:,1))),...
        (sum(Brick_Total(:,1))), (sum(Shared_Total(:,1))), (sum(Disorder_Total(:,1))), number_of_particles, 'VariableNames',colNames);

    
    
    cd(file_path);
    cd Excel_Files ;
    writetable(Output,[file_name_base, '_Percents', '.xlsx']); %writes analysis to excel per file.
    
    %write out parameters once
    
    colNames2 = {'apara', 'aperp1', 'aperp2', 'areamin', 'bandpass_bpass', 'bandpass_strel',...
        'bwareafilt_param', 'bwareaopen_param', 'crop',  ...
        'distmax_parallel', 'distmax_perp', 'distmin_parallel', 'distmin_perp',' imgauss_param', ...
        'malmin', 'perimin', 'x1', 'x2', 'y1', 'y2'};
    Output2 = table(apara, aperp1, aperp2, areamin, bandpass_bpass, bandpass_strel, ...
        bwareafilt_param, bwareaopen_param, crop, ...
        distmax_parallel, distmax_perp, distmin_parallel, distmin_perp, imgauss_param,...
        malmin, perimin, x1, x2, y1, y2, 'VariableNames',colNames2);
    writetable(Output2,[file_name_base, '_Parameters', '.xlsx']); %writes importable file to excel
    
    cd(file_path);
    disp('Done with file')
    
    
end
disp('Done with all files!')


