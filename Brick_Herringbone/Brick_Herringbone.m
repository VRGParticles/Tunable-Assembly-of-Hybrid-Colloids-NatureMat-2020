%% Analysis of Brick (cmm) vs. Herringbone (pgg) particle assemblies
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

disp('Brick_Herringbone.m Copyright (C) 2020 Veronica Grebe, Weck Group at New York University')
disp('This program comes with ABSOLUTELY NO WARRANTY.')
disp('This is free software released under GNU GPL-3.0,')
disp('and you are welcome to redistribute it under certain conditions;')
disp('see <http://www.gnu.org/licenses/>')

%   Use: Label an optical micrograph into cmm symmtery, pgg symmetry, boundaried
%       and disordered particle assemblies.

%   Author Contact: Veronica Grebe vrg234@nyu.edu, Weck Group, New York University.

%   Function calls: This script calls the function 'bpass.pro'during image segmentation and was released under the following license:
%       This function, bpass.pro, is copyright:
%       John C. Crocker and David G. Grier.  It should be considered 'freeware'- and may be
%       distributed freely in its original form when properly attributed.
%       https://es.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/42573/versions/1/previews/bpass.first_neighbor/index.html

%%  Input parameters **EDIT HERE**
clearvars %clear any current variables in the workspace.

%file location and type
file_path = '/where/your/files/live/'; %where your file(s) to be analyzed are located.
full_file_path = fullfile(file_path, '*.jpg'); % Change to whatever pattern your micrograph is.  Likely .jpg or .tif.  For one file, copy file name.  For multiple in the file_path, use a wildcard.


%Image Preprocessing.
imgauss_param = 2; %gaussian filter parameter;
bwareaopen_param = 30; %how large must the area be to keep the blob.
bwareafilt_param = [ 80 160 ]; %lower and upper area bound of kept blobs
bandpass_bpass = 1; %bandpass for bpass function
bandpass_strel = 1; %bandpass for oval/struture element

%distance cutoffs of centers of each particle.
% *** Note that in the script,
%constants are added to these numbers for the analysis of the images.
%These constants may need to be changed depending on the input images.***
distance_neighbors = 27; %Distance neighbors are apart, used for calculations for classification.
distmin_parallel = 0; %the min distance between particle centers  for parallel particles, initially
distmax_parallel = 25; %the max distance between particle centers for parallel particles, initially
distmin_perp = 0; %the min distance between particle centers  for perp particles
distmax_perp = 25; %the max distance between particle centers for perp particles


%angle cutoffs.
apara = 10; %parallel bounds, + or - about 0
aperp1 = 60;   %perpendicular lower bounds set 1
aperp2 = 120;  %perpendicular upper bounds set 2

%blob property cutoffs
malmin = 7; %minimum length a particle must be (in major axis)
areamin = 70; %minimum area a particle must have
perimin = 15; %minimum perimeter a particle must have
discent1min = 0; %minimum the center of one particle can be from the end of the other particle

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
    Brick_Total = zeros(number_of_particles,1); % Particles that are classfied as Brick
    cns = zeros(number_of_particles,2); %centroid (x,y) from regionprops
    Dif_Ang = zeros(number_of_particles); %orientation angle difference
    Disorder_Total= zeros(number_of_particles,1); %particles that are classified as disordered
    Herringbone_Total = zeros(number_of_particles,1); %particles that are classified as Herringbone
    mal = zeros(number_of_particles,1); %major axis length from regionprops
    Neighbors_Brick =  zeros(number_of_particles,1); %count of brick neighbors.
    Neighbors_Herringbone = zeros(number_of_particles,1); %count of herringbone neighbors.
    Neighbors_Perp = zeros(number_of_particles,1); %count of perp neighbors.
    ori = zeros(number_of_particles,1); %orientation from regionprops
    Parallel_Kept_0 =  zeros(number_of_particles,1); %particles that are classified as perp and are not disorder.
    peri = zeros(number_of_particles,1); %perimeter from regionprops
    Perp_Kept_0= zeros(number_of_particles,1); %particles that are kept meeting perp. criteria, intitally
    Perp_Kept_00 =  zeros(number_of_particles,1); %particles that are classified as perp and are not disorder.
    Perp_Kept_000 =  zeros(number_of_particles,1); %particles that are classified as perp and are not disorder, redefined.
    Shared_Total = zeros(number_of_particles,1); %particles that are classified as shared

    %% Index Region Stats
    %obtain centroids, orientation, major axis length, area and perimeter of
    %connected regions determined by 'stats'above
    
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
    
    
    [first_neighbor, second_neighbor] = find(dist < distance_neighbors + 13);   %determine which particles are neighboring, with some tolerance
    
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
    
    %% Perp Shape
    
    AllowedDistance_Perp = (distmin_perp < dist & dist < distmax_perp); %neighboring particles that meet distance criteria
    
    AllowedAngles_Perp = (aperp1 < Dif_Ang & Dif_Ang < aperp2) | (-aperp1 > Dif_Ang & Dif_Ang > -aperp2); %neighboring particles that meet angle criteria for perp
    
    KeptParticlesParam_Perp = (AllowedDistance_Perp == AllowedAngles_Perp) & (AllowedDistance_Perp ==1);  %determine if particles can be accepted by perp if they meet both thresholds
    
    Perp_Sum=sum(KeptParticlesParam_Perp,2);% move to one column: how many neighbors each particle meets the above perp criteria for.
    
    [first_neighbor, second_neighbor] = find(KeptParticlesParam_Perp == 1);  %find which particles meet the criteria.
    
    
    %Verify and count for each pair particle1, particle2, that they have
    %perp neighbors. If both particles have at least 1 Perp_Sum, then
    %count towards particle1 Perp_Kept_0.
    
    for particle = 1:length(first_neighbor)
        particle1 = first_neighbor(particle);
        particle2 = second_neighbor(particle);
        
        if Perp_Sum(particle1) > 0 && Perp_Sum(particle2) > 0
            Perp_Kept_0(particle1) = Perp_Kept_0(particle1) + 1; %count of how many particle pairs meet the perp criteria
        end
    end
    
    clear first_neighbor second_neighbor
    %%  Parallel particles
    
    AllowedDistance_Parallel = (distmin_parallel < dist & dist < distmax_parallel);  %neighboring particles that meet distance criteria
    
    AllowedAngles_Parallel = (apara > Dif_Ang & Dif_Ang > -apara) | ((180 - apara) < Dif_Ang & Dif_Ang < 180) ; %neighboring particles that meet angle criteria for parallel
    
    KeptParticlesParam_Parallel = (AllowedDistance_Parallel == AllowedAngles_Parallel) & (AllowedDistance_Parallel ==1); %determine if particles can be accepted by parallel if they meet both thresholds
    
    Parallel_Sum=sum(KeptParticlesParam_Parallel,2); % move to one column: how many neighbors each particle meets the above para criteria for.
    
    [first_neighbor, second_neighbor] = find(KeptParticlesParam_Parallel == 1); %find which particles meet the criteria.
    
    %Verify and count for each pair particle1, particle2, that they have
    %para neighbors.
    
    for particles = 1:length(first_neighbor)
        particle1 = first_neighbor(particles);
        particle2 = second_neighbor(particles);
        
        if Parallel_Sum(particle1) > 0 && Parallel_Sum(particle2) > 0    %Neighbors must also be parallel.
            Parallel_Kept_0(particle1) = Parallel_Kept_0(particle1)+ 1;
        end
    end
    
    Brick_Kept = (Parallel_Kept_0 >1 & Perp_Kept_0 == 0); %Particles are marked as brick if they have at least 2 parallel neighbors, and no perp neighbors.
    
    Parallel_Kept = Parallel_Kept_0 > 1; % Parallel_Kept is less stringent than Brick_Kept, will have more particles. Parallel particles, but with perp neighbos.
    %% Initial Classification.
    
    Herringbone_Kept = Parallel_Kept_0 > 0 & Perp_Kept_0 > 0; %loose, first classifiation of pgg particles.  They have at least 1 parallel and 1 perp neighbor.
    Disorder_Kept = Herringbone_Kept == 0 & Brick_Kept == 0; %Disorder particles are not members of structured classes.

    %% Remove 'Disorder' Neighbors for not meeting criteria for para.
    
    clear first_neighbor second_neighbor
    
    [first_neighbor, second_neighbor] = find(KeptParticlesParam_Perp == 1);
    
    %Particle pairs where particle1 and particle2 both have at least 1 Perp
    %neighbor, AND disorder is classified as disorder.  Removes disorder
    %particles being counted in the perp count.
    
    for particle = 1:length(first_neighbor)
        particle1 = first_neighbor(particle);
        particle2 = second_neighbor(particle);
        
        if Perp_Sum(particle1) > 0 && Perp_Sum(particle2) > 0 && Disorder_Kept(particle2) == 0 ...
                && Disorder_Kept(particle1) == 0  %count of how many particle pairs meet the perp criteria AND are not classified as disorder.
            Perp_Kept_00(particle1) = Perp_Kept_00(particle1) + 1;
        end
    end
    
    %Rewrite which particles are kept as perp, only counting those not
    %counted as disorder
    Perp_Kept = Perp_Kept_00 > 0;

    %Update classifications, now that disorder particles have been removed from
    %the perp count.
    Brick_Kept = Parallel_Kept_0 > 0 & Perp_Kept == 0;     %Update Brick_Kept for parallel particles that are not classfied as disorder, or perp.

    % Brick_Kept has at at least 1 particle pair that is para, and each member of the pair has at least one para neighbor. 
    Herringbone_Kept = Parallel_Kept_0 > 0 & Perp_Kept_00 >0; %Count herringbone particles as those that have at least 1 para and perp neighbor (above criteria)
    Shared_Kept = Herringbone_Kept == 1 & Brick_Kept == 1; %Count 'both'particles as both herrinbone and brick.

    %% Neighbor counts for each structure and new neighbor count conditions
    %herringbone must have at least 1 herringbone neighbor, and remove
    %chains of parallel brick particles.
    
    [second_neighbor,first_neighbor] = find(dist < distance_neighbors + 3); %Find particles that are 3 more than the distendsmin to check particle pairs
    
    for particle = 1:length(first_neighbor) %check each neighbor's properties.
        particle1 = first_neighbor(particle);
        particle2 = second_neighbor(particle);
        
        %For each class, count how many neighbors (with a tolerance of
        %3) are members of each class.
        
        if Herringbone_Kept(particle2) == 1
            Neighbors_Herringbone(particle1) = Neighbors_Herringbone(particle1) + 1;
        end
        
        if Brick_Kept(particle2) == 1
            Neighbors_Brick(particle1) = Neighbors_Brick(particle1) + 1;
        end
        
        if Perp_Kept(particle2) == 1
            Neighbors_Perp(particle1) = Neighbors_Perp(particle1) + 1;
        end
    end
    
    %Begin to assert neighbor rules.  The 'core'particles are beginning to
    %be located.  In this case, we remove particles below a neighbor
    %threshold.

    for particle = 1:number_of_particles
        if Herringbone_Kept(particle) == 1 && Neighbors_Herringbone(particle) < 1     %if herringbone has 0 herringbone neighbors, remove it.
            Herringbone_Kept(particle) = 0;
        end
        
        if Brick_Kept(particle) == 1 && Neighbors_Brick(particle) < 4 %if brick has 3 or fewer neighbors, remove it (removes lines of parallel particles).
            Brick_Kept(particle) = 0;
        end
        
    end
    
    Disorder_Kept = (Brick_Kept ==0 & Herringbone_Kept == 0 & Shared_Kept == 0); %redefine disorder.
    
    %% Recalculate perp based on disorder particles.
    [second_neighbor,first_neighbor] = find(dist < distance_neighbors + 3);
    
    %number_of_particlesw that disorder is redefined, we must redefine structured particles and remove disorder para particles.
    for particle = 1:length(first_neighbor)
        particle1 = first_neighbor(particle);
        particle2 = second_neighbor(particle);
        if Perp_Sum(particle1) > 0 && Perp_Sum(particle2) > 0 ...
                && dist(particle1,particle2) < (distance_neighbors -8) && Disorder_Kept(particle2) == 0% Only count perp structures that are NOT disorder.
            Perp_Kept_000(particle1) = Perp_Kept_000(particle1) + 1;
        end
    end
    
    Perp_Kept_2 = Perp_Kept_000 > 0; %Redefine Perp particles with newly updated disorder particles.
    
    %% Redefine herringbone neighbors: must have at least 2 herringbone neighbors.
    clear first_neighbor second_neighbor
    
    Neighbors_Herringbone=  zeros(number_of_particles,1); %reset neighbors
    
    [second_neighbor,first_neighbor] = find(dist < (distance_neighbors -2)); %find pairs.
    
    for particle = 1:length(first_neighbor) %check each neighbor's properties.
        particle1 = first_neighbor(particle);
        particle2 = second_neighbor(particle);
        
        if Herringbone_Kept(particle2) == 1
            Neighbors_Herringbone(particle1) = Neighbors_Herringbone(particle1) + 1; %count neighbor as herringbone
        end
    end
    
    %Assert more neighbor rules: hrrringbone must have at least 2
    %herringbone neighbors (bundle of 3)
    for particle = 1:number_of_particles
        if Herringbone_Kept(particle) == 1 && Neighbors_Herringbone(particle) < 2 %if herringbone has less than 2 neighbors, remove it.
            Herringbone_Kept(particle) = 0;
            Disorder_Kept(particle) = 1;
        end
    end
    
 
    %% Recalculate neighbors
    
    clear first_neighbor second_neighbor
    [first_neighbor,second_neighbor] = find(dist < distance_neighbors + 3);
    
    %reset all the neighbors so the addition of tmps are reset to 0.
    Neighbors_Herringbone=  zeros(number_of_particles,1);
    Neighbors_Brick=  zeros(number_of_particles,1);
    Neighbors_Parallel_Herringbone=  zeros(number_of_particles,1);
    
    
    for particle = 1:length(first_neighbor) %recalc neigbors
        particle1 = first_neighbor(particle);
        particle2 = second_neighbor(particle);
        if Herringbone_Kept(particle2) == 1
            Neighbors_Herringbone(particle1) = Neighbors_Herringbone(particle1)+ 1;
        end
        
        if Brick_Kept(particle2) == 1
            Neighbors_Brick(particle1) = Neighbors_Brick(particle1) + 1;
        end

        %For parallel neighbors of herringbone particles. 
        if Herringbone_Kept(particle2) == 1 && Herringbone_Kept(particle1) == 1 && Dif_Ang(particle1, particle2) < apara
            Neighbors_Parallel_Herringbone(particle1) = Neighbors_Parallel_Herringbone(particle1) + 1;
        end
    end

    
    %% Remove more particles based on neighbor count thresholds: Herringbone and Brick.  Must have at least 3 same class neighbors.

    %Assert that Herringbone and Brick particles must have at least 3
    %Herringbone or Brick neighbors, respectively.  This helps locates 'core'
    %structured particles.  This only counts particles that are in bundles
    %of 4 structured particles.
    
    for particle = 1:number_of_particles
        if Herringbone_Kept(particle) == 1 && Neighbors_Herringbone(particle) < 3 %if herringbone less than 3 neighbors, remove.
            Herringbone_Kept(particle) = 0;
        end
        
        if Brick_Kept(particle) == 1 && Neighbors_Brick(particle) < 3 %if brick less than 3 neighbors, remove.
            Brick_Kept(particle) = 0;
        end
    end
    %% Locate Shared Particles based on the above 'core'particles.
    Shared_Herringbone = Herringbone_Kept == 1 & Neighbors_Parallel_Herringbone > 3; %Assign herringbone particles with at least 4 para neighbors into likely boundaried particles

    %% Correct disorder particles that are para to core particles by reclassifying them, now that some core particles are classified as shared.
            % *** NOTE: apara is given more tolerance here.
            
    %Disorder to Brick
    Neighbors_Parallel = zeros(number_of_particles, 1); %reset neighbors
    
    for particle = 1:length(first_neighbor)
        particle1 = first_neighbor(particle);
        particle2 = second_neighbor(particle);

        % Disordered particles that are para to herringbone particles count
        % as parallel.
        if  Disorder_Kept(particle1) == 1 && Herringbone_Kept(particle2) == 1 && Dif_Ang(particle1, particle2) < apara + 5
            Neighbors_Parallel(particle1) = Neighbors_Parallel(particle1) + 1;
        end
    end
    
    
    % If a disordered particle has at least 2 para herringbone neighbors,
    % it is instead a Brick particle.
    for particle= 1:number_of_particles
        if Disorder_Kept(particle) == 1 && Neighbors_Parallel(particle) > 1 %If the disorder particle has at least 2 parallel, herringbone particles, it is a brick particle.
            Shared_Kept(particle) = 0;
            Herringbone_Kept(particle) = 0;
            Brick_Kept(particle) = 1;
            Disorder_Kept(particle) = 0;
        end
    end
    
    Neighbors_Parallel = zeros(number_of_particles, 1); %reset neighbors
    %% Find brick particles that are parallel to herringbone particles.
    % These particles run along the shared particles.
            % *** NOTE: apara is given more tolerance here.

    %Brick particles that are para to herringbone particles are counted.
    for particle = 1:length(first_neighbor)
        particle1 = first_neighbor(particle);
        particle2 = second_neighbor(particle);
        if  Brick_Kept(particle1) == 1 && Herringbone_Kept(particle2) == 1 && Dif_Ang(particle1, particle2) < apara + 5 %give a point if brick neighbors are parallel to herringbones
            Neighbors_Parallel(particle1) = Neighbors_Parallel(particle1) + 1;
        end
    end
        
    %If brick particles are para to herringbone particles, they are on the
    %Boundary of the two symmetry groups.
    
    %assign into Brick_Boundary 
    Brick_Boundary = (Brick_Kept == 1 & Neighbors_Parallel > 0);
    
%% Find more Herringbone particles that are actually shared particles.
%These particles are parallel to the Brick_Boundary particles, and are
%themselves the shared particles.

                % *** NOTE: apara is given more tolerance here.
    Neighbors_Parallel = zeros(number_of_particles, 1); %reset parallel neighbors

    for particle = 1:length(first_neighbor)
        particle1 = first_neighbor(particle);
        particle2 = second_neighbor(particle);
        if  Herringbone_Kept(particle1) == 1 && Brick_Boundary(particle2) == 1 && Dif_Ang(particle1, particle2) < apara + 5 %if herringbone is parallel to the shared particles give a point
            Neighbors_Parallel(particle1) =  Neighbors_Parallel(particle1) + 1;
        end
    end
    
 
    for particle= 1:number_of_particles
        if Herringbone_Kept(particle) == 1 && Neighbors_Parallel(particle) > 0 %assign Herringbone particles that are parallel to shared particles as shared particles.
            Shared_Kept(particle) = 1;
            Herringbone_Kept(particle) = 0;
            Brick_Kept(particle) = 0;
            Disorder_Kept(particle) = 0;
        end
    end
    
%% Overwrite Shared_Herringbone initial shared particles classification into Shared particles
% Move the herringbone particles with at least 4 para neighbors into
% Shared.

for particle = 1:number_of_particles
    if Shared_Herringbone(particle) ==1 %assign in tempherringbones
        Shared_Kept(particle) = 1;
        Herringbone_Kept(particle) = 0;
        Brick_Kept(particle) = 0;
        Disorder_Kept(particle) = 0;
    end
end

    %% Remove particles that are not in bundles of 4, INCLUDING SHARED PARTICLES.  Assert same  core neighbor rules.
    
    %Now that particles have moved to shared, assert the final neighbor
    %rules: particles must have at least 3 similarly structured particles,
    %or shared particles, neighbors.  This asserts that structured
    %particles must be in bundles of 4, but now includes shared particles.
    
    Neighbors_Brick = zeros(number_of_particles, 1); %Reset neighbors
    Neighbors_Herringbone = zeros(number_of_particles, 1);
    
    
    for particle = 1:length(first_neighbor)
        particle1 = first_neighbor(particle);
        particle2 = second_neighbor(particle);
        
        %Particles that are either herringbone or shared count as herringbone
        %neighbors.
        if Herringbone_Kept(particle2) == 1 || Shared_Kept(particle2) == 1
            Neighbors_Herringbone(particle1) = Neighbors_Herringbone(particle1) + 1;
        end
                %Particles that are either brick or shared count as brick
        %neighbors.
        if Brick_Kept(particle2) == 1 || Shared_Kept(particle2) == 1
            Neighbors_Brick(particle1) = Neighbors_Brick(particle1) + 1;
        end
    end
    
    %Only keep particles defined as herringbone and brick if they have at least 3
    %neighbors of the same class OR Boundary.  The final definition is that
    %structured particles kept must be in a bundle of at least 4 structured
    %particles.
    
    for particle = 1:number_of_particles
        if Herringbone_Kept(particle) == 1 && Neighbors_Herringbone(particle) < 3 %only keep if neighbors are 3+
            Herringbone_Kept(particle) = 0;
        end
        
        if Brick_Kept(particle) == 1 && Neighbors_Brick(particle) < 3
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
    %Herringbone, RED
    for particle = 1:number_of_particles
        if (Herringbone_Kept(particle) == 1) &&  (Brick_Kept(particle) == 0) && Shared_Kept(particle) == 0 %only take Perp particles
            %calculate line ot draw as red on major axis
            
            xMajor1 = stats(particle).Centroid(1) + (((stats(particle).MajorAxisLength)./2) * cosd(stats(particle).Orientation));
            xMajor2 = stats(particle).Centroid(1) - (((stats(particle).MajorAxisLength)./2) * cosd(stats(particle).Orientation));
            yMajor1 = stats(particle).Centroid(2) - (((stats(particle).MajorAxisLength)./2) * sind(stats(particle).Orientation));
            yMajor2 = stats(particle).Centroid(2) + (((stats(particle).MajorAxisLength)./2) * sind(stats(particle).Orientation));
            line([xMajor1 xMajor2],[yMajor1 yMajor2],'color','r', 'LineWidth',2)            
        end
        
        %Brick, BLUE
        if (Brick_Kept(particle) == 1) && (Herringbone_Kept(particle) == 0) && Shared_Kept(particle) == 0 %onyl take parallel
            %calculate line ot draw as blue on major axis

            xMajor1 = stats(particle).Centroid(1) + (((stats(particle).MajorAxisLength)./2) * cosd(stats(particle).Orientation));
            xMajor2 = stats(particle).Centroid(1) - (((stats(particle).MajorAxisLength)./2) * cosd(stats(particle).Orientation));
            yMajor1 = stats(particle).Centroid(2) - (((stats(particle).MajorAxisLength)./2) * sind(stats(particle).Orientation));
            yMajor2 = stats(particle).Centroid(2) + (((stats(particle).MajorAxisLength)./2) * sind(stats(particle).Orientation));
            line([xMajor1 xMajor2],[yMajor1 yMajor2],'color','blue', 'LineWidth',2);            
        end
        
        %SHARED, GREEN
        if (Shared_Kept(particle) ==1)  %Both must be equal to 1
            %calculate line ot draw as red on major axis

            xMajor1 = stats(particle).Centroid(1) + (((stats(particle).MajorAxisLength)./2) * cosd(stats(particle).Orientation));
            xMajor2 = stats(particle).Centroid(1) - (((stats(particle).MajorAxisLength)./2) * cosd(stats(particle).Orientation));
            yMajor1 = stats(particle).Centroid(2) - (((stats(particle).MajorAxisLength)./2) * sind(stats(particle).Orientation));
            yMajor2 = stats(particle).Centroid(2) + (((stats(particle).MajorAxisLength)./2) * sind(stats(particle).Orientation));
            line([xMajor1 xMajor2],[yMajor1 yMajor2],'color','green', 'LineWidth', 2);
        end
        
        %NEITHER, BLACK
        if  (Brick_Kept(particle) == 0) && (Herringbone_Kept(particle) == 0) && (Shared_Kept(particle) ==0 )  %both perp and parallel = 0
            plot(cns(particle,1),cns(particle,2),'blackX') %label as black X
        end
        
    end
    
    %Save labeled picture
        hold on
    labels = {'Herringbone','Brick','Shared', 'Disorder'};
    category_type = categorical(labels); %Categories to be labeled in histogram
    category_type = reordercats(category_type,labels); %Preserves order of category_type.
    
    legend_dummy(1) = plot(nan, nan, 'red');
    legend_dummy(2) = plot(nan, nan, 'blue');
    legend_dummy(3) = plot(nan, nan, 'green');
    legend_dummy(4) = plot(nan, nan, 'blackX');
    
    legend(legend_dummy,labels); %create a legend.
    cd(file_path);
    cd Image_Files;
    saveas(f, [file_name_base, '_labeled.png']);
    %% Calculate quantitative data for  all classes.
for particle = 1:number_of_particles
    if (Herringbone_Kept(particle) == 1) && (Shared_Kept(particle) == 0 )
        Herringbone_Total(particle) = 1;
    else
        Herringbone_Total(particle) = 0;
    end 
end 

for particle = 1:number_of_particles
    if (Brick_Kept(particle) == 1) && (Shared_Kept(particle) == 0)
        Brick_Total(particle) = 1;
    else
        Brick_Total(particle) = 0;
    end 
end

for particle = 1:number_of_particles
    if (Shared_Kept(particle) == 1)
        Shared_Total(particle) = 1;
    else
        Shared_Total(particle) = 0;
    end 
end 


for particle = 1:number_of_particles
    if (Brick_Kept(particle) == 0) && (Herringbone_Kept(particle) == 0) && (Shared_Kept(particle) ==0 )
        Disorder_Total(particle) = 1;
    else
        Disorder_Total(particle) = 0;
    end 
end 
    %% Output Quantitative Data
    %find percents
    Herringbone_Percent = ((sum(Herringbone_Total(:,1)/number_of_particles) * 100));
    Brick_Percent = ((sum(Brick_Total(:,1)/number_of_particles) * 100));
    Shared_Percent = ((sum(Shared_Total(:,1)/number_of_particles) * 100));
    Disorder_Percent = ((sum(Disorder_Total(:,1)/number_of_particles) * 100));
    
    %Write table
    colNames = {'Herringbone_Percent','Brick_Percent', 'Shared_Percent', 'Disorder_Percent', 'Herringbone_Total',  ...
        'Brick_Total', 'Shared_Total', 'Disorder_Total',  'Particles_Total'};
    Output = table(Herringbone_Percent,  Brick_Percent, Shared_Percent, Disorder_Percent,(sum(Herringbone_Total(:,1))), ...
        (sum(Brick_Total(:,1))),(sum(Shared_Total(:,1))) ,(sum(Disorder_Total(:,1))), number_of_particles, 'VariableNames',colNames);
    
    cd(file_path);
    cd Excel_Files ;
    writetable(Output,[file_name_base, '_Percents', '.xlsx']); %writes importable file to excel per file.  Can be used in
    %Excel or histogram script to analyze (Feb 2019).

    %write out parameters once
colNames2 = {'apara', 'aperp1', 'aperp2', 'areamin', 'bandpass_bpass', 'bandpass_strel',...
    'bwareafilt_param', 'bwareaopen_param', 'crop', 'distance_neighbors', ...
    'distmax_parallel', 'distmax_perp', 'distmin_parallel', 'distmin_perp','imgauss_param', ...
    'malmin', 'perimin', 'x1', 'x2', 'y1', 'y2'};
Output2 = table(apara, aperp1, aperp2, areamin, bandpass_bpass, bandpass_strel, ...
    bwareafilt_param, bwareaopen_param, crop,  distance_neighbors, ...
    distmax_parallel, distmax_perp, distmin_parallel, distmin_perp, imgauss_param,...
    malmin, perimin, x1, x2, y1, y2, 'VariableNames',colNames2);
writetable(Output2,[file_name_base, 'Parameters', '.xlsx']); %writes importable file to excel

    disp('Done with file')
    
    
    
end

disp('Done with all files!')
