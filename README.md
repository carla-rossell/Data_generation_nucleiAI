# 1. Quality Control Pipeline
This software performs a Quality Control on FFPE scans. The first step is the division of a scan into 10% overlapping smaller fields of view (FOVs) of 512x512 pixels. **Missing: Add Position within full scan**.
Afterwards, these FOVs go through a quality control check up to make sure only good images are being selected for further analysis. 
## Background Filter
With *filter_out_background.m* we make sure we are filtering out images that don't contain tissue or contain just stroma. The script measures the global intensity of an image or FOV and if it doesn't reach the 1000 intensity threshold it gets filtered out. 
Select the folder where the FOVs are located
```MATLAB
selpath_source=uigetdir('','Source Directory');
files=dir([selpath_source '/i*.tif']);
names_files={files(:).name};
```
Measure the intensities of each FOV in the folder:
```MATLAB
for i=1:numel(files)
    %Load Image
    test_image=imread(char(names_files(i)));
    %If the intensity isn't higher than 1000, condition=1 (The image will be filtered out)
    if mean(test_image(:))<1000
        condition=1;
    else
        condition=0;
    end
   conditions_background(i)=condition;
   clear condition
end
```
The ouput is variable that contains an array named *conditions_background*

## Blob Check and Filtering
Making sure images with large objects (that in the case of our datasets will be artifacts) get filtered out. The threshold used in this case is filtering out images that have a BLOB score larger than the mean of the dataset + 2,5 times the standard deviation.
To measure the blob score the script *blob.m* has to be saved in the same folder where the images are stored. The output will be a csv file (*blob.csv*).

```MATLAB
blob
```
The csv file should then be opened in order to perform the filtering
```MATLAB
blob=import_blob('blob.csv');
blob_score= blob.Blob_score;
a=1;
for y=1:numel(files)
   if conditions_background(y)==0
       blob_score_filtered(a)=blob_score(y);
       names_filtered_blob{a}=names_files{y};
       a=a+1;
   end
end
%Calculate the threshold
threshold=mean(blob_score_filtered)+2.5*std(blob_score_filtered);
%Count how many FOVs will be filtered out
which_blob_true=blob_score_filtered>mean(blob_score_filtered)+2.5*std(blob_score_filtered);
sum(which_blob_true)
```
## Focus Check and Filtering
Applying the focus measuring script we filter out images that have a score larger than the 0.7 threshold, meaning that they are out of focus. The script should be run in the following way:

Iterate over the dataset and measure focus on each FOV:
```MATLAB
for i=1:numel(files)
 I=imread(files(i).name);
 focus_score{i}=focus_minimal(I);
end
```
Filter out the FOVs that have a value higher than the threshold:
```MATLAB
b=1;
a=1;
for y=1:numel(files)
   if conditions_background(y)==0
       blob_score_filtered(a)=blob_score(y);
       a=a+1;
       %Filter out images that have blobs
       if blob_score(y)<threshold==1
           %Perform a extra filtering step to filter out images with a focus score higher than 0.7 
           if focus_score{1,y}<0.7==1
               filtered_images(b)=focus_score{1,y};
               filtered_names{b}=names_files{y};
               b=b+1;
           end
       end
   end
end
```

## Quality Control Output
After performing the quality control the images that pass the QC will be saved in a folder named *Filtered_dataset*.
 ```MATLAB
psource = '/Volumes/BigMichi/Digital_pathology/Lung&Brain_cancer/Splitted_images_512/2/'
pdest   = '/Volumes/BigMichi/Digital_pathology/Lung&Brain_cancer/Splitted_images_512/2/Filtered_dataset/';
pattern = 'i*.tif';
sourceDir = dir(fullfile(psource, pattern));
for k = 1:numel(sourceDir)
    Index = find(contains(filtered_names,sourceDir(k).name));
    if isempty(Index)==0
    sourceFile = fullfile(psource, sourceDir(k).name);
    destFile   = fullfile(pdest, sourceDir(k).name);
    copyfile(sourceFile, destFile);
    end
end
save('workspace_values.mat');
 ```
# 2. Ilastik Training
30 images that have passed the quality control are selected to manually perform training with Ilastik (selecting nuclei vs background). Afterwards prediction is applied to the rest of the samples in the dataset to obtain a segmentation mask of each FOV.
# 3. Watershed Segmentation
A watershed algorithm is applied to the segmentation mask in order to have a more precise segmentation of the nuclei. 
In MATLAB the code is the following. We have to select both the source folder (**segmented_#**) and the output folder (**watershed_#**)
```MATLAB
 %Define a source and a target directory

selpath_source=uigetdir('','Source Directory');
selpath_target=uigetdir('','Target Directory');

%Read all the objects in the source directory and delete the first 3 rows
%of the array

list_files_input=dir([selpath_source '/*.h5']);
% list_files_input=list_files_input(3:end);

%Start the loop through all the .h5 files - first the loop reads the .h5
%files, then it reshape them based on their original format, then it rotate
%the images of 90° left, then it flips vertically the images

for k=1:length(list_files_input)
    %image1 =hdf5info('IDC12879_001_Segmented_mask.h5');
    file_name=strcat([list_files_input(k).folder '/' list_files_input(k).name]);
    image1 = hdf5read(file_name,'exported_data');
    %image1=reshape(image1, [width heigth ]
    image_size=size(image1);
    width=image_size(3);
    height=image_size(2);
  image1 = reshape(image1, [width height]);
    image1 = imrotate(image1,90);
    image1 = flipud(image1);
    image1(image1==510)=0;
    image1(image1==255)=1;
    image1=logical(image1);
%   imshow(image1);
   

%To remove small objects from the image
bw2 = bwareaopen(image1,20);
image1=bw2;
%imshow(bw2);

%Fill holes in the image
bw2=imfill(image1,'holes');
image1=bw2;
%imshow(bw2);

%Distance transform
D = -bwdist(~image1);
%imshow(D,[]);

%Run watershed on the masks

Ld = watershed(D);
%bw2 = image1;
bw2 (Ld == 0) = 0;
%imshow(bw2)

mask = imextendedmin(D,0.75);
%imshowpair(image1,mask,'blend')

D2 = imimposemin(D,mask);
Ld2 = watershed(D2);
bw3 = image1;
bw3(Ld2 == 0) = 0;
%check the watershed
%imshow(bw3)

%Save files as HDF5 format
%timage = reshape(uint8(bw3),[1, height, width]);
timage = reshape(uint8(bw3),[1, height, width]);
%h5create([selpath_target,'/',list_files_input(k).name], '/exported_watershed_masks',[1 height width],'Datatype','uint8');
h5create([selpath_target,'/',list_files_input(k).name], '/exported_watershed_masks',[1 height width],'Datatype','uint8');
h5write([selpath_target,'/',list_files_input(k).name], '/exported_watershed_masks', timage);

end         
```
#4. Nuclei Properties Extraction into .txt table
This part of the Pipeline extracts morphometric properties (**area, perimeter, eccentricity, circularity**) of the nuclei as well as the **total intensity and mean intensity**.

# Procedure
After having both the watershed masks and the .tif tiles we run the script to extract nuclei properties. For the morphometric properties, the script relies only on the information provided by the watersheded binary mask (.h5). For the intensity, both information from the watershed mask and the DAPI staining tile (.tif) will be used.
```
%Select the path where the watershed images are located
selpath_source=uigetdir('','Source Directory');

%Read all the objects in the source directory and delete the first 3 rows
%of the array
pattern = 'i*';
files= dir(fullfile(selpath_source, pattern));

%Iterate over the images to extract the values of the nuclei
for i=1:numel(files)
    %Open Ilastik output (file obtained from Michele)
    %data=hdf5read(files(i).name,'exported_data');
   file_name=strcat([files(i).folder '/' files(i).name]);
    %data=hdf5read(files(i).name,'exported_data');
    data=hdf5read(file_name,'/exported_watershed_masks');
    clearvars file_name
    image_size=size(data);
    width=image_size(3);
    height=image_size(2);
    im= reshape(data, [height width]);
    %im=im';
    %Binarize image
    im(im==1)=1;
    im(im==2)=0;
    im=logical(im);
    %Find nuclei boundaries
    %nuclei_boundaries=bwboundaries(im,'nohole');
    nuclei_boundaries=bwboundaries(im);
    %save the boundaries to plot afterwards the masks
    nuclei_boundaries_all{i}=nuclei_boundaries;
    %Filter out those images that have no cells
    if isempty(nuclei_boundaries)
        continue
    else
    
    %Find nuclei properties
    %nuclei_prop=regionprops(im,'Area','Centroid','ConvexHull');
    nuclei_prop=regionprops(im,'Area','Centroid','Perimeter','MajorAxisLength','MinorAxisLength','Eccentricity','ConvexHull');
    properties_all_images{i,1}=nuclei_prop;
    clear vars  ratio_all_cells area_all_cells perimeter_all_cells  nuclei_prop properties_all_cells
    end
end
```
Put properties in a table

```
%This script should be run after the properties are extracted 
%Iterate over the properties to create a table of all_nuclei
%Empty arrays for area and circularity
all_areas=[]
all_circ=[]
all_perimeter=[]
for i=1:numel(files)
    list_area_fov=[];
    %Concatenate areas in a list
    list_area_fov=vertcat(properties_all_images{i,:}.Area);
    %Concatenate perimeter in a list
    list_per_fov=vertcat(properties_all_images{i,:}.Perimeter);
    %Calculate the formula for circularity
    explist_per=list_per_fov.^2;
    circularity=4*pi*(list_area_fov./explist_per);
    %Append to list
    all_areas=vertcat(all_areas,list_area_fov);
    all_circ=vertcat(all_circ,circularity);
    all_perimeter=vertcat(all_perimeter,list_per_fov);
end

%Get values for eccentricity
all_eccentricity=[]
for i=1:numel(files)
    list_area_fov=[];
    list_area_fov=vertcat(properties_all_images{i,:}.Eccentricity);
    all_eccentricity=vertcat(all_eccentricity,list_area_fov);
end

%Get x and y coordinates for each nuclei
all_xcoord=[]
all_ycoord=[]
im_names=[]
for i=1:numel(files)
    [a,other_col]=strtok(files(i).name,'c');
    image_name=files(i).name;
    other_col=erase(other_col,".h5");
    column=erase(other_col,"c");
    column=str2double(column);
    [a,other_row]=strtok(files(i).name,'r');
    [other_row,a]=strtok(other_row,'_');
    row=erase(other_row,"r");
    row=str2double(row);
    add_value_x=512*column;
    add_value_y=512*(row+1);
    list_coord_fov=[];
    list_coord_fov=vertcat(properties_all_images{i,:}.Centroid);
    xcoord=list_coord_fov(:,1);
    xcoord=xcoord+add_value_x;
    all_xcoord=vertcat(all_xcoord,xcoord);
    ycoord=list_coord_fov(:,2);
    ycoord=add_value_y-ycoord;
    all_ycoord=vertcat(all_ycoord,ycoord);

end

%Extract the names of each FOV
all_names=[]
for i=1:numel(files)
  [a,name]=strtok(files(i).name,'.');   
  name=erase(name,".h5");
  name=erase(name,"._");
  %[name]*numel(properties_all_images{i,1});
  test=repmat({name},numel(properties_all_images{i,1}),1);
  all_names=vertcat(all_names,test);
end
all_name_mat = string(all_names);


```

Get intensities
Select folder where Tiff files are located: **Filtered_#** 
Select folder where Watershed files are located: **watershed_#**
```
%Current folder should be patient's folder
%Get pixel intensities
%Read .tif file
selpath_tiff=uigetdir('','Tiff Filtered Directory');
selpath_h5=uigetdir('','Watershed Directory');
%Read all the objects in the source directory and delete the first 3 rows
%of the array
listing_tif=dir([selpath_tiff '/*.tif']);
listing_h5=dir([selpath_h5 '/*.h5']);

%get files from get_table

files.name
%for h5 files use this value
%for tif files delete

for i=1:numel(listing_h5)
current_file=erase(files(i).name,'.h5'); 
file_name=fullfile(listing_tif(i).folder, strcat(current_file,'.tif'));
splitted_image=imread(file_name,'tif'); %change command?
file_name2=fullfile(listing_h5(i).folder, files(i).name);
h5_image=hdf5read(file_name2,'/exported_watershed_masks');
im=reshape(h5_image,512,512);
    %Binarize image
    im(im==1)=1;
    im(im==2)=0;
    im=logical(im);
%imshow(a)
all_intensities=regionprops(im,splitted_image,'MeanIntensity');
pixel_values=regionprops(im,splitted_image,'PixelValues');
pixcell=struct2cell(pixel_values);
intensities_all_images{i,1}=all_intensities;
cell_intensities{i,1}= cellfun(@sum,pixcell);
end

all_mean_intensities=[]
all_areas=[]
for i=1:numel(intensities_all_images)
    list_int_fov=[];
    %Concatenate intensities in a list
    list_int_fov=vertcat(intensities_all_images{i,:}.MeanIntensity);
    list_pixels_fov=vertcat(properties_all_images{i,:}.Area);
    
    %Append to list
    all_mean_intensities=vertcat(all_mean_intensities,list_int_fov);
end

total_intensities=[]
for i=1:numel(cell_intensities)
    list_total_int=vertcat(cell_intensities{i,:});
    total_intensities=vertcat(total_intensities,list_total_int');
end

%Add intensity 

%Join all the variables into an array and convert it to a table
general_table=[all_name_mat,all_xcoord,all_ycoord,all_areas,all_perimeter,all_eccentricity, all_circ, all_mean_intensities, total_intensities];
general_table=array2table(general_table);

writetable(general_table,'nuclei_Patient.txt','Delimiter','\t')
```
