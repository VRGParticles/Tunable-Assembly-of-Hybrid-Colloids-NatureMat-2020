# Tunable-Assembly-of-Hybrid-Colloids-NatureMat-2020
This is a repository for the three custom algorithms used in: 

*Tunable Assembly of Hybrid Colloids Induced by Regioselective Depletion*, NatureMat, 2020.  

This repository provides image analysis for micrographs of polymorphic particle/colloid assemblies.  
Particles are classified into positional/rotational order or plane symmetry groups.  The summarized data is found in Figure 4, as well as Supplementary Tables 1-3 of the manuscript.

![Image](https://dl.boxcloud.com/api/2.0/internal_files/677676466958/versions/719227641758/representations/png_paged_2048x2048/content/1.png?access_token=1!6tR27cEQRYmO4oIB67Q89ZlfeFalJkz_SaX7hWMZ3XXtZdYqRFgUdCONctbaY4ObzJk9azLyrarxuy1xSbPiFRvCdQqcGHxolYTKJGvrGeDLaBd1cxi6_1naPjZSHRA17vTJ51Ron09__QiI2v102Gcn_0n8iiZK7CoHInR2o2BZISvmcC5IwHEV0ItBYi6-s3Qanc8Me1Jjwo5WsIh8YZv1NgfNF7cfycGyuLc_Mu_j_uMey2f3cCKpaPnehb58N-pfn4yYUyJtpedaaA4lWXq-KYdUqzZlAnDaHbIW0B3cqdBVcY-_014FAKOEKSz3tI76nDf8ByBgTecvdEJpKotRZ3nh7AduXj9eeldaSq6rlCMdWpKyzy-9WOrVYUo8CHQjqUgXiDMp_Hdv8Q9Q7ilJFm8aNq_X5tSMI0Uur01RhL4EcPZHH7MwNT7YL8BwQy-A_MFzVbzm6jvJBST42Aczo7YDcO6FwWFM3EIA1hzGu5OVXEkekYNEW15Dvgn_cGWih9dg40fkRxc5zwn7CEhMysXSLyyoPU4nCfeFHCsg7kY4BrVYZMvMg9oMPUug9dHr&shared_link=https%3A%2F%2Fnyu.app.box.com%2Fv%2FNM-Regioselective-Depletion&box_client_name=box-content-preview&box_client_version=2.43.0)
# What is Included #
There are three folders, one for each algorithm / particle classifier.


**MATLAB Files (.m):**

*Brick_Amorphous:* Amorphous (positional order only),  Brick (*cmm* plane symmetry group), Shared (boundary of both classes), and Disorder (none of the 3 previous classes).

*Brick_Herringbone:* Herringbone (*pgg* plane symmetry group), Brick (*cmm* plane symmetry group),  Shared (boundary of both classes), and Disorder (none of the 3 previous classes).

*Brick_Square:* Square (*p4m* plane symmetry group), Brick (*cmm* plane symmetry group),  Shared (boundary of both classes), and Disorder (none of the 3 previous classes).


Associated **parameter files** are included to replicate the results found in the manuscript (.csv).

### General Algorithm Overview
These algorithms classify particles based on classification conditions determined from interparticle distance and difference in particle orientation relationships. 
These two parameters are used to create further conditions, such as an amount of neighbors meeting a certain angle criteria or class criteria, to further classify particles.
Each part of the script is sectioned and commented to specify what conditions are being used to classify the particles.  For each classification algorithm,
the conditions are customized to the desired output, for example, herringbone vs. brick particles, and are named accordingly.

## Getting Started
*MATLAB Requirements*

To use these files, you must have MATLAB installed.  Scripts were written using MATLAB R2019b. These algorithms use the Image Processing Toolbox and the Statistics and Machine Learning Toolbox.

*Other MATLAB Prerequisites*

Each script calls upon the function, [bpass](http://www.physics.emory.edu/faculty/weeks/idl/kit/bpass.pro), copyright 1997, John C. Crocker and David G. Grier.

*Input Files*

The raw images/micrographs are located in their [associated folders](https://nyu.app.box.com/v/NM-Regioselective-Depletion/folder/115393668444) as part of Data Availability.
In each of these folders, the 10 input images are labeled in order as found in the Supplementary Tables 1-3.  The resulting labeled figures, as well as a figure legend, are also available.
The obtained labeled images, as well as the output analysis data, can be replicated by downloading these input files and running the file(s) with the associated parameters found in the corresponding repository folder.  The parameter files are labeled for which images they correspond to, if applicable.

## Running the Scripts

All input variables are located at the top of the script, below the 'Input Parameters' heading. 
To run the script, the path/directory of the image must be edited into the 'file_path' input variable, and the image file type into the 'full_file_path' variable. 
All other parameters below these two variables are pre-loaded with one of the parameters found in the script folder repository.  For a new image, these parameters will have to be tweaked until appropriate for the image.
```
Tip: For analysis of your own image (not provided by this repository), use the MATLAB function
'imshow' below each image segmentation step to find reasonable input parameters for your image.
```
When the script is then run, the image(s) are read into the workspace and analyzed.  If a wildcard (*.filetype) is used and there are multiple files of the same type in the directory, the script will loop through these images and analyze each.


### Analysis Outputs
Each algorithm will first open a MATLAB figure, and will then color code the input image into the four particle classes (or cropped image, if input variable 'cropped' is set equal to 1).  
Then, a .png file of this figure is saved to the created Image_Files directory.  The amount of particles in each class is summated and converted to a percent. 
This data, along with the total amount of detected particles, is output as an Excel File into the created Excel_Files directory.
For record-keeping, an Excel file of each of the input parameters is also output per image into the Excel_Files directory.


### More Advanced Script Customization 

Each of these algorithms were created specifically to the input files noted in the *Getting Started* section above.  For more advanced image analysis, 
i.e. customizing these algorithms more toward a user-provided image instead of the provided input files, the user can edit the script
below the 'Input parameters' section.   Each section is commented to provide the user with details as to what the classifier is doing, 
and ultimately what the final classification criteria is for each class.  The user can edit the classification criteria at a step (ex, the structure element, number of neighbors, or the constant 
added to a distance threshold) for a differently customized algorithm.  Similarly, the names of classes or color labels of particles can be changed.
```
Tip: Once the desired parameters are found/testing parameters for the image is no longer needed, the 
figure creation within MATLAB can be suppressed by changing 'figure;' to 'figure('visible', 'off');'
If a labeled figure is not needed, the two figure creation sections can be commented.
If the legend is not desired, the legend section can similarly be commented.
```


### Troubleshooting
*If a user-provided image is used instead of the provided input files*
Input parameters can be customized  by adding 'imshow' functions below each image segmentation step.
```
Tip: To go step-by-step during segmentation, the line 'return'
can be added to stop the script from running at that line. 
```
Each line can also be edited, if need be, to fit the image.
MATLAB related issues can be solved by viewing the [MATLAB documentation.](https://www.mathworks.com/help/index.html)

### Author

* **Veronica Grebe** - [Weck Group, New York University](http://weckresearch.com/home)

For a full list of manuscript authors who contributed to the overall project, see the [manuscript](PLACEHOLDER).

### License

This project is licensed under the GNU GPL-3.0 License - see the [LICENSE](LICENSE) file for details

### What is Next?
Stay tuned for more particle image analysis!
