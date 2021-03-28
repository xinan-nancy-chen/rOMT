# Regularized Optimal Mass Transport (rOMT) Applied in Glymphatic-Lymphatic fluid flows 
This project includes code for <br />
(1) rOMT algorithm <br />
(2) Lagrangian representation of Glymphatic Dynamics (GLaD) analysis <br />
(3) Neighborhood-based Cosine Analysis (NCA) <br />

where (1) runs the main rOMT model on the dataset and (2-3) post-processes the results from (1).

![pipeline](pipeline.png)

## System Requirements
The code was mainly written and ran in Matlab (R2018a for rOMT algorithm and R2019b for post-processing), with a small section in GLaD analysis ran with Python 3.7.3. 

### Dependencies
#### Matlab
(a) NIfTI_analyze https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image for loading and exporting nifty files <br />
(b) vtkwrite https://www.mathworks.com/matlabcentral/fileexchange/47814-vtkwrite-exports-various-2d-3d-data-to-paraview-in-vtk-file-format for exporting mat to vtk format <br />
#### Python
numpy, scipy, dipy
#### Recommended Software for Visualization
(a) Amira 6.5.0 for speed map <br />
(b) VisIt 3.0.2 for flux vectors <br />

## Demo
The original 3D MRI dataset (each of size 100x106x100 and in total 11 frames) is too time and memory-consuming to run on a typical desktop computer, so usually we put it on a CPU cluster with 40 cores. Even though it may lose lots of details and information, for the purpose of demonstration, we downsized the original sample data by 0.5 and reduced the data frames to 7. It takes about 70 minutes to run the sample data with default parameters on a computer with 2.6 GHz Intel Core i7 and 16 GB memory. <br />

### Instructions

To run the rOMT algorithm:<br />
(1) Set model parameters, data directory in ```getParams.m```. You can skip this step by using the default parameters within.<br />
(2) Run ```driver.m```.<br />

You are expected to see output like the following:<br />
```
=============== rOMT Starts ===============
______________________________________________

 tag:		C294
 dataset:	CAA
 sigma:		0.0020
 gamma:		0.0080
 beta:		0.0001
 nt:		10
 dt:		0.40
 pcg:		20
 mask:		1
 do_resize:	1
 resize_factor:	0.50
 start frame:	31
 end frame:	52
 frame jump:	3
```



The paramters and instructions for running on the original large dataset can be found in getParams_original.m and the whole dataset can be available upon request.

