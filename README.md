# regularized optimal mass transport (rOMT) applied in Glymphatic-Lymphatic fluid flows 
This project includes code for <br />
(1) rOMT algorithm <br />
(2) Lagrangian representation of Glymphatic Dynamics (GLaD) analysis <br />
(3) Neighborhood-based Cosine Analysis (NCA) <br />
where (1) runs the main rOMT model on the dataset and (2-3) post-processes the results from (1).

## System Requirements
The code was mainly written and ran in Matlab (R2018a for rOMT algorithm and R2019b for post-processing), with a small section in GLaD analysis ran with Python 3.7.3. 

### Dependencies
#### Matlab
(a) NIfTI_analyze https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image for loading and exporting nifty files <br />
(b) vtkwrite https://www.mathworks.com/matlabcentral/fileexchange/47814-vtkwrite-exports-various-2d-3d-data-to-paraview-in-vtk-file-format for visualization in vtk format <br />
#### Python
numpy, scipy, dipy
#### recommended software for visualization
(a) Amira 6.5.0 for speed map
(b) VisIt 3.0.2 for flux vectors
## Demo
The dataset is too large to run on a typical desktop computer, so usually we run on a cluster with 40 cores. However, for the purpose of demostration, we resampled the original sample data by 0.5 and reduced the data frames included. The paramters for running on the original large dataset can be found in getParams_original.m and the dataset can be downloaded at xxx.

