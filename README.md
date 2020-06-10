# SPART
The SPART model: a soil-plant-atmosphere radiative transfer model for satellite measurements in the solar spectrum
## Introduction
The model uses three computationally efficient RTMs for soil (BSM), vegetation canopies (PROSAIL) and atmosphere (SMAC), respectively. The sub-models are coupled by using the four-stream theory and the adding method. The resulting `Soil-Plant-Atmosphere Radiative Transfer model' (SPART) simulates directional TOA spectral observations, with all major effects included, such as sun-observer geometries and non-Lambertian reflectance of the land surface.

## References
The users are recommended to read the paper that describe the developments and usage of the model. 

Yang, P., van der Tol, C., Yin, T., & Verhoef, W. (2020). The SPART model: A soil-plant-atmosphere radiative transfer model for satellite measurements in the solar spectrum. Remote Sensing of Environment, 247, 111870.

For the details of the radiative transfer modelling 

Yang, P., Verhoef, W., & van der Tol, C. (2017). The mSCOPE model: A simple adaptation to the SCOPE model to describe reflectance, fluorescence and photosynthesis of vertically heterogeneous canopies. Remote sensing of environment, 201, 1-11.
## Prerequisites
The model was tested using **Matlab 2017a**. 
The users are expected to have the **Microsoft** installed. The input data are structured in an excel spreadsheet. 
If there is not Microsoft, it is possible to define the inputs in a script. 

## Run the model 
The users can run **run_SPART.m** or any examples in **run_SPART_some_examples.m** to have a general idea of this model. 
