# ZDWBM-msnow

Hwang &amp; Devineni (2022) developed a parsimonious snow module based on Budyko's framework and applied it to Zhang's Dynamic Water Balance Model (ZDWBM; Zhang et al., 2008). As part of their work, a monthly parameterization was also proposed to reflect the seasonal behavior of catchment characteristics. The augmented model that includes snow with monthly parameterization (ZDWBM-msnow) showed significant improvements in simulating monthly streamflow. This git repository contains a sample Matlab code for ZDWBM-msnow from Hwang &amp; Devineni (2022) presented in the Water Resources Research. Please cite this paper if you use ZDWBM-msnow.

Hwang, J., &amp; Devineni, N. (2022). An improved Zhang’s Dynamic Water Balance Model using Budyko‐based snow representation for better streamflow predictions. Water Resources Research, e2021WR030203.  
Zhang, L., Potter, N., Hickel, K., Zhang, Y., &amp; Shao, Q. (2008). Water balance modeling over variable time scales based on the Budyko framework–Model development and testing. Journal of Hydrology, 360(1-4), 117-131.

## Description of the Files

This repository contains nine files:

[Discharge Data](4033000.AMM)

[Potential Evapotranspiration Data](4033000.PET)

[Precipitation Data](4033000.PRE)

[Daily Max Temperature Data](4033000.TMAX)

[Daily Min Temperature Data](4033000.TMIN)

[Geographical Data](region4_Smax100.txt)

[Main Code](Sample_code.m)

[ZDWBM-msnow](zhang_model_snow.m)

[Objective Function](zhangModelError.m)

[Nash-Sutcliffe Efficiency](nash_sutcliffe.m)

## Getting started
This sample code estimates monthly runoff for station 4033000 based on ZDWBM-msnow. Please see Hwang & Devineni (2022) for detailed descriptions of the model.

## Requirements
We tested the code using Matlab R2018b.
