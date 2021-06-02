# Clear sky Models

These clear sky models can be used to estimate the surface shortwave broadband global/direct/diffuse irradiance at a given location and time. The models varii strongly in complexity and requirements of input. 
This package is a python converter of clear sky modules coded in R acquired from Jamie Brights GitHub: https://github.com/JamieMBright/clear-sky-models

# Requirements:
 - numpy
 - rpy2
 - matplotlib (only example.py)
 example will run with using the conda environment:
 `conda env create -f environment.yml`

# Install:
` pip install git+https://github.com/jonas-witthuhn/clear-sky-models.git#egg=clear_sky_models `

# Usage:
As an example check example.py

All models expect an input of numpy.datetime64 date and solar zenith angle (sza) as float or numpy.array. The sza serves as a control for the time and latitude/longitude coordinates, which are not explicitely required as input. 
For calculation of sza, see: https://gitea.tropos.de/deneke/trosat-base (sunpos)

Further inputs (e.g. altitude, pressure, ozone...) are expected as float or array with the same shape as sza. See models.py for individual model requirements.



# Forked from JamieMBright

---

# Welcome to the clear-sky irradiance model library
This Github repository contains many clear-sky irradiance models as coded in R. Occasionally a model will be available in Matlab, though this is not an objective of ours.
This repository was created to coincide with our research publication titled "Worldwide performance assessment of 75 global clear-sky irradiance models using Principal Component Analysis" published in the Journal of Renewable and Sustainable Energy Reviews and authored by Xixi Sun, Jamie M. Bright, Christian A. Gueymard, Brendan Acord, Peng Wang and Nicholas A. Engerer.

The R code available in this repository was written by Xixi Sun.

## About the models
The reader is referred to our publication and its accompanying supplementary material in order to fully understand the workings of these models. There you will find our interpretation and links to the original work.
PLACEHOLDER URL TO JOURNAL PAPER

The models all require different input data


|No. |Clear-sky Model|Eo|zen|h|alb|p|T|TL|aod|alp|beta|O3|NO2 |H2O|tau|Tot 
|---|------|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:
|1 | TJ || •| | | | | | | | | | | | | 1 
|2 | Schulze || •| | | | | | | | | | | | | 1 
|3 | DPP || •| | | | | | | | | | | | | 1 
|4 | Adnot || •| | | | | | | | | | | | | 1 
|5 | Biga ||•| | | | | | | | | | | | | 1 
|6 | ASHRAE ||•| | | | | | | | | | | | | 1 
|7 | Sharma | •|•| | | | | | | | | | | | | 2 
|8 | El Mghouchi | •|•| | | | | | | | | | | | | 2 
|9 | Yang & Walsh | •|•| | | | | | | | | | | | | 2 
|10 | HLJ |•|•|• ||||||| | | ||| 3 
|11 | Kumar |•|• |||• ||||||| | || 3 
|12 | Campbell |•|• |||• ||||||| | || 3 
|13 | Fu \| Rich | • |• |•|| | |||| | | ||| 3 
|14 | Atwater & Ball-1 |•|•|| |•| | | | | | | |• | |4 
|15 | KASM |•|•| | |•| | | | | | | |• | |4 
|16 | Capderou |•|•| •| |•| | | | | | | | | |4 
|17-22 | Kasten | •|• |• ||||• || | | | | | | 4 
|23-28 | Heliosat-1 | •|•|||• ||• ||| || | | | 4 
|29-34 | ESRA | •|• |•| | | |• | | | | | | | | 4 
|35-40 | Heliosat-2 |•|• |• ||||• ||| | | || | 4 
|41-46 | Ineichen & Perez | •|•|• ||||• || | | ||| | 4 
|47 | CLS |•|•||•|• ||| |||||• | |5 
|48 | King |•|•||•|• ||| ||•||| | |5 
|49 | Josefsson | •|•||•|• ||| |||||• | |5 
|50 | Badescu | • |• |||• ||| |||•||• | |5 
|51 | Simplified Solis | •|•|||• |||• |||||• | |5 
|52 | Advanced Solis | •|•|||• |||• |||||• | |5 
|53 | Perrin | •|•|||• ||| ||•|•||• | |6 
|54 | CEM | •|•||•|• ||| • | ||||• | |6 
|55 | Atwater & Ball-2 | •|• ||• |•|| |•|||||•| |6 
|56 | RSC | •|•||•|•|| | ||•|||•| |6 
|57 | PSIM | •|•||•|•|| | ||•|||•| |6 
|58 | Bashahu |•|•| | |• | | | |•|•| ||•| |6 
|59 | MMAC | •|•||•|•|||•|||||•| |6 
|60 | Yang | •|•|||•|||||•|•||•| |6 
|61 | Calinoiu |•|•||||||||•|•|•|• | |6 
|62 | Hoyt | •|•||•|•|||||•|•||•| | 7 
|63 | MAC | •|•||•|•||||•|•|||•||7 
|64 | METSTAT | •|•| |•|•|||•| | |•| |•| | 7 
|65 | PR | •|•| |•|•| | | | |•|•| |•| | 7 
|66 | Paulescu & Schlett | •|•|•| |•| | | | |•|•| |•| | 7 
|67 | MRM v5 | •|•||•|•|||||•|•||• | |7 
|68 | MRM v6.1 | •|•||•|•|||•| | |•||•|| 7 
|69 | Janjai |•|•|•| |||||•|• |• ||•| | 7 
|70 | Bird | •|•||• |•||| |•|• |• ||• | |8 
|71 | Iqbal-c | •|•||•|•||| |•|• |• ||• | |8 
|72 | Modified Iqbal-c |•|•|•|•| ||||•|•|•||•| |8 
|73 | REST2v5 | •|•||•|•||||• |• |• |• |• || 9 
|74 | REST2v9.1 |•|•||•|•||||•|•|•| |•|•| 9 
|75 | McClear | •|•|•|•|•|•|| •|•| | •| |•| | 10 


We obtained these input data from the MERRA2 reanalysis database. 
The reader is again referred to the data section of our paper in order to learn how to access these inputs.

## Disclaimer
We do not offer support for these models, however, we welcome suggested fixes or edits by users. 
The intellectual property of these models remains with the authors credited within the paper. 
All models made available are already in the public domain, notable methodologies are not such as REST2v9 which is proprietary and McClear which is a website.
