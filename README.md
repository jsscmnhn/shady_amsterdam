# Shady Amsterdam
This repository is created for the MsC Geomatics Synthesis Project 2024, for the clients Gemeente Amsterdam and MIT Senseable City Lab

The documentation can best be viewed on the [webpage](https://jsscmnhn.github.io/shady_amsterdam/).

The code in this repository makes it possible to create shade maps of any area in the Netherlands for any date and time,
use these shade maps and other datasets to identify and score cool spaces in the region, and to create a walkshed and
routing shady and fast routes.

The repository includes three steps: calculating and creating the shade maps, identifying and scoring the cool spaces and routing. 

## Shade Map (`/shade_calculation`)
This step 

For more information about the workings of the code, look at the page:  [Shade Map Calculation](docs/Shade-Map-Calculation.md)


## Cool Spaces (`/cool_place`)


For more information about the workings of the code, look at the page:  [Cool Space Process](docs/Cool-Spaces.md)

## Network (`/PedestrianNetwork`)

For more information about the workings of this code, look at the page:  [Network Process](docs/Network.md)

# Running the code 
The code can also run for all steps simultaneously.
For this use `main.py`

## API
Find here the [API](docs/api.md)

## Configuration file 
To run the code, configuration files have to be given as input. 
Find out how to do this at  [Setting up the Configuration File](docs/Configuration-setup.md)

## Requirements
The file containing the required packages for this repository can be found at

for a conda environment and for pip install. We recommend using a Conda environment as this simplifies installing GDAL. 

The required datasets for running Analysis Amsterdam region can be found here [Datasets](https://drive.google.com/drive/folders/1LsNp03WkUEMMzGZZci4n8d7l7EE5ZUVt)

## Example run (`/example_run`)
The input, outputs and configuration file settings of an example run of a small region in Amsterdam can be found 