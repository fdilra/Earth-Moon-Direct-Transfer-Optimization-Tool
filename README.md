## Introduction

!! This tool is a work in progress !!

Given an initial guess for the departure date, orbital elements of the transfer orbit, and given the desired lunar orbit height and inclination, this script finds and optimizes a possible trajectory from an Earth circular parking orbit to a circular lunar orbit.

Targeting and delta v optimization are both performed using the fmincon function.

To run this script the following kernels are required:
- de440.bsp
- gm_de440.tpc
- naif0012.tls
- pck00011.tpc

The kernels can be downloaded at https://naif.jpl.nasa.gov/pub/naif/generic_kernels/ and have to be stored in a folder named \kernels that has to be in the same directory as the script and all the functions.

![untitled](https://github.com/fdilra/Earth-Moon-Direct-Transfer-Optimization-Tool/assets/153423123/c57fed7b-4728-4f1b-bc22-560d3877d522)

## Model

The current model used by the script is a Sun-Earth-Moon-spacecraft four-body problem model that considers only point-mass gravity.

## Future Improvements

Future improvements will include an improved model with:
- higher order gravity perturbations
- solar radiation pressure
- higher precision
