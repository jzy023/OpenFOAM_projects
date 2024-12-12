# OpenFOAM projects

This repo contains case folders for ESI OpenFOAM.


## Table of Contents
1. [Multiphase boilingBubbles](#multiphase-boilingbubbles)
2. [Multiphase condensateEvaporate](#multiphase-condensateevaporate)
<!-- 3. [Urban testCase](#urban-testcase)
4. [CHAD validation cases](#chad-validation-cases) -->


## Multiphase boilingBubbles
solver: **icoReactingMultiphaseInterFoam**

The solver is based pn continuous-continuous phase interaction modelling (VoF) where the solver solves all $N$ phases and forces a 
sharp interface between phases. 

The surface tension between liquid and gas phase $\sigma$ are set to be $\sigma=0.055$, $\sigma=0.07$ and $\sigma=0.10$ from left to right.

![boilingBubblesGif](https://raw.githubusercontent.com/jzy023/gifs/main/OpenFOAM/multiphase_boilingBubbles/boilingBubbles.gif)


## Multiphase condensateEvaporate
solver: **interCondensatingEvaporatingFoam**

The solver is based on dispersed-continuous phase interaction modelling (Mixture model) where the solver solves $N-1$ phases without forceing a sharp interface between phases. 

The follwoing gif shows 3 case with $\text{cAlpha} = 0.0$ (no sharp interface in ./system/fvSolution) where the condensation and evaporation coefficients ($K_{C}$ and $K_{E}$) are, from left to right, to be:


$K_{C} = 50$, $K_{E} = 50$; 

$K_{C} = 250$, $K_{E} = 250$;

$K_{C} = 1e^{-16}$, $K_{E} = 50$  

<!-- The surface tension between liquid and gas phase $\sigma$ are set to be $\sigma=0.055$, $\sigma=0.07$ and $\sigma=0.10$ from left to right. -->

![condensateEvaperateGif](https://raw.githubusercontent.com/jzy023/gifs/main/OpenFOAM/multiphase_condensateEvaperate/condensateEvaperate.gif)

The follwoing gif shows 3 case with $\text{cAlpha} = 1.0$ (forcing sharp interface in ./system/fvSolution) where the condensation and evaporation coefficients ($K_{C}$ and $K_{E}$) are, from left to right, to be:


$K_{C} = 50$, $K_{E} = 50$; 

$K_{C} = 250$, $K_{E} = 250$;

$K_{C} = 1e^{-16}$, $K_{E} = 50$  

<!-- The surface tension between liquid and gas phase $\sigma$ are set to be $\sigma=0.055$, $\sigma=0.07$ and $\sigma=0.10$ from left to right. -->

![condensateEvaperateGif](https://raw.githubusercontent.com/jzy023/gifs/main/OpenFOAM/multiphase_condensateEvaperate/condensateEvaperate_sharp.gif)

## To Be Continued...




