# Reaction-Diffusion Model 

This repository provides the MATLAB implementation of the reaction-diffusion model for the vapor phase infiltration (VPI) process in our paper "[Reaction−Diffusion Transport Model to Predict Precursor Uptake and Spatial Distribution in Vapor-Phase Infiltration Processes][1]". 

`TMA_PMMA.m` contains the source code for the numerically solving the partial differential equations. `experiment.m` contains the code for inputting the required variables and then call the numerical solver in `TMA_PMMA.m`. Here is an example from `experiment.m` for using the code.

```
% time discretization 
t1 = linspace(0,5000,2000);
t2 = linspace(5000,125000,4000);
t = [t1 t2(2:end)];
% 8.7 Torr, 130 C, 483 nm
% define variables
l = 4.83E-5*2; % polymer thickness (cm)
df = 1.65E-10; % initial diffusivity (cm^2/s)
sc = 4.436E-3; % surface concentration (mol/cm^3)
pc = 5.656E-3; % polymer concentration (mol/cm^3)
hd = 1150; % hindering factor (cm^3/mol)
k = 1; % reaction rate (cm^3/mol s)
% call numerical solver
mass = TMA_PMMA(t, l, df, sc, pc, hd, k);
```
[1]:https://pubs.acs.org/doi/10.1021/acs.chemmater.1c01283