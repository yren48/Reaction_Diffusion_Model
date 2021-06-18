# Reaction-Diffusion Model 

This repository provides the MATLAB implementation of the reaction-diffusion model for the vapor phase infiltration (VPI) process from Ren et al. (2021). 

The model is characterized by the following system of the partial differential equations,

$$
\left\{\begin{array}{l}
    \frac{\partial C_{free}}{\partial t} = D\frac{\partial^2 C_{free}}{\partial x^2} - kC_{free}C_{polymer} \\
    \frac{\partial C_{product}}{\partial t} = kC_{free}C_{polymer} \\
    D = D_0 \exp(-K^{\prime}C_{product}) \\
    \frac{\partial C_{polymer}}{\partial t} = -kC_{free}C_{polymer} \\
  \end{array}\right.
$$

with the following initial and boundary conditions,

$$
\left\{\begin{array}{lll}
    C_{free} = 0, & 0 < x < l, & t = 0 \\
    C_{product} = 0, & 0 < x < l, & t = 0 \\
    C_{polymer} = C_{polymer}^{0}, & 0 < x < l, & t = 0 \\
    \frac{\partial C_{free}}{\partial x} = 0, & x = 0, & t > 0 \\
    C_{free} = C_{s}, & x = l, & t > 0 \\
  \end{array}\right.
$$

where $C_{free} (\text{mol}/\text{cm}^3)$ is the concentration of the free diffusing vapor-phase precursor, $C_{polymer} (\text{mol}/\text{cm}^3)$ is the concentration of the accessible reactive polymeric functional groups, $C_{product} (\text{mol}/\text{cm}^3)$ is the concentration of immobilized product from the reaction between the free diffusing vapor-phase precursor and the polymeric functional groups. There are five parameters $\theta = \{D_{0},C_{s},C_{polymer}^{0},K^{\prime},k\}$ that we need to provide for the MATLAB code, where $D_{0} (\text{cm}^2/\text{s})$ is initial diffusivity of the free diffusing vapor-phase precursor, $C_{s} (\text{mol}/\text{cm}^3)$ is the surface concentration of the free diffusing vapor-phase precursor, $C_{polymer}^{0} (\text{mol}/\text{cm}^3)$ is the initial concentration of accessible reactive polymeric functional groups, $K^{\prime} (\text{cm}^3/\text{mol})$ is the hindering factor describing how immobilized product $C_{product}$ slows down the diffusivity of free diffusing vapor, and $k (\text{cm}^3/\text{mol}\cdot\text{s})$ is the associated reaction rate. There are three more operational parameters, polymer thickness $l$, temperature, and vapor pressure, that we can control in the experiment, and the polymer thickness $l$ are also an input for the MATLAB code.

`TMA_PMMA.m` contains the source code for the numerically solving the partial differential equations. `experiment.m` contains the code for inputting the required variables and then call the numerical solver in `TMA_PMMA.m`. Here is an example from `experiment.m` for using the code.

```
% time discretization 
t1 = linspace(0,5000,2000);
t2 = linspace(5000,125000,4000);
t = [t1 t2(2:end)];
% 8.7 Torr, 130 C, 483 nm
% define variables
l = 4.83E-5*2; % polymer thickness (cm)
df = 1.65E-10; % initial diffusivity, D_0 (cm^2/s)
sc = 4.436E-3; % surface concentration, C_s (mol/cm^3)
pc = 5.656E-3; % polymer concentration, C_polymer^0 (mol/cm^3)
hd = 1150; % hindering factor, K' (cm^3/mol)
k = 1; % reaction rate, k (cm^3/mol s)
% call numerical solver
mass = TMA_PMMA(t, l, df, sc, pc, hd, k);
```