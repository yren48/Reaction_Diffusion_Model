clc;
clear;
close all;

%   t: time points at which a solution is requested
%   l: polymer thickness (cm)
%   df: initial diffusivity of the free diffusing 
%       vapor-phase precursor (cm^2/s)
%   sc: surface concentration of the free diffusing 
%       vapor-phase precurosr (mol/cm^3)
%   pc: initial concentration of accessible reactive 
%       polymeric functional groups ( mol/cm^3)
%   hd: hindering factor describing how immobilized product slows down 
%       the diffusivity of free diffusing vapor (cm^3/mol)
%   k: associated reaction rate (cm^3/mol s)

t1 = linspace(0,5000,2000);
t2 = linspace(5000,125000,4000);
t = [t1 t2(2:end)];

% 8.7 Torr, 130 C, 483 nm 
% define variables
l = 4.83E-5*2;
df = 1.65E-10;
sc = 4.436E-3;
pc = 5.656E-3;
hd = 1150;
k = 1;
% call pde solver
disp("running the model for 8.7 Torr, 130 C, 483 nm...")
mass = TMA_PMMA(t, l, df, sc, pc, hd, k);
% visualization
h = figure;
drawnow;
hold on
plot(t.^0.5, mass, '.-b', 'MarkerSize', 10);
title("8.7 Torr, 130 C, 483 nm");
hold off

waitfor(h);

% 10.5 Torr, 130 C, 607 nm 
% define variables
l = 6.07E-5*2; % Changed thickness
sc = sc * 10.5 / 8.7; % adjust for the different pressure
% call pde solver
disp("running the model for 10.5 Torr, 130 C, 607 nm...")
mass = TMA_PMMA(t, l, df, sc, pc, hd, k).*0.9621; %Applied the correction factor for small change in surface area
% visualization
h = figure;
drawnow;
hold on
plot(t.^0.5, mass, '.-b', 'MarkerSize', 10);
title("10.5 Torr, 130 C, 483 nm");
hold off
