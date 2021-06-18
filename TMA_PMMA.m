function mass = TMA_PMMA(t, l, df, sc, pc, hd, k)
%
% Input:
%   t: time points at which a solution is requested
%   l: Thickness (e.g. 480nm)
%   df: Diffusivity (e.g. 2.2e-10 cm^2/s)
%   sc: Surface Concentration (e.g. 3.708e-3 mol/cm^3)
%   pc: Accessible Polymer Concentration (e.g. 5.758e-3 mol/cm^3)
%   hd: Hindering Factor (e.g. 1200 cm^3/mol)
%   k: Reaction Rate (e.g. 1 cm^3/mol s)
% Output:
%   mass: mass uptake over time
%
    
    pde_options = [];
    x_dim = 1001;
    x = linspace(-l/2, l/2, x_dim); 
    u = pdepe(0, @pdefun, @pdeic, @pdebc, x, t, pde_options, df, sc, pc, hd, k);
    mass = mean((u(:,:,1)+u(:,:,2)), 2) .* l;
    % mass = reshape(mean(u,2), size(t,2), 3) .* l;
    % mass = u .* l / x_dim;
    mass = mass.*(72*10^9)./2; % adjust the unit
end

%------------------------------------------------------------
% PDE initial conditions
% u = [C,S,B]
% C + B -> S
% C: Concentration of free diffusing species
% B: Concentration of polymer functional groups;
% S: Concentration of immobilized product;
function u0 = pdeic(x, df, sc, pc, hd, k)
    u0=[0;0;pc];
end

%------------------------------------------------------------
% PDE components
function [c,f,s] = pdefun(x, t, u, dudx, df, sc, pc, hd, k)
    c=[1;1;1];
    f=[df*exp(-hd*u(2));0;0].*dudx;
    s=[-k*u(1)*u(3);k*u(1)*u(3);-k*u(1)*u(3)];
end

%------------------------------------------------------------
% PDE boundary conditions
function [pl,ql,pr,qr] = pdebc(xl, ul, xr, ur, t, df, sc, pc, hd, k)
    if t <= 62500
        pl=[ul(1)-sc;ul(2);ul(3)];
        ql=[0;0;0];
        pr=[ur(1)-sc;ur(2);ur(3)];
        qr=[0;0;0];
        
    elseif 62500<t && t<62560
        pl=[ul(1)-(sc-sc/60*(t-62500));ul(2);ul(3)];
        ql=[0;0;0];
        pr=[ur(1)-(sc-sc/60*(t-62500));ur(2);ur(3)];
        qr=[0;0;0];

    else
        pl=[ul(1);ul(2);ul(3)];
        ql=[0;0;0];
        pr=[ur(1);ur(2);ur(3)];
        qr=[0;0;0];
    end
end