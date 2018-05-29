clear
clc
clf
%% cfc reaction solver

%elevation
T = 295;%K


%rate constants
k1 = 1e-12;%{s^-1}
k2 = 1e-33;%{cm^6/molecule^2-s}
k3 = 1e-3;%{s^-1}
k4 = (8e-12)*exp(-2060/T);%{cm^3/molecule-second}
k5 = 1e-7;%:{s^-1}
k6 = 2.1e-11;%{cm^3/molecule-s}
k7 = 3.8e-11;%{cm^3/molecule-s}
na = 5.0e+16;%density of air {molecules/cm^3}
%initial values
%O3 = 3e+12;%molecules of O_3/cm^3 (iv for y3)
O = 0.2e+15;%atoms oxygen/cm^3 (iv for y1)
CF2Cl2 = 4e+12;%  iv for y4
Cl = 1e+15; % iv for y5
ClO = 2e+10;%iv for y6
%O2 = O3*k3/(k2*na*O);%iv for y2 -- equation given on pg 3: http://www.columbia.edu/itc/chemistry/chem-c2407/hw/ozone_kinetics.pdf
O2 = 0.2*na;
O3 = O2*na*O*k2/k3;%molecules of O_3/cm^3 (iv for y3)
%O = (2*k1*O2 + k3*O3)/(k2*O2*na + k4*O3);%atoms oxygen/cm^3 (iv for y1)
%ClO2 = 3e+6;
appxmass = O+CF2Cl2 + Cl + ClO + O2 + O3
KineticODE = @(t,y) [2*k1*y(2) - na*k2*y(1)*y(2) + k3*y(3) - k4*y(1)*y(3);
                     na*k2*y(1)*y(2) - k3*y(3) - k4*y(1)*y(3);
                     0;
                     -k5*y(4);
                     k5*y(4) - k6*y(5)*y(3) + k7*y(6)*y(1);
                     k6*y(5)*y(3) - k7*y(6)*y(1)];
                 
iv = [O ; O2 ; O3; CF2Cl2 ; Cl ; ClO];
tspan = [0 1e+6];
[t,y] = ode45(KineticODE , tspan, iv);


                     
                     
                     






