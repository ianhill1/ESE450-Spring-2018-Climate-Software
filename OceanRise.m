function [dH_dt,H,t] = OceanRise(t_years,feedback_initial,feedback_final,atmDepth,oceanDepth)

%parameters for ocean rise model
%a=  5.56;%sea level sensitivity {mm/yr-K)
a=  16.4;%sea level sensitivity {mm/yr-K)
a = a*(1/365.25)*(1/24)*(1/3600);%mm/yr-K *(1 yr/365.25 days)*(1 day/24 hrs)*(1 hr/3600 s) = mm/s-K
%b = -66;%fast response term {mm/K}
b = -49;%fast response term {mm/K}
%Tind = -0.43;%pre industrial equilibrium temperature 
Tind = -0.375;%pre industrial equilibrium temperature

[Ta,To,t,dTa_dt] = BoxModel(t_years,feedback_initial,feedback_final,atmDepth,oceanDepth);

%initial values;
%keep initial values set to zero (see reason in BoxModel)
H0 = 0;
dH_dt0 = 0;

H = [H0];
dH_dt = [dH_dt0];
%convert time to seconds for solving ode's, such that we're operating in basic SI units: 
tspan = t_years*(365.25/1)*(24/1)*(3600/1);%units of coefficients: 365.25 days/1 yr * 24 hrs/1 day * 3600 s / 1 hr 
dt = 10000;%time step (s)
t = 0:dt:tspan;
n = length(t)-1;% #loops to solve ode's = #time steps -1
                % the "-1" is to account for fact that
                % temp arrays already contain one entry (the initial
                % values)
             

for i = 1:n 
    dH_dt(i+1) = a*(Ta(i) - Tind) + b*dTa_dt;
    H(i+1) = H(i) + dt*(a*(Ta(i) - Tind) + b*dTa_dt);
end

dH_dt = dH_dt*3600*24*365.25;%convert back to mm/yr after calculations in mm/s
t = t*(1/3600)*(1/24)*(1/365.25);%{s}*(1hr/3600s)*(1 day/24 hrs)*(1 yr/365.25 days) = years

end
