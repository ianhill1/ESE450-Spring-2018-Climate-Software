function [dH_dt,H,t] = OceanRise(t_years,feedback_initial,feedback_final,atmDepth,oceanDepth,alphaType,alphaJumpTime)

%parameters for ocean rise model
a=  5.56;%sea level sensitivity {mm/yr-K)
a = a*(1/365.25)*(1/24)*(1/3600);%mm/yr-K *(1 yr/365.25 days)*(1 day/24 hrs)*(1 hr/3600 s) = mm/s-K
b = -49;%fast response term {mm/K}
Tind = -0.375;%pre industrial equilibrium temperature

secondsPerYear = (365.25/1)*(24/1)*(3600/1);%{years}*(365.25 days/yr) * (24 hrs/day) * (3600 s /hr) = {seconds} 
YearsPerSecond = (1/3600)*(1/24)*(1/365.25);%{s}*(1hr/3600s)*(1 day/24 hrs)*(1 yr/365.25 days) = {years}
[Ta,~,t,dTa_dt,~,dt] = BoxModel(t_years,feedback_initial,feedback_final,atmDepth,oceanDepth,alphaType,alphaJumpTime);

%initial values;
%keep initial values set to zero (see reason in BoxModel)
H0 = 0;
dH_dt0 = 0;

H = [H0];
dH_dt = [dH_dt0];
%convert time to seconds for solving ode's, such that we're operating in basic SI units: 
t= t*secondsPerYear;%since BoxModel returns time in years
n = length(t)-1;% #loops to solve ode's = #time steps -1
                % the "-1" is to account for fact that
                % temp arrays already contain one entry (the initial
                % values)
             

for i = 1:n 
    dH_dt(i+1) = a*(Ta(i) - Tind) + b*dTa_dt(i);
    H(i+1) = H(i) + dt*(a*(Ta(i) - Tind) + b*dTa_dt(i));
end

dH_dt = dH_dt*3600*24*365.25;%convert back to mm/yr after calculations in mm/s
t = t*YearsPerSecond;%{s}*YearsPerSecond = years

end
