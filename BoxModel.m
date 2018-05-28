

function [Ta,To,t,dTa_dt] = BoxModel(t_years,alpha_initial,alpha_final,atmDepth,oceanDepth)
%AtmDepth = 10000;%m
%OceanDepth = 1000;%m
Ca = 750*atmDepth;%J/m^2-K
Co = 3.35e+6*oceanDepth;
cao = 14.2;
Q = 4.0;%W/m^2


%initial values
%keep set to zero. The models predict temperature rise, ocean rate rise,
%and sea level rise, not the absolute magnitude of any of these properties.
%temperature. Morover, the model is unlikely to work due to being highly
%sensitive to initial conditions. 
To0 = 0;%2010 avg ocean temp: https://www.neefusa.org/nature/water/warming-ocean 
Ta0 = 0;%2010 avg temp (https://www.ncdc.noaa.gov/sotc/global/201013)
dTa_dt0 = 0;%atm temp i.v.
H0 = 0;
dH_dt0 = 0;
Ta = [Ta0];
To = [To0];
H = [H0];
dH_dt = [dH_dt0];
dTa_dt = [dTa_dt0];
tspan = t_years*(365.25/1)*(24/1)*(3600/1);%365.25 days/yr * 24 hrs/day * 3600 s /hr 
dt = 1000;%time step (s)
t = 0:dt:tspan;
n = length(t)-1;% #temp data entries = #time steps -1, to account for fact that
                % temp arrays already contain one entry (the initial
                % values)

alpha = linspace(alpha_initial,alpha_final,length(t));
for i = 1:n %start @ 2 since first element is initial vals
    %Ta(i+1) = Ta(i) + dt*( -(-f(i)/Ca + cao/Ca)*Ta(i) + (cao/Ca)*To(i) + Q/Ca);
    Ta(i+1) = Ta(i) + dt*( -(-alpha(i)/Ca + cao/Ca)*Ta(i) + (cao/Ca)*To(i) + Q/Ca);
    To(i+1) = To(i) + dt*( (cao/Co)*Ta(i) - (cao/Co)*To(i) );
    dTa_dt = -(-alpha(i)/Ca + cao/Ca)*Ta(i) + (cao/Ca)*To(i) + Q/Ca;
end

dH_dt = dH_dt*3600*24*365.25;%convert back to mm/yr after calculations in mm/s
t = t*(1/3600)*(1/24)*(1/365.25);%{s}*(1hr/3600s)*(1 day/24 hrs)*(1 yr/365.25 days) = years

end


