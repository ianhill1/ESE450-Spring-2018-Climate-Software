

function [Ta,To,t,dTa_dt] = BoxModel(t_years,feedback_initial,feedback_final)
OceanDepth = 100;%m
AtmDepth = 10000;%m
Ca = 750*AtmDepth;%J/m^2-K
Co = 3.35e+6*OceanDepth;
cao = 14.2;
Q = 4.0;


%initial values;
To0 = 0;
Ta0 = 0;
H0 = 0;
dH_dt0 = 0;
Ta = [Ta0];
To = [To0];
H = [H0];
dH_dt = [dH_dt0];
tspan = t_years*(365.25/1)*(24/1)*(3600/1);%365.25 days/yr * 24 hrs/day * 3600 s /hr 
dt = 1000;%time step (s)
t = 0:dt:tspan;
n = length(t)-1;% #temp data entries = #time steps -1, to account for fact that
                % temp arrays already contain one entry (the initial
                % values)

f = linspace(feedback_initial,feedback_final,length(t));
for i = 1:n %start @ 2 since first element is initial vals
    Ta(i+1) = Ta(i) + dt*( -(f(i)/Ca + cao/Ca)*Ta(i) + (cao/Ca)*To(i) + Q/Ca);
    To(i+1) = To(i) + dt*( (cao/Co)*Ta(i) - (cao/Co)*To(i) );
    dTa_dt = -(f(i)/Ca + cao/Ca)*Ta(i) + (cao/Ca)*To(i) + Q/Ca;
end

dH_dt = dH_dt*3600*24*365.25;%convert back to mm/yr after calculations in mm/s
t = t*(1/3600)*(1/24)*(1/365.25);%{s}*(1hr/3600s)*(1 day/24 hrs)*(1 yr/365.25 days) = years

end


