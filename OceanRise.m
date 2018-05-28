function [dH_dt,H,t] = OceanRise(t_years,feedback_initial,feedback_final)

%parameters for ocean rise model
a=  5.56;%sea level sensitivity {mm/yr-K)
a = a*(1/365.25)*(1/24)*(1/3600);%mm/yr-K *(1 yr/365.25 days)*(1 day/24 hrs)*(1 hr/3600 s) = mm/s-K
b = -66;%fast response term {mm/K}
Tind = -0.43;%pre industrial equilibrium temperature 


[Ta,To,t,dTa_dt] = BoxModel(t_years,feedback_initial,feedback_final);

%initial values;
H0 = 0;
dH_dt0 = 0;

H = [H0];
dH_dt = [dH_dt0];
tspan = t_years*(365.25/1)*(24/1)*(3600/1);%365.25 days/yr * 24 hrs/day * 3600 s /hr 
dt = 1000;%time step (s)
t = 0:dt:tspan;
n = length(t)-1;% #temp data entries = #time steps -1, to account for fact that
                % temp arrays already contain one entry (the initial
                % values)
             

%fb = 3.75; fw = linspace(-1.7/2,-1.7,n); fi = -0.8; fd = 0.75;
%f = fb+fw+fi+fd;
f = linspace(feedback_initial,feedback_final,length(t));
for i = 1:n %start @ 2 since first element is initial vals
    dH_dt(i+1) = a*(Ta(i) - Tind) + b*dTa_dt;
    H(i+1) = H(i) + dt*(a*(Ta(i) - Tind) + b*dTa_dt);
end

dH_dt = dH_dt*3600*24*365.25;%convert back to mm/yr after calculations in mm/s
t = t*(1/3600)*(1/24)*(1/365.25);%{s}*(1hr/3600s)*(1 day/24 hrs)*(1 yr/365.25 days) = years

end
