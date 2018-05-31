

function [Ta,To,t,dTa_dt,alpha,dt] = BoxModel(t_years,alpha_initial,alpha_final,atmDepth,oceanDepth,alphaType,alphaJumpTime)
Ca = 750*atmDepth;%J/m^2-K
Co = 3.35e+6*oceanDepth;
cao = 14.2;
Q = 4.0;%W/m^2

%unit convertions
secondsPerYear = (365.25/1)*(24/1)*(3600/1);%{years}*(365.25 days/yr) * (24 hrs/day) * (3600 s /hr) = {seconds} 
yearsPerSecond = (1/3600)*(1/24)*(1/365.25);%{s}*(1hr/3600s)*(1 day/24 hrs)*(1 yr/365.25 days) = {years}
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
tspan = t_years*secondsPerYear;%convert time to seconds 
dt = 1000;%time step (s)
t = 0:dt:tspan;
n = length(t)-1;% #temp data entries = #time steps -1, to account for fact that
                % temp arrays already contain one entry (the initial
                % values)
%alpha = linspace(alpha_initial,alpha_final,length(t));                
%% define alpha
alpha = zeros(1,length(t));
if alphaType == 1
    alpha = linspace(alpha_initial,alpha_final,length(t));     
end

if alphaType == 0
    for i  = 1:length(alpha)
        if t(i)*yearsPerSecond < alphaJumpTime%if time (years) < 21, use first value of alpha
            alpha(i) = alpha_initial;
        else %t > 21 years
            alpha(i) = alpha_final;
        end %endif
    end %endfor
end %endif

%define number of data entries per year for moving average
GridsPerYear = n/t_years;
NumberYearsMovingAverage = 9;
if alphaType == 2
    for i  = 1:length(alpha)
        if t(i)*yearsPerSecond < alphaJumpTime%if time (years) < 21, use first value of alpha
            alpha(i) = alpha_initial;
        else %t > 21 years
            alpha(i) = alpha_final;
        end %endif
    end %endfor
    alpha = movmean(alpha,NumberYearsMovingAverage*GridsPerYear);
end %endif
    
   
%%
for i = 1:n %start @ 2 since first element is initial vals
    %Ta(i+1) = Ta(i) + dt*( -(-f(i)/Ca + cao/Ca)*Ta(i) + (cao/Ca)*To(i) + Q/Ca);
    Ta(i+1) = Ta(i) + dt*( -(-alpha(i)/Ca + cao/Ca)*Ta(i) + (cao/Ca)*To(i) + Q/Ca);
    To(i+1) = To(i) + dt*( (cao/Co)*Ta(i) - (cao/Co)*To(i) );
    dTa_dt(i) = -(-alpha(i)/Ca + cao/Ca)*Ta(i) + (cao/Ca)*To(i) + Q/Ca;
end
t = t*yearsPerSecond;%{s}*YearsPerSecond = years
end


