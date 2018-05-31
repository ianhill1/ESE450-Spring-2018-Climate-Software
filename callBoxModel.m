%% calls BoxModel
clear all
close all

%% Model inputs
%give initial and final feedback parameters over observed time scale
alpha_init = -2.37;%CMIP5 mean net feedback 0-20 years {W * m^-2 * K^-1}
alpha_final = -1.33;%CMIP5 mean net feedbac k 21-150 years {W * m^-2 * K^-1}
alphaJumpTime = 21;%time when alpha changes (year) according to model

%choose a method for varying alpha between two experimental apha values
%note: if you do not want to change alpha, just set alpha_final=alpha_initial
%the following options my be used for varying alpha with respect to the chosen time interval:
%alphaType = 1 for linear change
%alphaType = 0 for step @ t = alphaJumpTime
%alphaType = 2 for modified moving average step change
alphaType = 2;


%ocean and atmosphere depth (used to calculate respective specific heats)
OceanDepth = 1000;% m
AtmDepth = 10000;% m

% time range
year_begin = 2010;
year_end = 2110;


%% call models
years = year_end - year_begin;
[Ta,To,t,dTa_dt,alpha] = BoxModel(years,alpha_init,alpha_final,AtmDepth,OceanDepth,alphaType,alphaJumpTime);
[dH_dt,H,t1] = OceanRise(years,alpha_init,alpha_final,AtmDepth,OceanDepth,alphaType,alphaJumpTime);


%% plots
alphaInitStr = num2str(alpha_init);
alphaFinalStr = num2str(alpha_final);
lw = 2;%line width
xTix = year_begin+10:10:year_end;%axis labels (+10 makes first tick @year + 10...)
x_labels = strsplit(num2str(xTix));


figure(1)
plot(t,Ta,'r.')
hold on
plot(t,To,'k--','Linewidth',lw)
xlim([t(1500)*1.1,t(end)])%gets rid of initial spike from zero
xlabel('Year','fontsize',15)
ylabel('Temperature (°C)','fontsize',15)
xticklabels(x_labels)
legend('Atm temp', 'Ocean temp','Location','southeast')
titleStr = strcat('Temperature Rise From ' ,'  ', num2str(year_begin));
title(titleStr,'fontsize',15)
grid on

figure(2)
subplot(2,1,1)
plot(t1,dH_dt,'Linewidth',lw)%plot time in years and ocean rise rate in mm/yr
xlabel('Year','fontsize',15)
ylabel('Ocean Rise Rate (mm/yr)','fontsize',15)
titleStr = strcat('Sea Level Rise/ Rise Rate From ' ,'  ', num2str(year_begin));
title(titleStr,'fontsize',15)
ylim([0 , dH_dt(end)*1.02])
xTix = year_begin:10:year_end;%axis labels (+10 makes first tick @year + 10...)
x_labels = strsplit(num2str(xTix));
xticklabels(x_labels)
grid on
subplot(2,1,2)
plot(t1,H,'Linewidth',lw)
xlabel('Year','fontsize',15)
ylabel('Ocean Rise (mm)','fontsize',15)
xticklabels(x_labels)
grid on

figure(3)
plot(t,alpha,'m.')
title('Feedback Parameter Vs. Time','fontsize',15)
xlabel('Year','fontsize',15)
ylabel('Feedback {\alpha} (W/m^2-K)','fontsize',15)
xticklabels(x_labels)
grid on







