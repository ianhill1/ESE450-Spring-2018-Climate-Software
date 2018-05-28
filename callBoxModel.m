%% calls BoxModel
clear all
close all

%% Model inputs
%give initial and final feedback parameters over observed time scale
alpha_init = -1.25;%CMIP5 mean net feedback 0-20 years {W * m^-2 * K^-1}
alpha_final = -1.33;%CMIP5 mean net feedback 21-150 years {W * m^-2 * K^-1}

%ocean and atmosphere depth (used to calculate respective specific heats)
OceanDepth = 1000;% m
AtmDepth = 10000;% m

% time range
year_begin = 2010;
year_end = 2110;


%% call models
years = year_end - year_begin;
[Ta,To,t,dTa_dt] = BoxModel(years,alpha_init,alpha_final,AtmDepth,OceanDepth);
[dH_dt,H,t1] = OceanRise(years,alpha_init,alpha_final,AtmDepth,OceanDepth);


%% plots
alphaInitStr = num2str(alpha_init)
alphaFinalStr = num2str(alpha_final)
lw = 1;%line width
xTix = year_begin+10:10:year_end;%axis labels (+10 makes first tick @year + 10...)
x_labels = strsplit(num2str(xTix));


figure(1)
plot(t,Ta,'r.')
hold on
plot(t,To,'k--','Linewidth',lw)
xlim([t(1500)*1.1,t(end)])%gets rid of initial spike from zero

xlabel('Year')
ylabel('Temperature (°C)')
xticklabels(x_labels)
legend('Atm temp', 'Ocean temp')
titleStr = strcat('Temperature Rise From ' ,'  ', num2str(year_begin))
titleStr2 = strcat('({\alpha} varies linearly from   ',num2str(alpha_init),' to    ' ,num2str(alpha_final),' W/m^2-K  )')
title({titleStr; titleStr2},'fontsize',15)
grid on

figure(2)
subplot(2,1,1)
plot(t1,dH_dt,'Linewidth',lw)%plot time in years and ocean rise rate in mm/yr
xlim([t(1500)*1.1,t(end)])%gets rid of initial spike from zero
xlabel('Year')
ylabel('Ocean Rise Rate (mm/yr)')
titleStr = strcat('Sea Level Rise/ Rise Rate From ' ,'  ', num2str(year_begin))
%titleStr2 = strcat('({\alpha} varies linearly from   ',num2str(alpha_init),' to    ' ,num2str(alpha_final),' W/m^2-K  )')
%title({titleStr; titleStr2},'fontsize',15)
title(titleStr,'fontsize',15)
%ylim([dH_dt(2) , dH_dt(end)*1.02])
xticklabels(x_labels)
grid on
subplot(2,1,2)
plot(t1,H,'Linewidth',lw)
xlabel('Year')
ylabel('Ocean Rise (mm)')
xticklabels(x_labels)
xlim([t(1500)*1.1,t(end)])%gets rid of initial spike from zero

grid on


