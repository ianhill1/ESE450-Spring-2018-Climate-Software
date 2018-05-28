%% calls BoxModel
clear all;
close all;

%give initial and final feedback parameters over observed time scale
f_init = 2.03;
f_final = 1.92;

%give time range
year_begin = 2010;
year_end = 2100;
years = year_end - year_begin;

%call models
[Ta,To,t,dTa_dt] = BoxModel(years,f_init,f_final);
[dH_dt,H,t] = OceanRise(years,f_init,f_final);


%% plots
lw = 1;%line width
xTix = year_begin:10:year_end;%axis labels
x_labels = strsplit(num2str(xTix));


figure(1)
plot(t,Ta,'r--')
hold on
plot(t,To,'k.')
xlabel('Year')
ylabel('Temperature (°C)')
xticklabels(x_labels)
grid on

figure(2)
subplot(2,1,1)
plot(t,dH_dt,'Linewidth',lw)%plot time in years and ocean rise rate in mm/yr
xlabel('Year')
ylabel('Ocean Rise Rate (mm/yr)')
ylim([dH_dt(2) , dH_dt(end)*1.02])
xticklabels(x_labels)
grid on
subplot(2,1,2)
plot(t,H,'Linewidth',lw)
xlabel('Year')
ylabel('Ocean Rise (mm)')
xticklabels(x_labels)
grid on


