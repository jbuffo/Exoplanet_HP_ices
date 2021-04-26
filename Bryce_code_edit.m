%Pressure Profile
tic
 
clc;
% clear all;
close all;
 
%set(groot,'defaultFigureVisible','off'); %commented out plots so no need
addpath('/Users/jacobbuffo/Documents/MATLAB/SeaFreeze-master/Matlab')

rho1=917; %density of ice I
rho2=1170; %density of ice II
rho5=1230; %density of ice V
rho6=1310; %density of ice VI
 
g=23.7; %m/s^2
 
h1=13.804; %height of glacier up to Ice I and Ice II boundary
h2=21.0166; %height of glacier between Ice II and Ice V boundary
h3=41.5991; %height of glacier between Ice V and Ice VI boundary
h4=80.2502; %height of glacier between Ice VI and Ice VII boundary

%Depth (in km) Resolution Goes Below:
 
start_height=1;
resolution=1;
final_height=80;
Height_list=[];
Pressure_list=[];
Phase_list=[];
Density_list=[];
Specific_Heat_list=[];
Thermal_Conductivity_list=[];
Melting_Temperature_list=[];
 
% i=1;
Height_list=start_height:resolution:final_height;

% jb: reduced to 1 loop over a predefined height_list, then index with 'i'.
% This should speed things up. Could speed up even more if you predefine
% other lists as zero arrays of the right size and populate them
for i=1:length(Height_list)
    
h=Height_list(i);
 
%Heights = table(h)
 
%Height_list(end+1)=h; 
 
%F(h,g,dz) -> [P(z),n(z),c(z),k(z),rho(z),Tm(z)]
 
%For each depth (h) in km with a given resolution (dh) gives:
 
%Pressure (P) 
%Ice type (n) 
%Specific heat (c)
%Thermal conductivity (k) 
%Density (rho)
%Melting temperature (Tm) 
 
%P=rho*g*h
%h=P/(rho*g)
    
if (h<13.804)
    
    P1=(rho1*g*h)*1000; %pressure of glacier up to Ice I and Ice II boundary
 
elseif (h>=13.804) && (h<21.0166)
    
    P1=((rho1*g*h1) + (rho2*g*(h-h1)))*1000; %pressure of glacier between Ice II and Ice V boundary
    
elseif (h>=21.0166) && (h<41.5991)
    
    P1=((rho1*g*h1) + (rho2*g*(h2-h1))+(rho5*g*(h-h2)))*1000; %pressure of glacier between Ice V and Ice VI boundary
    
else %(h>=41.5991) && (h<80.2502);
    
    P1=((rho1*g*h1) + (rho2*g*(h2-h1))+(rho5*g*(h3-h2)+(rho1*g*(h-h3))))*1000; %pressure of glacier between Ice VI and Ice VII boundary
 
end

 
%Calculator to know ice phase
 
Pressure_list(i)=P1.*(10^-6); %give pressure in MPa
 
Pressure=Pressure_list(i);
 
T=235;   %give temperature in Kelvin-assume constant temp for now (surface temp of LHS 1140 b)

Phase= SF_WhichPhase({Pressure,T}); 

Phase_list(i)=Phase;
 
%Pressures = table(Pressure)
%Phases = table(Phase)
 
%output gives:
% 0=liquid
% 1= 'Ih' for ice Ih 
% 2= 'II' for ice II 
% 3= 'III' for ice III 
% 5= 'V' for ice V 
% 6 = 'VI' for ice VI
 
%SeaFreeze Calculator to get the specific heat and density

if Phase==0 
    PT = {Pressure,T};
    out=SeaFreeze(PT,'water1');
 
elseif Phase==1
    PT = {Pressure,T};
    out=SeaFreeze(PT,'Ih');
 
elseif Phase==2
    PT = {Pressure,T};
    out=SeaFreeze(PT,'II');
 
elseif Phase==3
    PT = {Pressure,T};
    out=SeaFreeze(PT,'III');
    
elseif Phase==5
    PT = {Pressure,T};
    out=SeaFreeze(PT,'V');
 
else %Phase==6
    PT = {Pressure,T};
    out=SeaFreeze(PT,'VI');
 
end
   

Densities=out.rho;
Specific_Heats=out.Cp;
 
%C = struct2table( out )
%Densities=C(1,6) %(kg/m^3)
 
%Specific_Heats=C(1,7) %(J/kg/K)
 

 
Density_list(i)=Densities;
 
Specific_Heat_list(i)=Specific_Heats;
 
% Cp = gives specific heat J/kg K
% rho = gives density in kg/m^3
 
%Thermal conductivity calculator(from Andersson 1 paper)

% jb: you could move these constants outside of the loop for speed
 
D1= 630;
D2= 695;
D3= 93.2;
D5= 38;
D6= 50.9;
x1= .995;
x2= 1.097;
x3= .822;
x5= .612;
x6= .612;
 
 

    
if (Phase == 0);
    Thermal_Conductivity=.60; %(W/(m*K))
 
elseif (Phase == 1) && (T>=40) && (T<=180);
    Thermal_Conductivity=D1*T^(-x1);
    
elseif (Phase == 2) && (T>=120) && (T<=240);
    Thermal_Conductivity=D2*T^(-x2);
 
elseif (Phase == 3) && (T>=180) && (T<=250);
    Thermal_Conductivity=D3*T^(-x3);
 
elseif (Phase == 5) && (T>=240) && (T<=270);
    Thermal_Conductivity=D5*T^(-x5);
 
elseif (Phase == 6) && (T>=135) && (T<=250);
    Thermal_Conductivity=D6*T^(-x6);
 
else
    Thermal="Not Available";
    Thermal_Conductivity= str2double(Thermal);
end

 
Thermal_Conductivity_list(i)=Thermal_Conductivity;
 
 
%Water Melting Points Diagram

minP=0; %MPa
stepP=10;
maxP=2300;
minT=200; %K
maxT=400;
T1 = [minT:maxT]; %temperature range
P1 = [minP:maxP]; %pressure range
out = SF_WhichPhase({minP:stepP:maxP,minT:maxT});   % careful redefining 'out' variable again, maybe 'out2'
% imagesc(minT:maxT,minP:maxP,out);
% ylabel('Pressure (MPa)')
% xlabel('Temperature (K)')
% hcb=colorbar;
% title(hcb,'Ice Phase')
% figure;
 
 
nT1=repmat(T1,231,1);
 
mT=nT1(find(out==0));
 
P1_short=minP:stepP:maxP;
 
nP1=reshape(repmat(P1_short',201,1)', 231,201);
 
mP=nP1(find(out==0));
 
% scatter(mT,mP);
% ylabel('Pressure (MPa)')
% xlabel('Temperature (K)')
% title('Melting Points for Water')
% hold on

%Get melting points of that pressure
 
%Psimple=round(P,-1)
%plot(mT,Psimple,'r*')

%plot(mT,Pressure,'r*');
 
tol = 5;
MeltingTemps=mT(abs(mP-Pressure) < tol);
Melting_Temperature=MeltingTemps(1);
 
Melting_Temperature_list(i)= Melting_Temperature;
 
%Melting_Temperatures = table(Melting_Temperature)
 
%i=i+1;
end

Pressure_Profile = cat(1,Height_list,Pressure_list, Phase_list, Density_list,Specific_Heat_list,...
    Thermal_Conductivity_list,Melting_Temperature_list);
 
 
Pressure_Profile_Final = array2table(Pressure_Profile,'RowNames',{'Height(km)','Pressure(MPa)',...
    'Phase','Density(kg/m^3)','Specific Heat((J/kg/K)','Thermal Conductivity(W/m/K)','Melting Temperature(K)'});

% Making a figure of all properties
figure
subplot(2,3,1)
plot(Pressure_list,-Height_list,'LineWidth',2)
title("Pressure vs. Depth Below Surface")
xlabel("Pressure (MPa)")
ylabel("Depth Below Surface (km)")
subplot(2,3,2)
plot(Phase_list,-Height_list,'LineWidth',2)
title("Phase vs. Depth Below Surface")
xlabel("Ice Phase")
ylabel("Depth Below Surface (km)")
xlim([0 7])
subplot(2,3,3)
plot(Density_list,-Height_list,'LineWidth',2)
title("Density vs. Depth Below Surface")
xlabel("Density (kg/m^3)")
ylabel("Depth Below Surface (km)")
subplot(2,3,4)
plot(Specific_Heat_list,-Height_list,'LineWidth',2)
title("Specific Heat vs. Depth Below Surface")
xlabel("Specific Heat (J/(kg K))")
ylabel("Depth Below Surface (km)")
subplot(2,3,5)
plot(Thermal_Conductivity_list,-Height_list,'LineWidth',2)
title("Thermal Conductivity vs. Depth Below Surface")
xlabel("Thermal Conductivity (W/(m K))")
ylabel("Depth Below Surface (km)")
subplot(2,3,6)
plot(Melting_Temperature_list,-Height_list,'LineWidth',2)
title("Melting Temperature vs. Depth Below Surface")
xlabel("Melting Temperature (K)")
ylabel("Depth Below Surface (km)")

toc