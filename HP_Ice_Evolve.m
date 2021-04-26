function T_np1=HP_Ice_Evolve(T_n,k,rho,c,dt,dz,T_surf,Base_flux)

% Function that takes the original vertical profiles of ice temperature (T),
% thermal conductivity (k), density (rho), and heat capacity (c), alongside
% time step (dt), vertical resolution (dz - in meters), surface temperature
% (in Kelvin), and basal heat flux (in W/m^2) and calculates the thermal 
% evolution of the ice.


% Creating shifted thermal conductivity vectors and up/down conductivity
% averages
k_jp1=circshift(k,1);
k_jp1(1)=k(1);
k_jm1=circshift(k,-1);
k_jm1(end)=k(end);

k_up_avg=(k_jp1+k)./2;
k_down_avg=(k_jm1+k)./2;

% Constant coefficient for solution matrix
C_0=dt./((dz^2)*rho.*c);

% Build diagonal vectors of solution matrix to put in Thomas Tridiagonal
% algorithm
a=zeros(1,length(T_n));
b=zeros(1,length(T_n));
c=zeros(1,length(T_n));

for i=1:length(T_n)
    a(i)=1+C_0(i)*(k_up_avg(i)+k_down_avg(i));
    b(i)=-C_0(i)*k_down_avg(i);
    c(i)=-C_0(i)*k_up_avg(i);
end

b(end)=[];
c(1)=[];

% Temp profile from previous time step (used in Thom_Trid algorithm) and
% inclusion of temperature boundary conditions
y=zeros(1,length(T_n));
y(1:end)=T_n;

% surface BC (surface temp)
y(1)=y(1)+C_0(1)*k_up_avg(1)*T_surf;
% basal BC (basal heat flux)
y(end)=y(end)+C_0(end)*k_down_avg(end)*(T_n(end)+Base_flux*dz/k_down_avg(end));

% Solve for new Temperature profile
T_np1=Thomas_Trid(a,b,c,y);





