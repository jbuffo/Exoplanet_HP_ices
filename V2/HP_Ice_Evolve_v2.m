function [T_np1,Phi_np1]=HP_Ice_Evolve_v2(T_n,Phi_n,k_i,rho_i,c_i,k_w,rho_w,...
    c_w,dt,dz,T_surf,Base_flux,Tm,L,TTol,PhiTol)

% Improved function that takes the original vertical profiles of ice
% temperature (T_n), porosity (Phi_n), thermal conductivity (k_i), density
% (rho_i), heat capacity (c_i), and melting temperature (Tm) alongside time
% step (dt), vertical resolution (dz - in meters), surface temperature (in
% Kelvin), basal heat flux (in W/m^2), latent heat of fusion (L), water
% properties (k_w, rho_w, c_w), and error tolerances for the finite
% difference iterator (TTol and PhiTol) and calculates the physical and
% thermal evolution of the ice.

% This function does thermal diffusion and phase change (via the enthalpy
% method)

matrix_dimension=length(T_n);
counter=0;

T_np1_km1=T_n;
phi_np1_km2=Phi_n;
T_evolve=[];
phi_evolve=[];
TErr=999;
PhiErr=999;

while  TErr>TTol || PhiErr>PhiTol
    counter=counter+1;

%% Enthalpy Method
Hs=c_i.*Tm;     % enthalpy of solid ice to compar in enthalpy method

% Value of Phi(n+1,k-1,j)
phi_np1_km1=phi_np1_km2;
    for i=1:matrix_dimension
% Test to see if melting ice or freezing brine
    En=c_i(i)*T_np1_km1(i)+L*phi_np1_km2(i);
        if En<Hs(i)
            phi_np1_km1(i)=0;
        elseif En>Hs(i)+L
            phi_np1_km1(i)=1;
        else
            phi_np1_km1(i)=((c_i(i)*T_np1_km1(i)+phi_np1_km2(i)*L)-Hs(i))/...
            L;
        end
    end

% Reassigning Phi(n+1,k-2) as Phi(n+1,k-1) for next
% iteration
phi_np1_km2=phi_np1_km1;


%% Volume averaging conductivity, density*specific heat
k_bar=phi_np1_km1.*k_w+(1-phi_np1_km1).*k_i;
rho_c_bar=phi_np1_km1.*rho_w.*c_w+(1-phi_np1_km1).*rho_i.*c_i;

% Creating shifted thermal conductivity vectors and up/down conductivity
% averages
k_jp1=circshift(k_bar,1);
k_jp1(1)=k_bar(1);
k_jm1=circshift(k_bar,-1);
k_jm1(end)=k_bar(end);

k_up_avg=(k_jp1+k_bar)./2;
k_down_avg=(k_jm1+k_bar)./2;

% Constant coefficient for solution matrix
C_0=dt./((dz^2)*rho_c_bar);

% Build diagonal vectors of solution matrix to put in Thomas Tridiagonal
% algorithm
a=zeros(1,length(T_np1_km1));
b=zeros(1,length(T_np1_km1));
c=zeros(1,length(T_np1_km1));

for i=1:length(T_np1_km1)
    a(i)=1+C_0(i)*(k_up_avg(i)+k_down_avg(i));
    b(i)=-C_0(i)*k_down_avg(i);
    c(i)=-C_0(i)*k_up_avg(i);
end

b(end)=[];
c(1)=[];

% Temp profile from previous time step (used in Thom_Trid algorithm) and
% inclusion of temperature boundary conditions
y=zeros(1,length(T_np1_km1));
y(1:end)=T_n-((L*rho_i./rho_c_bar).*(phi_np1_km1-Phi_n));

% surface BC (surface temp)
y(1)=y(1)+C_0(1)*k_up_avg(1)*T_surf;
% basal BC (basal heat flux)
y(end)=y(end)+C_0(end)*k_down_avg(end)*(T_np1_km1(end)+Base_flux*dz/k_down_avg(end));

% Solve for new Temperature profile
T_np1_km1=Thomas_Trid(a,b,c,y)';

    %% Appending value to matrix to check for convergence
    T_evolve=[T_evolve T_np1_km1'];
    phi_evolve=[phi_evolve phi_np1_km1'];
    
    if counter==1
        TErr=999;
        PhiErr=999;
    else
        TErr=max(abs(T_evolve(:,counter)-T_evolve(:,counter-1)));
        PhiErr=max(abs(phi_evolve(:,counter)-phi_evolve(:,counter-1)));
    end
end

%% OUTSIDE OF ITERATIVE LOOP - advective mixing of any present water pockets
top=0;
size=0;
T_hold=0;
depth=1;

while depth <= matrix_dimension
    if phi_np1_km1(depth) == 0
        depth=depth+1;
        top=depth;
    elseif phi_np1_km1(depth) > 0
        test=phi_np1_km1(top);
        while test > 0
            T_hold=T_hold+T_np1_km1(top+size);
            size=size+1;
            if top+size > matrix_dimension
                test=0;
            else
                test=phi_np1_km1(top+size);
            end
        end
        T_np1_km1(top:top+size-1)=T_hold/size;
        depth=top+size;
        size=0;
        T_hold=0;
        top=depth;
    end
end
    


%% Final temperature and porosity profiles
T_np1=T_np1_km1';
Phi_np1=phi_np1_km1';

