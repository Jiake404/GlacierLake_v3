function [T,enthalpy,lambda,grid_profile,total_grid_num,c,Io,k,SW_prop,tau] = add_water(ex_water,grid_profile_temp,...
    T_temp,enthalpy_temp,lambda_temp,tt,hydro_T,ro_water,ice_grid_z,c_ice,...
    T_melt,Lf,c_water,c_temp,Io_temp,k_temp,SW_prop_temp,tau_temp)

T1 = zeros(ex_water,1);
enthalpy1 = zeros(ex_water,1);
lambda1 = zeros(ex_water,1);
grid_profile1 = zeros(ex_water,1);
c1 = zeros(ex_water,1);
Io1 = zeros(ex_water,1);
k1 = zeros(ex_water,1);
SW_prop1 = zeros(ex_water,1);
tau1 = zeros(ex_water,1);

T1(:,tt) = hydro_T;
enthalpy1(:,tt) = ro_water.*ice_grid_z.*...
        (c_ice*T_melt + Lf + c_water.*hydro_T - c_water.*T_melt);
lambda1(:,tt) = 1;
grid_profile1(:,tt) = ice_grid_z;

T = [T1;T_temp];
enthalpy = [enthalpy1;enthalpy_temp];
lambda = [lambda1;lambda_temp];
grid_profile = [grid_profile1;grid_profile_temp];
c = [c1;c_temp];
Io = [Io1;Io_temp];
k = [k1;k_temp];
SW_prop = [SW_prop1;SW_prop_temp];
tau = [tau1;tau_temp];

total_grid_num = length(grid_profile);

