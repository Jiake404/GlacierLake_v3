function [model_stage,enthalpy,lambda,T,grid_profile,k,c,SW_prop,Io,tau,q,lake_depth,...
    surface_lake_depth,sub_lake_depth,lake_albedo,lid_thick,lake_T_av,lf_upper,lf_lower,...
    d_enthalpy,sur_lake_T_av,snow_z,snow_mwe,snow_lambda,snow_ice,snow_k,snow_c,snow_T,...
    snow_e,snow_ro,snow_l_ro,track_a,track_b] = creat_array(total_grid_num,tau_ice,...
    ice_grid_num,ice_grid_z,deep_ice_grid_z)
model_stage=3;
enthalpy = zeros(total_grid_num, 1);        %(J) enthalpy of cells
lambda = zeros(total_grid_num, 1);          %water content (0 = no water, 1 = all water)
T = zeros(total_grid_num, 1);               %(K) temperature
grid_profile = zeros(total_grid_num, 1);
grid_profile(1:ice_grid_num,1) = ice_grid_z; 
grid_profile(ice_grid_num + 1:total_grid_num,1) = deep_ice_grid_z; 
k = zeros(total_grid_num, 1);               %(W/(m.K)) thermal conductivity
c = zeros(total_grid_num, 1);               %(J/K) heat capacity
SW_prop = zeros(total_grid_num, 1);         %SW propagation through ice
Io = zeros(total_grid_num, 1);              %(/m) fraction of shortwave that can propagate beyond surface
tau = ones(total_grid_num, 1)*tau_ice;      %(/m) bulk shortwave extinction
q = zeros(1, 1);                            %(J) energy transfer at the surface
lake_depth = zeros(1, 1);                   %(m) depth of main lake
surface_lake_depth = zeros(1, 1);           %(m) depth of lake that forms on lid surface
sub_lake_depth = zeros(1, 1);               %(m) depth of subsurface lake (beneath lid)
lake_albedo = zeros(1, 1);                  %empty array to record lake albedo change
lid_thick = zeros(1, 1);                    %(m) thickness of ice lid
lake_T_av = zeros(1, 1);                    %(K) average (core) lake temperature
lf_upper = zeros(1, 1);                     %(J) lake flux upper per time step
lf_lower = zeros(1, 1);                     %(J) lake flux lower per time step
d_enthalpy = zeros(1, 1);                   %(J) total flux leaving lake core
sur_lake_T_av = zeros(1, 1);                %(K) average (core) surface lake temperature
snow_z = zeros(1, 1);                       %(m) overall snow depth including refrozen water
snow_mwe = zeros(1, 1);                     %(mwe) depth of each layer in mwe
snow_lambda = zeros(1, 1);                  %lambda of each layer of the snow pack model (0 = no water, 1 = all water)
snow_ice = zeros(1, 1);                     %of snow pack layer that is not water, how much is snow and how much is snow-ice (0 = all snow, 1 = all snow-ice)
snow_k = zeros(1, 1);                       %(W/(m.K)) snow thermal conductivity of each snow layer
snow_c = zeros(1, 1);                       %(J/K) heat capacity of each snow layer
snow_T = zeros(1, 1);                       %(K) snow temperature
snow_e = zeros(1, 1);                       %(J) snow enthalpy
snow_ro = zeros(1, 1);                      %(kg/m^3)density of each snow layer
snow_l_ro = zeros(1, 1);                    %(kg/m^3)overall density of snow layer, incl water and snow ice
track_a = zeros(1, 1);                      %tracking array, keep everything in this one and add rows as required.
track_b = zeros(4,1);                           %second tracking array, this one does not record across whole 1