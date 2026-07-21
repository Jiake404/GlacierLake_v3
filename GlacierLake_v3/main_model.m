clear;clc
tic

%USER INPUTS
AWS = 'L3';            % Station name L1,L2
run_start = [2021,01,01];     % Model initiation time L1:2021; L2:2020
run_end = [2021,12,31];       % Model termination time
gcnet_num = 0;                % Weather station from 0: PROMICE or 1: GC-NET

%add directories
working_directory = cd;
addpath(strcat(working_directory,'\Data'))

% [clip_hour1,clip_hour2,clip_day1,clip_day2,run_day,run_hour] = readtime(AWS,run_start,run_end);

%===================ERA5 input=====================
run_day = day(datetime(run_end),'dayofyear');
demname = 'L32dem.tif'; % L1:L12dem.tif; L2:L22dem.tif
filename = 'lake_depth_L3.xlsx';  
sheetname = '2021'; % L1:2021; L2:2020
data_zhou = readtable(filename, 'Sheet', sheetname);
area_zhou = data_zhou{:,2};
vols_zhou = data_zhou{:,3};
areaz = area_zhou(2:end,:)./1000000;
volsz = vols_zhou(2:end,:)./1000000000;
doy_zhou = data_zhou{:,4};
doyz = doy_zhou(2:end,:);
%==================================================

% Parameter setting
[t_step, model_stage, albedo_mult, precip_mult, albedo_add, ice_grid_num, ice_grid_z,...
    deep_ice_grid_num, deep_ice_grid_z, lower_boundary, low_plot_T, hi_plot_T, ice_lim,...
    snow_threshold, slush_lid_num, new_snow_albedo, snow_albedo_min, max_hydro_input,...
    breakup_threshold, lake_turb_threshold, lake_prof_threshold, slush_lid_threshold,...
    slush_mech_threshold, stage_print_hysteresis, refresh_albedo, tau_snow_cold, tau_snow_melt,...
    input_figs, output_figs, hydro_input, sp_year, data_import, import_save, AWS_albedo, basal_SW_distribute,...
    ice_albedo, Io_ice, Io_water, Io_slush, tau_ice, tau_water, tau_slush, tau_snow, J, ro_snow_initial,...
    b_exp, ro_melt_max, ro_cold_max, tau_ro, k_snow_ice_max, T_melt, Lf, ro_water, c_ice, c_water,...
    k_water, k_ice, k_air, e_ice, T_ro_max, ro_ice, ies80, t_total, t_num, total_grid_num, grid_profile,...
    lake_grid_num,hydro_T,ice_T_bottom,ice_T_surface] = cons_output(run_day);

if strcmp(AWS, 'L1')
    % import AWS data
    [SW_down_out, LW_down_out, air_T_out, rh_out, hum_out, pressure_out,...
        wind_speed_out, albedo_out, s_data, time_out] = import_AWS1(AWS,clip_hour1,clip_hour2,run_start,gcnet_num);

    % import precipitation data
    [precip_out,precip_s] = import_snow1(AWS,clip_day1,clip_day2,run_start,run_end,t_num);
else
    load('ERA5_L3_2021.mat')
    s_data = [SW_down_out,albedo_out,LW_down_out,air_T_out,rh_out,hum_out,pressure_out,wind_speed_out];
    load('precip_L3_2021.mat')
    precip = precip_L3_2021;
    load('snowfall_L3_2021.mat')
    sf = snowfall_L3_2021;
    precip_out=[precip';sf';precip'-sf'];
    precip_out(precip_out<0)=0;
    precip_out=precip_out./1000;
    original_length = length(precip_out);
    precip_out = imresize(precip_out, [3 run_day*24], 'bilinear'); %reshape to the correct size
    precip_out = precip_out*(original_length/(run_day*24)); %correct for stretching
    precip_s = precip_out;
end

% spin_up
[model_stage, ini_condtions] = spin_up_era5(s_data, precip_s, t_step, ice_T_bottom, ice_T_surface,...
    albedo_mult, model_stage, ice_grid_num, ice_grid_z, deep_ice_grid_num, deep_ice_grid_z,...
    lower_boundary, snow_threshold, new_snow_albedo, snow_albedo_min, refresh_albedo,...
    tau_snow_cold, tau_snow_melt, AWS_albedo, ice_albedo, Io_ice, Io_water, Io_slush,...
    ro_snow_initial, b_exp, ro_melt_max, ro_cold_max, tau_ro, k_snow_ice_max,...
    T_melt, Lf, ro_water, c_ice, c_water, k_water, k_ice, k_air, e_ice, ro_ice, sp_year, gcnet_num);

% [hydro_out] = import_runoff(AWS,clip_day1,clip_day2,run_start,run_end,t_num);
[hydro_out,melt_ice,~,Q_mice,~,T_s,T_b,tr_a,s_rol,s_ro,s_mwe] = runoff_prod_era5(SW_down_out, LW_down_out, air_T_out, rh_out, hum_out,...
    pressure_out, wind_speed_out, albedo_out, precip_out, run_day, t_step, albedo_mult, model_stage,...
    ini_condtions, ice_grid_num, ice_grid_z, deep_ice_grid_num, deep_ice_grid_z, lower_boundary,...
    snow_threshold, new_snow_albedo, snow_albedo_min, refresh_albedo, tau_snow_cold, tau_snow_melt,...
    AWS_albedo, ice_albedo, Io_ice, Io_water, Io_slush, ro_snow_initial, b_exp, ro_melt_max,...
    ro_cold_max, tau_ro, k_snow_ice_max, T_melt, Lf, ro_water, c_ice, c_water, k_water, k_ice,...
    k_air, e_ice, ro_ice, gcnet_num, AWS);
s_tnum = find(hydro_out>0,1,'first');

%==========================dem input===============================
[dem,R] = geotiffread(demname);
dem(dem<0)=nan;
[dem]=filling(dem);
% Create a binary mask
BW = ~isnan(dem);
% find dem's boundary
B = bwperim(BW); 
% find index of boundary
[row, col] = find(B); 
% dem of boundary
edgeValues = dem(B);  
% pour point
[minValue, minIndex] = min(edgeValues, [], 'omitnan');
spilloverRow = row(minIndex);
spilloverCol = col(minIndex);
% dem resolution
cell_size = 100;
dem_size = size(dem);

%============================runoff_loss_rate==============================
cc=0.67;%L1:0.5; L2:0.8
lake_vols=[];lake_areas=[];
for cci = 1:length(cc)
    state_flag = zeros([dem_size,t_num-1]);
    state_flag(dem<0) = nan;
    s_flag = zeros(dem_size);
    external_water = zeros([dem_size,t_num-1]);
    dem_pro = zeros([dem_size,t_num-1]);
    dem_pro(:,:,1)=dem;
    lake_f = zeros(dem_size);% Record whether this grid is the first time to run the lake module

    % creating arrays
    Enthalpy = cell(dem_size);
    T_total = cell(dem_size);
    Grid_profile = cell(dem_size);
    Track_b = cell(dem_size);
    Lambda = cell(dem_size);
    Model_stage = zeros([dem_size,t_num-1]);
    Track_a = zeros([dem_size,t_num-1]);
    Snow_T = zeros([dem_size,t_num-1]);
    Snow_e = zeros([dem_size,t_num-1]);
    Snow_lambda = zeros([dem_size,t_num-1]);
    Snow_mwe = zeros([dem_size,t_num-1]);
    Snow_ro = zeros([dem_size,t_num-1]);
    Snow_l_ro = zeros([dem_size,t_num-1]);
    Snow_z = zeros([dem_size,t_num-1]);
    Snow_ice = zeros([dem_size,t_num-1]);
    Lake_depth = zeros([dem_size,t_num-1]);
    Lake_albedo = zeros([dem_size,t_num-1]);
    Lake_ind = zeros([dem_size,t_num-1]);
    bm = zeros([dem_size,t_num-1]);
    lid_depth = zeros([dem_size,t_num-1]);
    Total_grid_num = zeros(dem_size);
    Total_grid_num(:) = ice_grid_num+deep_ice_grid_num;
    New_r = repmat((1:dem_size(1))' * ones(1, dem_size(2)), [1, 1, t_num-1]);
    New_c = repmat(ones(dem_size(1), 1) * (1:dem_size(2)), [1, 1, t_num-1]);
    Remaining_dist_matrix = zeros([dem_size,t_num-1]);
    fill_flag = zeros(dem_size);


    cc_one = cc(cci);
    hydro_out_new = hydro_out.*(1-cc_one);
    s_tnum=(data_zhou{1,4}-1)*24+14;
    for tt = s_tnum:t_num

        if tt == s_tnum
            water_depth = vols_zhou(1,1)/10000;
            [filled_dem,fill_flag,external_water,s_flag,total_external_input]=flow_speed0(dem,tt,dem,external_water,s_flag,water_depth,fill_flag,spilloverRow,spilloverCol);
        else
            meltseason = 1;
            if meltseason == 1%&&~any(Model_stage(:,:,tt) == 4, 'all')

                [filled_dem,fill_flag,New_r,New_c,Remaining_dist_matrix,external_water,s_flag,total_external_input]=flow_speed2(dem_pro(:,:,tt-1),...
                    s_rol(tt),ro_water,s_mwe(tt),t_step,hydro_out_new,tt,...
                    s_tnum,New_r,New_c,Remaining_dist_matrix,fill_flag,dem,external_water,s_flag,tr_a(tt),cell_size,spilloverRow,spilloverCol);%第二个demdem_pro(:,:,s_tnum)
            end
        end
        dem_pro(:,:,tt) = filled_dem;
        ex_water = external_water(:,:,tt);
        state_flag(:,:,tt) = s_flag;
        % submerged grids
        [idx] = find(s_flag~=0&~isnan(s_flag)); 
        [rown, coln] = ind2sub(dem_size, idx);

        for i = 1:length(rown)
            m = rown(i);
            n = coln(i);

            if lake_f(m,n)==0
                % creat array
                total_grid_num = Total_grid_num(m,n);
                [model_stage,enthalpy,lambda,T,grid_profile,k,c,SW_prop,Io,tau,q,lake_depth,...
                    surface_lake_depth,sub_lake_depth,lake_albedo,lid_thick,lake_T_av,...
                    lf_upper,lf_lower,d_enthalpy,sur_lake_T_av,snow_z,snow_mwe,...
                    snow_lambda,snow_ice,snow_k,snow_c,snow_T,snow_e,snow_ro,snow_l_ro,...
                    track_a,track_b] = creat_array(total_grid_num,tau_ice,ice_grid_num,ice_grid_z,deep_ice_grid_z);
                % input initial condtions
                ice_T_surface = T_s(tt-1);
                ice_T_bottom = T_b(tt-1);
                %initialise arrays
                T(:,1) = linspace(ice_T_surface, ice_T_bottom, total_grid_num)';
                enthalpy(:,1) = c_ice.*ro_water.*T(:,1).*grid_profile(:,1);    
                lambda(:,1) = 0;                                                   
                track_a(1,1) = model_stage;
                lake_ind = 1:ex_water(m,n);
                lake_f(m,n) = 1;
            else
                model_stage = Model_stage(m,n,tt-1);
                total_grid_num = Total_grid_num(m,n);
                snow_ro_e = Snow_ro(m,n,tt-1);
                snow_mwe_e = Snow_mwe(m,n,tt-1);
                snow_ice_e = Snow_ice(m,n,tt-1);
                snow_lambda_e = Snow_lambda(m,n,tt-1);
                enthalpy = Enthalpy{m,n}(:,tt-1);
                T = T_total{m,n}(:,tt-1);
                grid_profile = Grid_profile{m,n};
                lambda = Lambda{m,n}(:,tt-1);
                track_a = Track_a(m,n,tt-1);
                track_b = Track_b{m,n};
                snow_T = Snow_T(m,n,tt-1);
                snow_e = Snow_e(m,n,tt-1);
                snow_lambda = Snow_lambda(m,n,tt-1);
                snow_mwe = Snow_mwe(m,n,tt-1);
                snow_ro = Snow_ro(m,n,tt-1);
                snow_l_ro = Snow_l_ro(m,n,tt-1);
                snow_z = Snow_z(m,n,tt-1);
                snow_ice = Snow_ice(m,n,tt-1);
                lake_depth = Lake_depth(m,n,tt-1);
                lake_albedo = Lake_albedo(m,n,tt-1);
                c = zeros(total_grid_num,1);
                Io = zeros(total_grid_num,1);
                k = zeros(total_grid_num,1);
                SW_prop = zeros(total_grid_num,1);
                tau = zeros(total_grid_num,1);
                lake_ind = 1:Lake_ind(m,n,tt-1);
            end

            if model_stage == 1||3||5
                [T,enthalpy,lambda,grid_profile,total_grid_num,c,Io,k,SW_prop,tau]...
                    = add_water(ex_water(m,n),grid_profile,T,enthalpy,lambda,1,hydro_T,...
                    ro_water,ice_grid_z,c_ice,T_melt,Lf,c_water,c,Io,k,SW_prop,tau);
            end

            snow_mwe = snow_mwe+precip_out(2,tt);

            if model_stage == 3 || model_stage == 5
                snow_mwe = 0;
            elseif model_stage == 1
                if snow_mwe>= snow_threshold
                    model_stage = 2;
                    fprintf('Moving from stage %1.0f to %1.0f. Day %3.0f. \n \n', track_a(1), model_stage, tt*t_step/24);
                end
            end

            %%==============================stage1=============================
            if model_stage == 1

                if gcnet_num == 1 %then calculate LW_in seperately
                    [LW_down_out(tt)] = longwave_in(rh_out(tt), air_T_out(tt), lambda(1), model_stage);
                end

                %surface flux
                [q] = surface_flux1(SW_down_out(tt), albedo_out(tt), air_T_out(tt),...
                    wind_speed_out(tt), pressure_out(tt), LW_down_out(tt), hum_out(tt),...
                    T(1), ice_albedo, AWS_albedo, Io_ice, e_ice,0);

                %apply incoming energy flux to enthalpy array
                enthalpy(1) = enthalpy(1) + q*3600*t_step; %*3600*t_step for amount in that timestep

                %calcualte temperature and lambda for cells based on enthalpy
                [T(:,1), lambda(:,1), model_stage] = temp_lambda_profile1(enthalpy(:,1),...
                    total_grid_num, grid_profile(:,1), T_melt, c_ice, Lf, model_stage,...
                    ro_water, c_water);

                %only continue if model stage still = 1
                if model_stage == 1

                    %update thermal conductivity and heat capacity to that of ice
                    k(:,1) = k_ice;
                    c(:,1) = c_ice;

                    %conduction calculations
                    [enthalpy(:,1)] = conduct_update_enthalpy_new(k(:,1), c(:,1), T(:,1),...
                        grid_profile(:,1), enthalpy(:,1), t_step, total_grid_num, ro_water, lower_boundary,...
                        model_stage);
                    track_a = 1;
                else
                    track_a(1) = 1;
                end
            end

            %%==============================stage2=============================
            if  model_stage == 2

                if gcnet_num == 1 %then calculate LW_in seperately
                    [LW_down_out(tt)] = longwave_in(rh_out(tt), air_T_out(tt), lambda(1,tt), model_stage);
                end

                if track_a ~= 2 %if coming from another stage then initialise

                    %initialise snow variables
                    [snow_lambda, snow_ro, snow_c, snow_T,snow_e, snow_k, snow_z,...
                        snow_ice, snow_l_ro] = initialise_snow1(snow_mwe, ro_snow_initial,...
                        c_ice, k_ice, k_air, ro_water, T_melt, air_T_out(tt));

                    snow_albedo = new_snow_albedo*albedo_mult; if snow_albedo > 1; snow_albedo = 1; end

                elseif track_a == 2 %if within stage 2 then update as normal

                    if tt == 1
                        [snow_c, snow_e, snow_lambda, snow_ro,...
                            snow_k] = update_snow_a1(c_ice, snow_lambda,...
                            precip_out(:,tt), snow_mwe, snow_e, air_T_out(tt), ro_water,...
                            T_melt, ro_snow_initial, snow_ro, snow_mwe, k_ice,...
                            ro_ice, b_exp, c_water, snow_l_ro, k_snow_ice_max);
                    else
                        %update snow variables
                        [snow_c, snow_e, snow_lambda, snow_ro,...
                            snow_k] = update_snow_a1(c_ice, snow_lambda,...
                            precip_out(:,tt), snow_mwe, snow_e, air_T_out(tt), ro_water,...
                            T_melt, ro_snow_initial, snow_ro_e, snow_mwe_e, k_ice,...
                            ro_ice, b_exp, c_water, snow_l_ro, k_snow_ice_max);
                    end
                end

                %surface flux

                [snow_albedo] = snow_albedo_calc(snow_albedo, precip_out(tt),...
                    refresh_albedo, t_step, new_snow_albedo, snow_albedo_min, tau_snow_cold,...
                    tau_snow_melt, snow_T, T_melt);
                snow_albedo = snow_albedo*albedo_mult; if snow_albedo > 1; snow_albedo = 1; end

                [q] = surface_flux1(SW_down_out(tt), albedo_out(tt), air_T_out(tt),...
                    wind_speed_out(tt), pressure_out(tt), LW_down_out(tt), hum_out(tt),...
                    snow_T, snow_albedo, AWS_albedo, Io_ice, e_ice,0);

                %apply surface flux to surface cell
                snow_e = snow_e + q*t_step*3600; %*t_step... to get energy over timestep

                %calcualte temperature and lambda for cells based on enthalpy
                [T(:,1), lambda(:,1), model_stage, snow_T] = temp_lambda_profile1(enthalpy(:,1),...
                    total_grid_num, grid_profile(:,1), T_melt, c_ice, Lf, model_stage,...
                    ro_water, c_water, snow_mwe, snow_e);
                if model_stage == 3
                    snow_mwe = 0;
                end
                %only continue if still model stage 2
                if model_stage == 2

                    %set thermal conductivity and heat capacity of ice cells
                    k(:,1) = lambda(:,1)*k_water + (1 - lambda(:,1))*k_ice;
                    c(:,1) = lambda(:,1)*c_water + (1 - lambda(:,1))*c_ice;

                    %conduction calculations
                    [enthalpy_temp] = conduct_update_enthalpy_new(k(:,1), c(:,1), T(:,1),...
                        grid_profile(:,1), enthalpy(:,1), t_step, total_grid_num, ro_water, lower_boundary,...
                        model_stage, snow_T, snow_k, snow_c, snow_z, snow_e);

                    %seperate enthalpy_temp
                    snow_e = enthalpy_temp(1);
                    enthalpy = enthalpy_temp(2:end);

                    %second round of snow updates
                    if tt > 1
                        [snow_ro, snow_z, snow_lambda, snow_T, snow_ice, snow_l_ro] =...
                            update_snow_b(snow_e, T_melt, snow_mwe, c_ice, ro_water,...
                            Lf, snow_ro, ro_melt_max, t_step, tau_ro, ro_cold_max,...
                            snow_ice_e, snow_lambda_e, snow_l_ro);
                    else
                        [snow_ro, snow_z, snow_lambda, snow_T, snow_ice, snow_l_ro] =...
                            update_snow_b(snow_e, T_melt, snow_mwe, c_ice, ro_water,...
                            Lf, snow_ro, ro_melt_max, t_step, tau_ro, ro_cold_max,...
                            0, 0, snow_l_ro);
                    end
                    track_a = 2;
                else
                    fprintf('Moving from stage 2 to %1.0f. Day %3.0f. \n \n', model_stage, tt*t_step/24)
                end
            end

            %%==============================stage3=============================
            %STAGE 3, lake
            if model_stage == 3 && track_a ~= 3 %if coming from another stage then only do conduction and necessary set up

                %recalcualte temperature and lambda for cells based on enthalpy
                [T(:,1), lambda(:,1), model_stage] = temp_lambda_profile1(enthalpy(:,1),...
                    total_grid_num, grid_profile(:,1), T_melt, c_ice, Lf, model_stage,...
                    ro_water, c_water);

                %calculate lake index for conduction
                [lake_ind, bot_lake_depth, top_lake_ind, top_lake_depth, slush_lid_top,...
                    slush_lid_bot] = lake_index_new(lambda(:,1), ice_lim, total_grid_num, grid_profile(:,1),...
                    tt, t_step, slush_lid_threshold, slush_lid_num);


                %record variables
                lake_depth = bot_lake_depth;
                lake_depth_plot = bot_lake_depth;
                lake_T_av = sum(T(lake_ind,1).*grid_profile(lake_ind))./sum(grid_profile(lake_ind));

                %calculate turbulent heat flux if required
                if lake_depth >= lake_turb_threshold
                    lake_mode = 1; %mode 1 = turbulent flux only, no convection

                    [enthalpy(:,1)] = turbulent_flux1(lake_ind, lake_T_av, air_T_out(tt),...
                        t_step, ro_water, c_water, J, T_melt, lake_depth, Lf,...
                        enthalpy(:,1), grid_profile(:,1), T(:,1), c_ice, lake_mode, model_stage);
                end

                %recalcualte temperature and lambda for cells based on enthalpy
                [T(:,1), lambda(:,1), model_stage] = temp_lambda_profile1(enthalpy(:,1),...
                    total_grid_num, grid_profile(:,1), T_melt, c_ice, Lf, model_stage,...
                    ro_water, c_water);

                %set thermal conductivity and heat capacity
                k(:,1) = lambda(:,1)*k_water + (1 - lambda(:,1))*k_ice;
                c(:,1) = lambda(:,1)*c_water + (1 - lambda(:,1))*c_ice;

                %conduction calculations
                [enthalpy(:,1)] = conduct_update_enthalpy_new(k(:,1), c(:,1), T(:,1),...
                    grid_profile(:,1), enthalpy(:,1), t_step, total_grid_num, ro_water, lower_boundary,...
                    model_stage);

                track_a = 3;

            elseif model_stage == 3 && track_a == 3 %main stage 3 section (not coming from another stage in same time step)

                [T(:,1), lambda(:,1), model_stage] = temp_lambda_profile1(enthalpy(:,1),...
                    total_grid_num, grid_profile(:,1), T_melt, c_ice, Lf, model_stage,...
                    ro_water, c_water, snow_mwe, snow_e, lake_ind);

                [lake_ind, bot_lake_depth, top_lake_ind, top_lake_depth, slush_lid_top,...
                    slush_lid_bot] = lake_index_new(lambda(:,1), ice_lim, total_grid_num, grid_profile(:,1),...
                    tt, t_step, slush_lid_threshold, slush_lid_num);

                if ~isempty(top_lake_ind)
                    lake_depth = top_lake_depth+bot_lake_depth;
                    lake_ind = [top_lake_ind;lake_ind];
                else
                    lake_depth = bot_lake_depth;
                end

                %calculate lake albedo using lake depth
                lake_albedo = 0.1911+exp(-1.0445.*lake_depth-0.9183);
                %albedo multiplier
                lake_albedo = lake_albedo*albedo_mult; if lake_albedo > 1; lake_albedo = 1; end

                %surface flux
                [q,~,~,~,~,~,~,~] = surface_flux1(SW_down_out(tt), albedo_out(tt), air_T_out(tt),...
                    wind_speed_out(tt), pressure_out(tt), LW_down_out(tt), hum_out(tt),...
                    T(1,1), lake_albedo, 0, Io_water, e_ice,3); %AWS_albedo = 0 as lake used

                %calculate Io and tau profiles
                [Io(:,1)] = Io_calc(Io_water, Io_ice, Io_slush, total_grid_num, lambda(:,1));
                [tau(:,1)] = tau_calc(tau_water, tau_slush, tau_ice, total_grid_num, lambda(:,1));

                %SW propagation through water. 0 is to prevent AWS albedo from
                %being used
                [SW_prop(:,1)] = SW_propagate1(SW_down_out(tt), grid_profile(:,1), albedo_out(tt),...
                    t_step, Io(:,1), tau(:,1), lake_albedo, 0, basal_SW_distribute);

                %apply enthalpy changes
                enthalpy(1,1) = enthalpy(1,1) + q*3600*t_step; %*3600*t_step for amount in that timestep
                enthalpy(:,1) = enthalpy(:,1) + SW_prop(:,1);

                [T(:,1), lambda(:,1), model_stage] = temp_lambda_profile1(enthalpy(:,1),...
                    total_grid_num, grid_profile(:,1), T_melt, c_ice, Lf, model_stage,...
                    ro_water, c_water, snow_mwe, snow_e, lake_ind);

                %
                [lake_ind, bot_lake_depth, top_lake_ind, top_lake_depth, slush_lid_top,...
                    slush_lid_bot] = lake_index_new(lambda(:,1), ice_lim, total_grid_num, grid_profile(:,1),...
                    tt, t_step, slush_lid_threshold, slush_lid_num);

                if ~isempty(top_lake_ind)%
                    lake_depth = top_lake_depth+bot_lake_depth;
                    lake_ind = [top_lake_ind;lake_ind];
                    %             lake_depth_plot(tt) = lake_depth(tt)-lake_insert_ind(tt)*ice_grid_z;
                else
                    lake_depth = bot_lake_depth;
                    %             lake_depth_plot(tt) = lake_depth(tt)-lake_insert_ind(tt)*ice_grid_z;
                end

                %otherwise calculate as normal
                lake_T_av = sum(T(lake_ind).*grid_profile(lake_ind))./sum(grid_profile(lake_ind));

                %if lake is above threshold for turbulence calculate flux, calculate fluxes
                %homogenize core temperature if above turbulence threshold, apply
                %turbulence fluxes and update enthalpy
                if lake_depth >= lake_turb_threshold

                    %calculate upper and lower turbulent heat fluxes following
                    %Buzzard (2017)
                    if min(lake_ind) == 1
                        lf_upper = sign(lake_T_av - air_T_out(tt))*ro_water*c_water*J*abs(lake_T_av - air_T_out(tt))^(4/3)*3600*t_step;
                    else
                        lf_upper = sign(lake_T_av - T_melt)*ro_water*c_water*J*abs(lake_T_av - T_melt)^(4/3)*3600*t_step;
                    end
                    lf_lower = sign(lake_T_av - T(max(lake_ind)+1))*ro_water*c_water*J*abs(lake_T_av - T(max(lake_ind)+1))^(4/3)*3600*t_step; %assumes ice-water interface at T_melt

                    %calculate change in temperature of core
                    dT = (-lf_upper - lf_lower)/...
                        (ro_water*c_water*lake_depth);

                    %apply turbulent fluxes as required
                    if min(lake_ind) > 1 %if statement to apply top flux if lid has formed
                        enthalpy(min(lake_ind) - 1) = enthalpy(min(lake_ind) - 1) + lf_upper;
                    end
                    enthalpy(max(lake_ind)+1) = enthalpy(max(lake_ind)+1) + lf_lower;

                    [T(:,1), lambda(:,1), model_stage] = temp_lambda_profile1(enthalpy(:,1),...
                        total_grid_num, grid_profile(:,1), T_melt, c_ice, Lf, model_stage,...
                        ro_water, c_water, snow_mwe, snow_e, lake_ind);

                    %update core temperature
                    lake_T_av = lake_T_av + dT;
                    T(lake_ind) = lake_T_av;
                end

                if lake_depth >= lake_prof_threshold

                    %calculate convection profile and assign to main temperature and
                    %enthalpy profiles
                    [T(lake_ind), enthalpy(lake_ind)] = convection1(T(lake_ind),...
                        grid_profile(lake_ind), ies80, ro_water, c_ice, T_melt, Lf, c_water);

                end

                %========================test=================================
                [lake_ind, bot_lake_depth, top_lake_ind, top_lake_depth, slush_lid_top,...
                    slush_lid_bot] = lake_index_new(lambda(:,1), ice_lim, total_grid_num, grid_profile(:,1),...
                    tt, t_step, slush_lid_threshold, slush_lid_num);
                
                if ~isempty(top_lake_ind)
                    lake_depth(1) = top_lake_depth+bot_lake_depth;
                    lake_ind = [top_lake_ind;lake_ind];
                    %             lake_depth_plot(tt) = lake_depth(tt)-lake_insert_ind(tt)*ice_grid_z;
                else
                    lake_depth(1) = bot_lake_depth;
                    %             lake_depth_plot(tt) = lake_depth(tt)-lake_insert_ind(tt)*ice_grid_z;
                end
                %========================test=================================

                %set thermal conductivity and heat capacity
                k(:,1) = lambda(:,1)*k_water + (1 - lambda(:,1))*k_ice;
                c(:,1) = lambda(:,1)*c_water + (1 - lambda(:,1))*c_ice;

                %conduction calculations
                [enthalpy(:,1)] = conduct_update_enthalpy_new(k(:,1), c(:,1), T(:,1),...
                    grid_profile(:,1), enthalpy(:,1), t_step, total_grid_num, ro_water, lower_boundary,...
                    model_stage, snow_T(1), snow_k(1), snow_c(1), snow_z(1), snow_e(1));

                track_a = 3;
                %display stage change if sufficient days have passed
                if model_stage ~= 3
                    stage_1_end = find(track_a(1,:) == 1, 1, 'last');
                    if ((tt - stage_1_end)*t_step)/24 >= stage_print_hysteresis
                        fprintf('Moving from stage 3 to %1.0f. Day %3.0f. \n \n', model_stage, tt*t_step/24)
                    end
                end
            end

            %%==============================stage4=============================
            %STAGE 4, lake with lid
            if model_stage == 4

                %when the moment comes, and everything is in place, unleash the snow
                if snow_mwe(1)>= snow_threshold && track_b(2) == 0

                    %initialise snow variables
                    [snow_lambda(1), snow_ro(1), snow_c(1), snow_T(1),...
                        snow_e(1), snow_k(1), snow_z(1), snow_ice(1), snow_l_ro(1)] = initialise_snow1(snow_mwe(1), ...
                        ro_snow_initial,c_ice, k_ice, k_air, ro_water, T_melt, air_T_out(tt));

                    snow_albedo = new_snow_albedo*albedo_mult; if snow_albedo > 1; snow_albedo = 1; end
                    track_b(2) = 1; %update tracker

                elseif snow_mwe(1)>= snow_threshold && track_b(2) == 1 %if the snow has been unleashed, proceed to handle it, but proceed with care #bangers (https://olivercoates.bandcamp.com/album/shelleys-on-zenn-la)

                    %update snow variables
                    [snow_c(1), snow_e(1), snow_lambda(1), snow_ro(1),...
                        snow_k(1)] = update_snow_a1(c_ice, snow_lambda(1),...
                        precip_out(:,tt), snow_mwe(1), snow_e(1), air_T_out(tt), ro_water,...
                        T_melt, ro_snow_initial, snow_ro_e, snow_mwe_e, k_ice,...
                        ro_ice, b_exp, c_water, snow_l_ro(1), k_snow_ice_max);
                end


                %surface flux based on snow occurance
                if track_b(2) == 0
                    [q(1)] = surface_flux1(SW_down_out(tt), albedo_out(tt), air_T_out(tt),...
                        wind_speed_out(tt), pressure_out(tt), LW_down_out(tt), hum_out(tt),...
                        T(1,1), ice_albedo, AWS_albedo, Io_ice, e_ice,0);
                elseif track_b(2) == 1
                    [snow_albedo] = snow_albedo_calc(snow_albedo, precip_out(tt),...
                        refresh_albedo, t_step, new_snow_albedo, snow_albedo_min, tau_snow_cold,...
                        tau_snow_melt, snow_T(1), T_melt);
                    snow_albedo = snow_albedo*albedo_mult; if snow_albedo > 1; snow_albedo = 1; end

                    [q(1)] = surface_flux1(SW_down_out(tt), albedo_out(tt), air_T_out(tt),...
                        wind_speed_out(tt), pressure_out(tt), LW_down_out(tt), hum_out(tt),...
                        snow_T(1), snow_albedo, AWS_albedo, Io_ice, e_ice,0);
                end

                %apply surface flux to surface cell
                if track_b(2) == 0
                    enthalpy(1,1) = enthalpy(1,1) + q(1)*t_step*3600; %*t_step... to get energy over timestep
                elseif track_b(2) == 1
                    snow_e(1) = snow_e(1) + q(1)*t_step*3600; %*t_step... to get energy over timestep
                else; fprintf('Error in main, stage 4, surface flux application. \n \n');
                end


                %calcualte temperature and lambda for cells based on enthalpy
                [T(:,1), lambda(:,1), model_stage, snow_T(1)] = temp_lambda_profile1(enthalpy(:,1),...
                    total_grid_num, grid_profile(:,1), T_melt, c_ice, Lf, model_stage,...
                    ro_water, c_water, snow_mwe(1), snow_e(1),lake_ind);

                if model_stage == 5; track_a(1) = 5; track_b(2) = 0; end %to leapfrog temp_lambda_profile that shoots it back to model_stage = 4

                %calculate lake index for conduction
                [lake_ind, bot_lake_depth, top_lake_ind, top_lake_depth, slush_lid_top,...
                    slush_lid_bot] = lake_index_new(lambda(:,1), ice_lim, total_grid_num, grid_profile(:,1),...
                    tt, t_step, slush_lid_threshold, slush_lid_num);

                %                 %record depth and average temperature
                %                 lake_depth(1) = bot_lake_depth;
                if ~isempty(top_lake_ind)
                    lake_depth(1) = top_lake_depth+bot_lake_depth;
                    lake_ind = [top_lake_ind;lake_ind];
                else
                    lake_depth(1) = bot_lake_depth;
                end
                lake_T_av(1) = sum(T(lake_ind,1).*grid_profile(lake_ind,1))./sum(grid_profile(lake_ind,1));
                lid_thick(1) = sum(grid_profile(1:min(lake_ind) - 1,1));

                %calculate turbulent heat flux if required
                if lake_depth(1) >= lake_turb_threshold && lake_depth(1) < lake_prof_threshold
                    lake_mode = 1; %mode 1 = turbulent flux only, no convection

                    [enthalpy(:,1)] = turbulent_flux1(lake_ind, lake_T_av(1), air_T_out(tt),...
                        t_step, ro_water, c_water, J, T_melt, lake_depth(1), Lf,...
                        enthalpy(:,1), grid_profile(:,1), T(:,1), c_ice, lake_mode, model_stage);


                    %recalcualte temperature and lambda for cells based on enthalpy
                    [T(:,1), lambda(:,1), model_stage, snow_T(1)] = temp_lambda_profile1(enthalpy(:,1),...
                        total_grid_num, grid_profile(:,1), T_melt, c_ice, Lf, model_stage,...
                        ro_water, c_water, snow_mwe(1), snow_e(1),lake_ind);

                elseif lake_depth(1) >= lake_prof_threshold %calculate convection profile if required
                    lake_mode = 2; %mode 2 = turbulent flux and convection

                    %calculate turbulent heat flux
                    [~, dT, lf_lower(1), lf_upper(1)] = turbulent_flux1(lake_ind, lake_T_av(1), air_T_out(tt),...
                        t_step, ro_water, c_water, J, T_melt, lake_depth(1), Lf,...
                        enthalpy(:,1), grid_profile(:,1), T(:,1), c_ice, lake_mode, model_stage);

                    %apply temperature change to lake
                    T(lake_ind,1) = T(lake_ind,1) + dT;

                    %calculate convection profile and assign to main temperature and
                    %enthalpy profiles
                    [T(lake_ind,1), enthalpy(lake_ind,1)] = convection1(T(lake_ind,1),...
                        grid_profile(lake_ind,1), ies80, ro_water, c_ice, T_melt, Lf, c_water);

                    %apply upper and lower convective flux
                    enthalpy(max(lake_ind) + 1,1) = enthalpy(max(lake_ind) + 1,1) + lf_lower(1);
                    if min(lake_ind) - 1 >= 1 %to prevent if no lid occurs just before stage 4 switch
                        enthalpy(min(lake_ind) - 1,1) = enthalpy(min(lake_ind) - 1,1) + lf_upper(1);
                    end

                    %recalcualte temperature and lambda for cells based on enthalpy
                    [T(:,1), lambda(:,1), model_stage] = temp_lambda_profile1(enthalpy(:,1),...
                        total_grid_num, grid_profile(:,1), T_melt, c_ice, Lf, model_stage,...
                        ro_water, c_water, snow_mwe(1), snow_e(1),lake_ind);
                    
                end

                if track_a(1) == 5; model_stage = 5;  track_b(2) = 0; snow_z(1) = 0; end %complete leapfrog from initial temp_lambda_profile
                if all(lambda == 0, 'all'); model_stage = 1;end
                %set thermal conductivity and heat capacity of ice cells
                k(:,1) = lambda(:,1)*k_water + (1 - lambda(:,1))*k_ice;
                c(:,1) = lambda(:,1)*c_water + (1 - lambda(:,1))*c_ice;

                %conduction calculations. Model stage held at 4 to prevent a
                %repetitive initialisation section at the start of stage 5
                [enthalpy_temp] = conduct_update_enthalpy_new(k(:,1), c(:,1), T(:,1),...
                    grid_profile(:,1), enthalpy(:,1), t_step, total_grid_num, ro_water, lower_boundary,...
                    4,  snow_T(1), snow_k(1), snow_c(1), snow_z(1), snow_e(1));

                %seperate enthalpy_temp and apply snow updates if required
                if track_b(2) == 0
                    enthalpy(:,1) = enthalpy_temp;
                elseif track_b(2) == 1 %so snow is occuring
                    snow_e(1) = enthalpy_temp(1);
                    enthalpy(:,1) = enthalpy_temp(2:end);

                    %second round of snow updates
                    [snow_ro(1), snow_z(1), snow_lambda(1), snow_T(1), snow_ice(1), snow_l_ro(1)] =...
                        update_snow_b(snow_e(1), T_melt, snow_mwe(1), c_ice, ro_water,...
                        Lf, snow_ro(1), ro_melt_max, t_step, tau_ro, ro_cold_max,...
                        snow_ice_e, snow_lambda_e, snow_l_ro(1));

                else; fprintf('Error in main, stage 4 enthalpy allocation. \n \n');
                end

                track_a = 4; %model stage = 2


                if model_stage == 5
                    fprintf('Moving from stage 4 to %1.0f. Day %3.0f. \n \n', model_stage, tt*t_step/24)

                    %reset snow trackers
                    track_b(1:2) = 0;
                end
            end

            %%==============================stage5=============================
            %STAGE 5, lid breakup
            if model_stage == 5

                %calculate tau profile
                [tau(:,1)] = tau_calc(tau_water, tau_slush, tau_ice, total_grid_num, lambda(:,1));

                %calculate lake albedo
                lake_albedo(1) = 0.1911+exp(-1.0445.*lake_depth(1)-0.9183);                %albedo multiplier
                lake_albedo(1) = lake_albedo(1)*albedo_mult; if lake_albedo(1) > 1; lake_albedo(1) = 1; end

                %surface flux
                [q(1)] = surface_flux1(SW_down_out(tt), albedo_out(tt), air_T_out(tt),...
                    wind_speed_out(tt), pressure_out(tt), LW_down_out(tt), hum_out(tt),...
                    T(1,1), lake_albedo(1), 0, Io_water, e_ice,3); %AWS_albedo = 0 as lake used

                %calculate Io_profile
                [Io(:,1)] = Io_calc(Io_water, Io_ice, Io_slush, total_grid_num, lambda(:,1));

                %SW propagation through water. 0 is to prevent AWS albedo from
                %being used
                [SW_prop(:,1)] = SW_propagate1(SW_down_out(tt), grid_profile(:,1), albedo_out(tt),...
                    t_step, Io(:,1), tau(:,1), lake_albedo(1), 0, basal_SW_distribute);

                %apply incoming energy flux to enthalpy array
                enthalpy(1,1) = enthalpy(1,1) + q(1)*3600*t_step;
                enthalpy(:,1) = enthalpy(:,1) + SW_prop(:,1);

                %calcualte temperature and lambda for cells based on enthalpy
                [T(:,1), lambda(:,1), model_stage, snow_T(1)] = temp_lambda_profile1(enthalpy(:,1),...
                    total_grid_num, grid_profile(:,1), T_melt, c_ice, Lf, model_stage,...
                    ro_water, c_water, snow_mwe(1), snow_e(1),lake_ind);

                %calculate lake index for conduction
                %
                [bot_lake_ind, bot_lake_depth, top_lake_ind, top_lake_depth, slush_lid_top,...
                    slush_lid_bot] = lake_index_new(lambda(:,1), ice_lim, total_grid_num, grid_profile(:,1),...
                    tt, t_step, slush_lid_threshold, slush_lid_num);

                %see if lid breakup should occur if it hasn't already
                if track_b(3) == 0
                    [track_b(3)] = lid_instability(lambda(:,1), grid_profile(:,1), slush_lid_top,...
                        slush_lid_bot, slush_mech_threshold, breakup_threshold, ro_water, ro_ice);
                end
                if isempty(top_lake_ind)
                    lake_ind = bot_lake_ind;
                else
                    lake_ind = [top_lake_ind;bot_lake_ind];
                end
                %lid may also breakup if lake_index subfunction combines lakes into
                %one. If this is the case, prompt lid breakup and reobtain lake
                %index with a greater slush_lid_threshold
                if (min(lake_ind) - 1) <= 1
                    track_b(3) = 1;

                    [lake_ind, bot_lake_depth, top_lake_ind, top_lake_depth, slush_lid_top,...
                        slush_lid_bot] = lake_index_new(lambda(:,1), ice_lim, total_grid_num, grid_profile(:,1),...
                        tt, t_step, 0.5, slush_lid_num);
                end

                %record lake values
                surface_lake_depth(1) = top_lake_depth; %depth of the surface lake
                lake_depth(1) = bot_lake_depth; %depth of the main lake
                sur_lake_bot = max(top_lake_ind) + 1; %bottom of the surface lake
                lake_top = min(lake_ind) - 1; %top of the main lake
                lake_bot = max(lake_ind) + 1; %bottom of the main lake

                %if lid is unstable
                if track_b(3) == 1

                    %lid breakup and lake combine
                    [model_stage, lake_ind, enthalpy(:,1), track_b, track_a(1),...
                        lake_depth(1)] = lid_breakup(track_b, enthalpy(:,1),...
                        lake_ind, grid_profile(:,1));

                    %updates
                    fprintf('Moving from stage 5 to 3. Day %3.0f. \n \n', tt*t_step/24)
                    %                 T(:,tt + 1) = T(:,tt); %skip temperature calculation and just update

                else %calculate average temperatures as normal

                    sur_lake_T_av(1) = sum(grid_profile(top_lake_ind,1).*T(top_lake_ind,1))./...
                        sum(grid_profile(top_lake_ind,1)); %surface lake average temperature
                    lake_T_av(1) = sum(grid_profile(lake_ind,1).*T(lake_ind,1))./...
                        sum(grid_profile(lake_ind,1)); %main lake average temperature
                end

                snow_z(1) = 0;
                snow_mwe(1) = 0;
                %only continue if lid has not yet disintegrated
                if model_stage == 5 || model_stage == 4

                    %calculate turbulent heat flux for surface lake if required.
                    %convection profile is not considered here
                    if surface_lake_depth(1) >= lake_turb_threshold

                        %calculate upper and lower turbulent heat fluxes following
                        %Buzzard (2017)
                        if min(top_lake_ind) == 1
                            lf_upper(1) = sign(sur_lake_T_av(1) - air_T_out(tt))*ro_water*c_water*J*abs(sur_lake_T_av(1) - air_T_out(tt))^(4/3)*3600*t_step;
                        else
                            lf_upper(1) = sign(sur_lake_T_av(1) - T_melt)*ro_water*c_water*J*abs(sur_lake_T_av(1) - T_melt)^(4/3)*3600*t_step;
                        end
                        lf_lower(1) = sign(sur_lake_T_av(1) - T(max(lake_ind)+1,1))*ro_water*c_water*J*abs(sur_lake_T_av(1) - T(max(lake_ind)+1,1))^(4/3)*3600*t_step; %assumes ice-water interface at T_melt

                        %calculate change in temperature of core
                        dT = (-lf_upper(1) - lf_lower(1))/...
                            (ro_water*c_water*surface_lake_depth(1));

                        %apply turbulent fluxes as required
                        if min(top_lake_ind) > 1 %if statement to apply top flux if lid has formed
                            enthalpy(min(top_lake_ind) - 1,1) = enthalpy(min(top_lake_ind) - 1,1) + lf_upper(1);
                        end
                        enthalpy(max(top_lake_ind)+1,1) = enthalpy(max(top_lake_ind)+1,1) + lf_lower(1);

                        %calcualte temperature and lambda for cells based on enthalpy
                        [T(:,1), lambda(:,1), model_stage] = temp_lambda_profile1(enthalpy(:,1),...
                            total_grid_num, grid_profile(:,1), T_melt, c_ice, Lf, model_stage,...
                            ro_water, c_water, snow_mwe(1), snow_e(1),lake_ind);

                        %update core temperature
                        sur_lake_T_av(1) = sur_lake_T_av(1) + dT;
                        T(top_lake_ind,1) = sur_lake_T_av(1);

                    end

                    %calculate turbulent heat flux for main lake if required
                    if lake_depth(1) >= lake_turb_threshold && lake_depth(1) < lake_prof_threshold
                        lake_mode = 1; %mode 1 = turbulent flux only, no convection

                        [enthalpy(:,1)] = turbulent_flux1(lake_ind, lake_T_av(1), air_T_out(tt),...
                            t_step, ro_water, c_water, J, T_melt, lake_depth(1), Lf,...
                            enthalpy(:,1), grid_profile(:,1), T(:,1), c_ice, lake_mode, model_stage);

                        %recalcualte temperature and lambda for cells based on enthalpy
                        [T(:,1), lambda(:,1), model_stage] = temp_lambda_profile1(enthalpy(:,1),...
                            total_grid_num, grid_profile(:,1), T_melt, c_ice, Lf, model_stage,...
                            ro_water, c_water, snow_mwe(1),snow_e(1),lake_ind);


                    elseif lake_depth(1) >= lake_prof_threshold %calculate convection profile for main lake if required
                        lake_mode = 2; %mode 2 = turbulent flux and convection

                        %calculate turbulent heat flux
                        [~, dT, lf_lower(1)] = turbulent_flux1(lake_ind, lake_T_av(1), air_T_out(tt),...
                            t_step, ro_water, c_water, J, T_melt, lake_depth(1), Lf,...
                            enthalpy(:,1), grid_profile(:,1), T(:,1), c_ice, lake_mode, model_stage);

                        %apply temperature change to lake
                        T(lake_ind,1) = T(lake_ind,1) + dT;


                        %calculate convection profile and assign to main temperature and
                        %enthalpy profiles
                        [T(lake_ind,1), enthalpy(lake_ind,1)] = convection1(T(lake_ind,1),...
                            grid_profile(lake_ind,1), ies80, ro_water, c_ice, T_melt, Lf, c_water);

                        %apply lower convective flux
                        enthalpy(max(lake_ind) + 1,1) = enthalpy(max(lake_ind) + 1,1) + lf_lower(1);

                        %a recalc_ind is used as if updating the whole profile the
                        %old enthalpy updates the new surface temperature
                        recalc_ind = [min(lake_ind) - 1; lake_ind; max(lake_ind + 1)];
                        if model_stage == 4; track_a(1) = 4; end %leapfrog stage 4 to avoid reset in temp_lambda_profile due to recalc_ind

                        %recalcualte temperature and lambda for lake cells based on enthalpy
                        [T(recalc_ind,1), lambda(recalc_ind,1), model_stage] = temp_lambda_profile1(enthalpy(recalc_ind,1),...
                            numel(recalc_ind), grid_profile(recalc_ind,1), T_melt, c_ice, Lf, model_stage,...
                            ro_water, c_water, snow_mwe(1), snow_e(1),lake_ind);

                        if track_a(1) == 4; model_stage = 4; end %copmlete leapfrog

                    end

                    %set thermal conductivity and heat capacity of ice cells
                    k(:,1) = lambda(:,1)*k_water + (1 - lambda(:,1))*k_ice;
                    c(:,1) = lambda(:,1)*c_water + (1 - lambda(:,1))*c_ice;
                    if model_stage == 4
                        track_b(2) = 0;
                    end
                    %conduction calculations
                    snowz = 0;
                    [enthalpy(:,1)] = conduct_update_enthalpy_new(k(:,1), c(:,1), T(:,1),...
                        grid_profile(:,1), enthalpy(:,1), t_step, total_grid_num, ro_water, lower_boundary,...
                        model_stage, snow_T(1), snow_k(1), snow_c(1), snowz, snow_e(1));
                end

                track_a = model_stage;
            end

            
            Track_b{m,n} = track_b;
            Grid_profile{m,n} = grid_profile;

            if size(Lambda{m,n},1)==length(lambda)
                Lambda{m,n}(:,tt) = lambda;
                Enthalpy{m,n}(:,tt) = enthalpy;
                T_total{m,n}(:,tt) = T;
            else
                lambda_add = nan(length(lambda)-size(Lambda{m,n},1),tt-1);
                Lambda{m,n} = [lambda_add;Lambda{m,n}];
                Lambda{m,n}(:,tt) = lambda;
                enthalpy_add = nan(size(lambda_add));
                Enthalpy{m,n} = [enthalpy_add;Enthalpy{m,n}];
                Enthalpy{m,n}(:,tt) = enthalpy;
                T_add = nan(size(lambda_add));
                T_total{m,n} = [T_add;T_total{m,n}];
                T_total{m,n}(:,tt) = T;
            end
            Track_a(m,n,tt) = track_a;
            Snow_T(m,n,tt) = snow_T;
            Snow_e(m,n,tt) = snow_e;
            Snow_lambda(m,n,tt) = snow_lambda;
            Snow_mwe(m,n,tt) = snow_mwe;
            Snow_ro(m,n,tt) = snow_ro;
            Snow_l_ro(m,n,tt) = snow_l_ro;
            Snow_z(m,n,tt) = snow_z;
            Snow_ice(m,n,tt) = snow_ice;
            Lake_depth(m,n,tt) = lake_depth;
            Lake_albedo(m,n,tt) = lake_albedo;
            Model_stage(m,n,tt) = model_stage;
            Total_grid_num(m,n) = total_grid_num;
            Lake_ind(m,n,tt) = length(lake_ind);
            Q(m,n,tt) = q;



            % bottom ablation
            if ismember(model_stage, [3 4 5])
                if tt == s_tnum
                    bm(m,n,tt) = lake_depth - ice_grid_z*ex_water(m,n);
                else
                    bm(m,n,tt) = lake_depth - Lake_depth(m,n,tt-1) - ice_grid_z*ex_water(m,n);
                end
            end

            % lid depth
            if model_stage == 4 && ~isempty(lake_ind)
                lid_depth(m,n,tt) = length(1:lake_ind-1)*ice_grid_z+snow_z;% If there are two lakes here, we will have to make changes later!!!
            elseif model_stage == 4 && isempty(lake_ind)
                lid_depth(m,n,tt) = ice_grid_z+snow_z;
            end
        end

        % lake area
        if any(ismember(unique(Model_stage(:,:,tt)), 3))
            lake_area(tt,1) = length(find(squeeze(Model_stage(:,:,tt))==3))*cell_size^2;
            lake_vol(tt,1) = sum(sum(squeeze(Lake_depth(:,:,tt)).*cell_size^2));
        end

        if mod(tt,24)==0
            disp(tt/24)
        end
        if any(Model_stage(:,:,tt) == 4, 'all')
            break;
        end
    end

    toc; fprintf('\n \n')
    lake_vols(:,cci) = lake_vol(doyz.*24-10)./1000000000;
    lake_areas(:,cci) = lake_area(doyz.*24-10)./1000000;
    R2_vols(cci,1) = 1-(sum((lake_vols(:,cci)-volsz).^2)/sum((volsz-mean(volsz)).^2));
    R2_areas(cci,1) = 1-(sum((lake_areas(:,cci)-areaz).^2)/sum((areaz-mean(areaz)).^2));
    rmse_vols(cci,1)=sqrt(mean((lake_vols(:,cci)-volsz).^2));
    rmse_areas(cci,1)=sqrt(mean((lake_areas(:,cci)-areaz).^2));
end

% INDEX = [cc',R2_areas,R2_vols,rmse_areas,rmse_vols];

% 
% % output
% % Projection information
% info = geotiffinfo(demname);
% for i = 1:length(doyz)
%         lakedepth = squeeze(Lake_depth(:,:,doyz(i)*24-10));
%         % 
%         outputFile = ['output\L3_',num2str(doyz(i)),'.tif'];
%         geotiffwrite(outputFile, lakedepth, R, ...
%             'CoordRefSysCode', info.GeoTIFFCodes.PCS);
% end
Lake_depth_plot = zeros([size(Lambda),tt]);
for i = 1:size(Lambda,1)
    for j = 1:size(Lambda,2)
        if ~isempty(Lambda{i,j})
        for z = 1:tt
            Lake_depth_plot(i,j,z) = nansum(Lambda{i,j}(:,z))*0.1;
        end
        end
    end
end
save('L3_2021.mat') %L1:L1_2021.mat; L2:L2_2020.mat




