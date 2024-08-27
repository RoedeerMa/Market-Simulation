function [srevenue, soutput, energy_price, genergy_ratio] = ed_run_SM(case_config,a_parameter)
    T = case_config.time_unit_resolution;
    day_num = case_config.time_unit_num;

    x_gen_wind = sdpvar(1,T);
    x_gen_pv = sdpvar(1,T);
    x_gen_coal_high = sdpvar(1,T);
    x_gen_coal_low = sdpvar(1,T);
    x_gen_gas = sdpvar(1,T);
    x_gen_nuclear = sdpvar(1,T);
    x_gen_hydro = sdpvar(1,T);
    x_gen_dc = sdpvar(1,T);
    x_load_shed = sdpvar(1,T);
    
    x_storage_charge = sdpvar(1,T);
    x_storage_discharge = sdpvar(1,T);
    x_storage_soc = sdpvar(1,T);
    x_p2h_charge = sdpvar(1,T);
    x_p2h_discharge = sdpvar(1,T);

    load_peak = case_config.load_peak(1);
    
    pv_capacity = case_config.pv_penetration * load_peak;
    wind_capacity = case_config.wind_penetration  * load_peak;
    coal_high_capacity = case_config.coal_high_penetration * load_peak;
    coal_low_capacity = case_config.coal_low_penetration * load_peak;
    gas_capacity = case_config.gas_penetration  * load_peak;
    nuclear_capacity = case_config.nuclear_penetration * load_peak;
    hydro_capacity = case_config.hydro_penetration * load_peak;
    storage_capacity = case_config.storage_penetration * load_peak;
    p2h_capacity = case_config.p2h_penetration * load_peak;
    dc_capacity = case_config.dc_penetration * load_peak;
    
    storage_eta = case_config.storage_eta;
    storage_duration = case_config.storage_duration;

    load_pu = a_parameter(:,1)';
    pv_pu = a_parameter(:,2)';
    wind_pu = a_parameter(:,3)';
    
    
    reserve_segment_upper = case_config.reserve_insufficient_penalty(:,2);
    reserve_segment_lower = case_config.reserve_insufficient_penalty(:,1);
    reserve_segment_penalty = case_config.reserve_insufficient_penalty(:,3);
    reserve_segment_num = size(reserve_segment_penalty,1);
    
    reserve_req = case_config.reserve_req_ratio * load_peak * load_pu;
    
    x_reserve_insufficient_segment = sdpvar(reserve_segment_num,T);
    x_reserve_storage = sdpvar(1,T);
    x_reserve_p2h = sdpvar(1,T);
    %% constraints
    constraints = [x_gen_wind + x_gen_pv + x_gen_coal_high + x_gen_coal_low + x_gen_gas + x_gen_hydro + x_gen_nuclear + x_gen_dc + x_load_shed + x_storage_discharge + x_p2h_discharge - x_p2h_charge - x_storage_charge - load_peak * load_pu == zeros(1,T)];
    
    constraints = [constraints,x_load_shed <= load_peak * load_pu];
    constraints = [constraints,x_load_shed >= 0];  
    % wind
    constraints = [constraints,x_gen_wind >= 0];
    constraints = [constraints,x_gen_wind <= wind_capacity .* wind_pu];
    % pv
    constraints = [constraints,x_gen_pv >= 0];
    constraints = [constraints,x_gen_pv <= pv_capacity .* pv_pu];
    % gas
    constraints = [constraints,x_gen_gas >= 0];
    constraints = [constraints,x_gen_gas <= gas_capacity];
    % coal 
    constraints = [constraints,x_gen_coal_high >= 0];
    constraints = [constraints,x_gen_coal_high <= coal_high_capacity];
    constraints = [constraints,x_gen_coal_low >= 0];
    constraints = [constraints,x_gen_coal_low <= coal_low_capacity];
    % nuclear
    constraints = [constraints,x_gen_nuclear >= case_config.nuclear_min * nuclear_capacity];
    constraints = [constraints,x_gen_nuclear <= nuclear_capacity];
    % hydro
    constraints = [constraints,x_gen_hydro >= 0];
    constraints = [constraints,x_gen_hydro <= hydro_capacity];
    % dc import
    constraints = [constraints,x_gen_dc >= 0];
    constraints = [constraints,x_gen_dc <= dc_capacity];
    % battery storage 
    day_start_idx = 1:24:T;
    day_end_idx = 24:24:T;
    constraints = [constraints,x_storage_charge >= 0];
    constraints = [constraints,x_storage_discharge >= 0];
    constraints = [constraints,x_storage_charge <= storage_capacity];
    constraints = [constraints,x_storage_discharge <= storage_capacity];  
    constraints = [constraints,x_storage_soc(2:end) == x_storage_soc(1:end-1) + x_storage_charge(1:end-1) * storage_eta - x_storage_discharge(1:end-1) / storage_eta];
    constraints = [constraints,x_storage_soc(day_start_idx) == x_storage_soc(day_end_idx) + x_storage_charge(day_end_idx) * storage_eta - x_storage_discharge(day_end_idx) / storage_eta];
    
    constraints = [constraints,x_storage_soc <= storage_duration * storage_capacity];
    constraints = [constraints,x_storage_soc >= 0];
    
    constraints = [constraints,x_reserve_storage <= storage_capacity - x_storage_charge]; 
    constraints = [constraints,x_reserve_storage <= storage_capacity - x_storage_discharge]; 
    constraints = [constraints,x_reserve_storage >= 0];
    
    % P2H storage
    constraints = [constraints,x_p2h_charge >= 0];
    constraints = [constraints,x_p2h_charge <= p2h_capacity];
    constraints = [constraints,x_p2h_discharge >= 0];
    constraints = [constraints,x_p2h_discharge <= p2h_capacity];
    constraints = [constraints,sum(x_p2h_charge,'all') == sum(x_p2h_discharge,'all')];
    constraints = [constraints,x_reserve_p2h <= p2h_capacity - x_p2h_discharge];
    constraints = [constraints,x_reserve_p2h <= p2h_capacity - x_p2h_charge];
    constraints = [constraints,x_reserve_p2h >= 0];

    % reserve
    reserve = (coal_high_capacity - x_gen_coal_high) + (coal_low_capacity - x_gen_coal_low) +(gas_capacity - x_gen_gas) + (hydro_capacity - x_gen_hydro) + x_reserve_storage + x_reserve_p2h;
    
    % ORDC 
    constraints = [constraints, reserve + ones(1,size(x_reserve_insufficient_segment,1)) * x_reserve_insufficient_segment >= reserve_req];
    constraints = [constraints, x_reserve_insufficient_segment >= 0];
    constraints = [constraints, x_reserve_insufficient_segment <= (reserve_segment_upper  - reserve_segment_lower) * reserve_req];
    
    %% objective
    gen_coal_high_cost = case_config.gen_cost(1,2);
    gen_coal_low_cost = case_config.gen_cost(2,2);
    gen_gas_cost = case_config.gen_cost(3,2);
    gen_nuclear_cost = case_config.gen_cost(4,2);
    gen_hydro_cost = case_config.gen_cost(5,2);
    p2h_cost = case_config.gen_cost(6,2);
    load_voll = case_config.load_voll;
    
    reserve_social_welfare_reduction = sum(diag(reserve_segment_penalty * load_voll) * x_reserve_insufficient_segment,'all');
    obj = sum (gen_coal_high_cost * x_gen_coal_high ,'all') + sum (gen_coal_low_cost * x_gen_coal_low ,'all') + sum (gen_nuclear_cost * x_gen_nuclear ,'all')  ...
        + sum (gen_hydro_cost * x_gen_hydro ,'all') + sum (gen_gas_cost * x_gen_gas ,'all') + sum(p2h_cost * x_p2h_discharge, 'all')...
        + load_voll * sum(x_load_shed, 'all') + reserve_social_welfare_reduction;

    opts = sdpsettings('solver','gurobi','verbose',0);
    optimize(constraints,obj,opts);
    
    energy_price = -dual(constraints(1))';
    
    wind_output = value(x_gen_wind);
    pv_output = value(x_gen_pv);
    
    soutput.wind = wind_output;
    soutput.pv = pv_output;
    
    srevenue_wind = wind_output * energy_price;
    srevenue_pv = pv_output * energy_price;
    
    srevenue.wind = srevenue_wind;
    srevenue.pv = srevenue_pv;
    
    genergy_ratio = (sum(value(x_gen_wind),'all') + sum(value(x_gen_pv),'all'))/(sum(load_pu) * load_peak);

    yalmip('clear');

end