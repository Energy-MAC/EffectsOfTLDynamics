using PowerSimulationsDynamics
using PowerSystems
using Parameters
using EffectsOfTLDynamics

const ETL = EffectsOfTLDynamics


function get_small_signal_results_from_samples(file_name, df, gfm_bus, gfl_bus, sm_bus, dyn_line, multiseg, fd)

    # Create new df columns 
    df[!,:stable] .= NaN
    df[!,:gfm_eta] .= NaN
    df[!,:gfl_eta] .= NaN
    df[!,:gfm_eta_t] .= NaN
    df[!,:gfl_eta_t] .= NaN
    df[!,:inv_share_actual] .= NaN
    df[!,:sm_eta] .= NaN
    df[!,:status] .= ""

    for n = 1:size(df)[1];
        # Retrieve parameter values 
        kq = df[n,:kq]
        Ta = df[n,:Ta]
        kd = df[n,:kd]
        kpv = df[n,:kpv]
        kiv = df[n,:kiv]
        kpc = df[n,:kpc]
        kic = df[n,:kic]
        inv_share = df[n,:inv_share]
        gfm_share = df[n, :gfm_share]
        load_scale = df[n,:load_scale]
        line_scale = df[n, :line_scale]

        # Create system  
        sys, p1, p3, perturbation1, perturbation3 = create_sys_with_params(file_name, kpv, kiv, Ta, kd, kq, kpc, kic, load_scale, line_scale, inv_share, gfm_share, gfm_bus, gfl_bus);
        if fd == true; # MSMB only - use p3 
            stab, gfm_eta, gfl_eta, sm_eta, status = run_analysis(sys, dyn_line, multiseg, p3, perturbation1,  gfm_bus, gfl_bus, sm_bus)
        else
            stab, gfm_eta, gfl_eta, sm_eta, status = run_analysis(sys, dyn_line, multiseg, p1, perturbation1,  gfm_bus, gfl_bus, sm_bus)
        end

        # Update df with results 
        df[n,:stable] = stab 
        df[n,:gfm_eta] = gfm_eta
        df[n,:gfl_eta] = gfl_eta
        df[n,:gfm_eta_t] = inv_share*gfm_share
        df[n,:gfl_eta_t] = inv_share*(1-gfm_share)
        df[n,:inv_share_actual] = gfm_eta + gfl_eta
        df[n,:sm_eta] = sm_eta 
        df[n,:status] = string(status)       
    end
 
    complete_df = df
    df = remove_nans_from_df(df)
    stable_df = filter(:stable => x -> x<0, df)
    unstable_df = filter(:stable => x -> x>0, df)
    
    return df, stable_df, unstable_df, complete_df

end

function run_analysis(sys, dyn_lines, multi_seg, p, perturbation, gfm_bus, gfl_bus, sm_bus)

    tspan = (0.0, 0.25)
    if multi_seg == false
        sys = ETL.build_new_impedance_model!(sys, p, dyn_lines, "")
    else
        sys = ETL.build_seg_model!(sys, p, dyn_lines, "")
    end
    sim = build_sim(sys, tspan, perturbation, dyn_lines, p);
    gfm_eta, gfl_eta, sm_eta = calculate_etas(sim.sys, gfm_bus, gfl_bus, sm_bus)
    stab = save_max_nonzero_eig!(sim, [])[1]

    return stab, gfm_eta, gfl_eta, sm_eta, sim.status 
end

function create_sys_with_params(file_name, kpv, kiv, Ta, kd, kq, kpc, kic, load_scale, line_scale, inv_share, gfm_share, gfm_bus, gfl_bus)

    sys = get_system(file_name)
    # Get generators 
    gfm = get_component(Generator, sys, gfm_bus)
    gfl = get_component(Generator, sys, gfl_bus)

    p1 = p1_9bus;
    p3 = p3_9bus;

    # Adjust load scale 
    p1.load_scale = load_scale;
    p3.load_scale = load_scale;

    # Adjust line scale 
    p1.line_scale = line_scale;
    p3.line_scale = line_scale;
    
    # Scale loads according to load scale 
    total_p_load = 0.0;
    for l in get_components(PSY.StandardLoad, sys)
        transform_load_to_constant_impedance(l)
        l.impedance_active_power = l.impedance_active_power * p1.load_scale 
        l.impedance_reactive_power = l.impedance_reactive_power * p1.load_scale 
        total_p_load += l.impedance_active_power
    end

    total_p_load = total_p_load*sys.units_settings.base_value;

    gfm_share_MW = inv_share*gfm_share*total_p_load; # in MW 
    gfl_share_MW = inv_share*(1.0 - gfm_share)*total_p_load; # MW 

    #Adjust GFM gains  
    gfm.dynamic_injector.outer_control.reactive_power_control.kq = kq
    gfm.dynamic_injector.inner_control.kpc = kpc
    gfm.dynamic_injector.inner_control.kpv = kpv
    gfm.dynamic_injector.inner_control.kic = kic
    gfm.dynamic_injector.inner_control.kiv = kiv
    gfm.dynamic_injector.outer_control.active_power_control.kd = kd
    gfm.dynamic_injector.outer_control.active_power_control.Ta = Ta

    #Adjust power injections 
    gfm.active_power = gfm_share_MW/gfm.base_power;
    gfl.active_power = gfl_share_MW/gfl.base_power;

    p1.perturbation_params.crc_params = CRCParam(DynamicGenerator, sm_bus, :V_ref, 0.95)
    p3.perturbation_params.crc_params = CRCParam(DynamicGenerator, sm_bus, :V_ref, 0.95)

    perturbation1 = choose_disturbance(sys, p1.perturbation, p1)
    perturbation3 = choose_disturbance(sys, p3.perturbation, p3)

    return sys, p1, p3, perturbation1, perturbation3

end

function remove_nans_from_df(df)
    filter(:stable => x -> !any(f -> f(x), (ismissing, isnan)), df)
end

function calculate_etas(sys_m, gfm_bus, gfl_bus, sm_bus)
    gfm = get_component(Generator, sys_m, gfm_bus);
    sm = get_component(Generator, sys_m, sm_bus);
    gfl = get_component(Generator, sys_m, gfl_bus);
    # Account for the fact that each machine has a different base power 
    total_gen = gfm.active_power*gfm.base_power + gfl.active_power*gfl.base_power + sm.active_power*sm.base_power # in MW 
    eta_gfm = gfm.active_power*gfm.base_power/total_gen; 
    eta_gfl = gfl.active_power*gfl.base_power/total_gen;
    eta_sm = sm.active_power*sm.base_power/total_gen;
    return eta_gfm, eta_gfl, eta_sm
end

function get_system(file_name)
    sys = System(joinpath(pwd(), file_name));
    return sys 
end 

## PARAMETERS ## 

# Simulation parameters
sim_p = SimParams(
    abstol = 1e-13,
    reltol = 1e-10,
    maxiters = Int(1e10),
    dtmax = 1e-4,
    solver = "Rodas4",
    t_max = 0.1,
)

# Line data 
impedance_csv = "../data/cable_data/dommel_data.csv"
capacitance_csv = "../data/cable_data/dommel_data_C.csv"

M = 1;
factor_z = 1.0; 
factor_y = 1.0;
z_km_1, y_km_1, Z_c_abs_1, z_km_ω_1, z_km_ω_5_to_1_1, Z_c_5_to_1_abs_1 = get_line_parameters(impedance_csv, capacitance_csv, M, factor_z, factor_y);

M = 3;
z_km_3, y_km_3, Z_c_abs_3, z_km_ω_3, z_km_ω_5_to_1_3, Z_c_5_to_1_abs_3 = get_line_parameters(impedance_csv, capacitance_csv, M, factor_z, factor_y);

perturbation_type = "CRC"
perturbation_params = PerturbationParams(t_fault=0.25)
perturbation_params.crc_params = CRCParam(DynamicInverter, "generator-2-1", :V_ref, 0.95)

p1_9bus = ETL.ExpParams(
    nothing, # N 
    1, # M
    100.0, # gets overwritten
    10.0, # lseg  
    Z_c_abs_1, 
    z_km_1,
    y_km_1,
    z_km_ω_1,
    z_km_ω_5_to_1_1,
    Z_c_5_to_1_abs_1,
    default_9_bus_line_dict,
    sim_p, 
    perturbation_type, 
    perturbation_params,
    0.0, #pload 
    0.0, # qload 
    1.0, # line_scale 
    1.0 # load_scale 
);

p3_9bus = ETL.ExpParams(
    nothing, # N 
    3, # M
    100.0, # gets overwritten
    10.0, # lseg  
    Z_c_abs_3, 
    z_km_3,
    y_km_3,
    z_km_ω_3,
    z_km_ω_5_to_1_3,
    Z_c_5_to_1_abs_3,
    default_9_bus_line_dict,
    sim_p, 
    perturbation_type, 
    perturbation_params,
    0.0, # pload 
    0.0, # qload 
    1.0, # line_scale 
    1.0 # load_scale 
);

function bin_data(df, sym)
    # bins data into 5 bins between 0.0 and 1.0 
    x1 = filter(sym => n -> n < 0.2, df)
    x1[!,sym] = 0.1*ones(size(x1)[1]);
    x2 = filter(sym => n -> 0.2 < n < 0.4, df)
    x2[!,sym] = 0.3*ones(size(x2)[1]);
    x3 = filter(sym => n -> 0.4 < n < 0.6, df)
    x3[!,sym] = 0.5*ones(size(x3)[1]);
    x4 = filter(sym => n -> 0.6 < n < 0.8, df)
    x4[!,sym] = 0.7*ones(size(x4)[1]);
    x5 = filter(sym => n -> 0.8 < n < 1.0, df)
    x5[!,sym] = 0.9*ones(size(x5)[1]);

    mod_df = vcat(x1,x2,x3,x4,x5)
    return mod_df
end
