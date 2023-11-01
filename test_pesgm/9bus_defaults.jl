using PowerSimulationsDynamics
using PowerSystems
using Parameters
using EffectsOfTLDynamics

const ETL = EffectsOfTLDynamics

## FUNCTIONS ##

function get_small_signal_results_from_samples(file_name, df, gfm_bus, gfl_bus, sm_bus, dyn_line, multiseg, fd)

    # Create new df columns 
    df[!,:stable] .= NaN
    df[!,:gfm_eta] .= NaN
    df[!,:gfl_eta] .= NaN
    df[!,:sm_eta] .= NaN
    df[!,:status] .= ""

    for n = 1:size(df)[1];
        # Retrieve parameter values 
        kpv = df[n,:kpv]
        kiv = df[n,:kiv]
        kq = df[n,:kq]
        kpc = df[n,:kpc]
        kic = df[n,:kic]
        inv_share = df[n,:inv_share]
        gfm_share = df[n, :gfm_share]
        load_scale = df[n,:load_scale]

        # Create system  
        sys, p1, p3, perturbation1, perturbation3 = create_sys_with_params(file_name, kpv, kiv, kq, kpc, kic, load_scale, inv_share, gfm_share, gfm_bus, gfl_bus);
        if fd == true; # MSMB only - use p3 
            stab, gfm_eta, gfl_eta, sm_eta, status = run_analysis(sys, dyn_line, multiseg, p3, perturbation1,  gfm_bus, gfl_bus, sm_bus)
        else
            stab, gfm_eta, gfl_eta, sm_eta, status = run_analysis(sys, dyn_line, multiseg, p1, perturbation1,  gfm_bus, gfl_bus, sm_bus)
        end

        # Update df with results 
        df[n,:stable] = stab 
        df[n,:gfm_eta] = gfm_eta
        df[n,:gfl_eta] = gfl_eta
        df[n,:sm_eta] = sm_eta 
        df[n,:status] = string(status)       
    end
 
    stable_df = filter(:stable => x -> x<0, df)
    unstable_df = filter(:stable => x -> x>0, df)
    
    return df, stable_df, unstable_df

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

function pesgm_build_sim(file_name, alg_df, dyn_df, mssb_df, msmb_df, gfm_bus, gfl_bus, sm_bus)

    # Sample parameters 
    kpv = rand(kpv_nrel_range); # inner loop voltage gain 
    kiv = rand(kiv_nrel_range); # inner loop integral gain 
    kq = rand(gfm_kq_range); # outer loop reactive power control gain

    kpc = rand(kpc_range); # Inner current loop gain - GFM 
    kic = rand(kic_range); # Inner current loop gain - GFM

    load_scale = rand(load_scale_range)
    inv_share = rand(inverter_share_range)
    gfm_share = rand(gfm_share_range)
    
    # Algebraic 
    sys, p1, p3, perturbation1, perturbation3 = create_sys_with_params(file_name, kpv, kiv, kq, kpc, kic, load_scale, inv_share, gfm_share, gfm_bus, gfl_bus);
    stab, gfm_eta, gfl_eta, sm_eta, status = run_analysis(sys, false, false, p1, perturbation1,  gfm_bus, gfl_bus, sm_bus)
    push!(alg_df, [status, load_scale, kq, kpv, kiv, kpc, kic, stab, gfm_eta, gfl_eta, sm_eta, inv_share, gfm_share])

    # Dynamic 
    sys, p1, p3, perturbation1, perturbation3 = create_sys_with_params(file_name, kpv, kiv, kq, kpc, kic, load_scale, inv_share, gfm_share, gfm_bus, gfl_bus);
    stab, gfm_eta, gfl_eta, sm_eta, status = run_analysis(sys, true, false, p1, perturbation1,  gfm_bus, gfl_bus, sm_bus)
    push!(dyn_df, [status, load_scale, kq, kpv, kiv, kpc, kic, stab, gfm_eta, gfl_eta, sm_eta, inv_share, gfm_share])

    # # MSSB 
    # sys, p1, p3, perturbation1, perturbation3 = create_sys_with_params(file_name, kpv, kiv, kq, kpc, kic, load_scale, inv_share, gfm_share);
    # stab, gfm_eta, gfl_eta, sm_eta, status = run_analysis(sys, true, true, p1, perturbation1)
    # push!(mssb_df, [status, load_scale, kq, kpv, kiv, kpc, kic, stab, gfm_eta, gfl_eta, sm_eta, inv_share, gfm_share])

    # # MSMB 
    # sys, p1, p3, perturbation1, perturbation3 = create_sys_with_params(file_name, kpv, kiv, kq, kpc, kic, load_scale, inv_share, gfm_share);
    # stab, gfm_eta, gfl_eta, sm_eta, status = run_analysis(sys, true, true, p3, perturbation1)
    # push!(msmb_df, [status, load_scale, kq, kpv, kiv, kpc, kic, stab, gfm_eta, gfl_eta, sm_eta, inv_share, gfm_share])

end

function create_sys_with_params(file_name, kpv, kiv, kq, kpc, kic, load_scale, inv_share, gfm_share, gfm_bus, gfl_bus)


    sys = get_system(file_name)
    # Get generators 
    gfm = get_component(Generator, sys, gfm_bus)
    gfl = get_component(Generator, sys, gfl_bus)

    p1 = p1_9bus;
    p3 = p3_9bus;

    # # Update params 
    gfm_p = inv_share*gfm_share; # Create p setpoints 
    gfl_p = inv_share*(1.0-gfm_share); #

    # Adjust load scale 
    p1.load_scale = load_scale;
    p3.load_scale = load_scale;
    #Adjust GFM gains  
    gfm.dynamic_injector.outer_control.reactive_power_control.kq = kq
    gfm.dynamic_injector.inner_control.kpc = kpc
    gfm.dynamic_injector.inner_control.kpv = kpv
    gfm.dynamic_injector.inner_control.kic = kic
    gfm.dynamic_injector.inner_control.kiv = kiv

    #Adjust power injections 
    gfm.active_power = gfm_p * load_scale
    gfl.active_power = gfl_p * load_scale
    
    # Scale loads according to load scale 
    for l in get_components(PSY.StandardLoad, sys)
        transform_load_to_constant_impedance(l)
        l.impedance_active_power = l.impedance_active_power * p1.load_scale 
        l.impedance_reactive_power = l.impedance_reactive_power * p1.load_scale 
    end

    p1.perturbation_params.crc_params = CRCParam(DynamicGenerator, sm_bus, :V_ref, 0.95)
    p3.perturbation_params.crc_params = CRCParam(DynamicGenerator, sm_bus, :V_ref, 0.95)

    perturbation1 = choose_disturbance(sys, p1.perturbation, p1)
    perturbation3 = choose_disturbance(sys, p3.perturbation, p3)

    return sys, p1, p3, perturbation1, perturbation3

end



function remove_nans_from_df(df)
    filter(:stable => x -> !any(f -> f(x), (ismissing, isnan)), df)
end

function build_sim_from_file(file_name::String, dyn_lines::Bool, multi_segment::Bool, p::ExpParams, load_bus::String, inv_share, rating_scale, update_reference_bus::Bool)

    # build system
    sys = System(joinpath(pwd(), file_name));

    if update_reference_bus
        b1 = get_component(Bus, sys, "BUS 1")
        b1.bustype = BusTypes.PV

        b2 = get_component(Bus, sys, "BUS 2")
        b2.bustype = BusTypes.REF
    end

    # Simulation time span
    tspan = (0.0, p.sim_params.t_max);
    perturbation = ETL.choose_disturbance(sys, p.perturbation, p);

    # Add load
    load = StandardLoad(
        name = "load1",
        available = true,
        bus = get_component(Bus, sys, load_bus),
        base_power = 100.0,
        constant_active_power = 0.0,
        constant_reactive_power = 0.0,
        impedance_active_power = p.p_load*p.load_scale,
        impedance_reactive_power = p.q_load*p.load_scale,
        current_active_power = 0.0,
        current_reactive_power = 0.0,
        max_constant_active_power = 0.0,
        max_constant_reactive_power = 0.0,
        max_impedance_active_power = p.p_load*p.load_scale,
        max_impedance_reactive_power = p.q_load*p.load_scale,
        max_current_active_power = 0.0,
        max_current_reactive_power = 0.0,
    );
    add_component!(sys, load);

    for g in get_components(PSY.Generator, sys)
        set_rating!(g, g.rating * rating_scale) # adjust the ratings of all generators
        set_active_power!(g, g.rating * inv_share); # adjust the share of all generators 
        # if g.bus.name == inv_bus # at the inverter bus 
        #     set_active_power!(g, g.rating * inv_share); # adjust P setpoint to be correct fraction of rating 
        # end
    end
    
    # build segments model
    if (multi_segment == true)
        sys = ETL.build_seg_model!(sys, p, dyn_lines, "");
    else
        sys = ETL.build_new_impedance_model!(sys, p, dyn_lines, "");
    end

    sim = build_sim(sys, tspan, perturbation, dyn_lines, p);
    return sim 

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
