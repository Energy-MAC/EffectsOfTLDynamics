cd(@__DIR__)
using PowerSystems
using PowerSimulationsDynamics
using Plots
using Revise
using EffectsOfTLDynamics
using LaTeXStrings
using DataFrames
using StatsPlots

const ETL = EffectsOfTLDynamics;
const PSY = PowerSystems;
const PSID = PowerSimulationsDynamics;

include("defaults.jl")
include("../test/eig_analysis_funcs.jl")

file_name = "../data/json_data/9bus_VSM_SM_GFL_.json"
#file_name = "../data/json_data/9bus_slackless.json"

# define buses 
gfm_bus = "generator-1-1";
gfl_bus = "generator-3-1";
sm_bus = "generator-2-1";

N = 50; # number of samples 

# results DF
col_names = ["STATUS", "load_scale", "kq", "kpv", "kiv", "kpc", "kic", "stable", "gfm_eta", "gfl_eta", "sm_eta", "sum_eta", "inv_share", "gfm_share"]
df = DataFrame([name => [] for name in col_names])
df_dyn = DataFrame([name => [] for name in col_names])

# Set up 
sys = get_system(file_name)

p = p1_9bus;
tspan = (0.0, p.sim_params.t_max)
perturbation = choose_disturbance(sys, p.perturbation, p)
alg_line_name = p.perturbation_params.branch_trip_params.line_to_trip

for n = 1:N;
    # Sample GFM params 
    kpv = rand(kpv_markovic_range); # inner loop voltage gain 
    kiv = rand(kiv_markovic_range); # inner loop integral gain 
    kq = rand(gfm_kq_range); # outer loop reactive power control gain
    
    kpc = rand(kpc_range); # Inner current loop gain - GFM 
    kic = rand(kic_range); # Inner current loop gain - GFM

    # Sample system params 
    load_scale = rand(load_scale_range); # Load scaling 
    inv_share = rand(inverter_share_range); # share of generation from IBRs
    gfm_share = rand(gfm_share_range); # share from GFM
    gfm_p = inv_share*gfm_share*load_scale; # Create p setpoints 
    gfl_p = inv_share*(1.0-gfm_share)*load_scale;

    # Create system 
    sys = get_system(file_name) # to be safe, rebuild these each time 
    gfm = get_component(Generator, sys, gfm_bus)
    gfm.base_power = 100;
    gfl = get_component(Generator, sys, gfl_bus)
    gfl.base_power = 100;

    p = p1_9bus;

    # Update params 
    p.load_scale = load_scale;
    # GFM 
    gfm.dynamic_injector.outer_control.reactive_power_control.kq = kq
    gfm.dynamic_injector.inner_control.kpc = kpc
    gfm.dynamic_injector.inner_control.kpv = kpv
    gfm.dynamic_injector.inner_control.kic = kic
    gfm.dynamic_injector.inner_control.kiv = kiv
    gfm.active_power = gfm_p 
    # GFL 
    gfl.dynamic_injector.inner_control.kpc = kpc
    gfl.dynamic_injector.inner_control.kic = kic 
    gfl.active_power = gfl_p 

    # Scale loads according to load scale 
    for l in get_components(PSY.StandardLoad, sys)
        transform_load_to_constant_impedance(l)
        l.impedance_active_power = l.impedance_active_power * p.load_scale 
        l.impedance_reactive_power = l.impedance_reactive_power * p.load_scale 
    end

    sys_alg = ETL.build_new_impedance_model!(sys, p, false, alg_line_name)
    sim_alg = build_sim(sys_alg, tspan, perturbation, false, p);
    gfm_eta, gfl_eta, sm_eta = calculate_etas(sim_alg.sys, gfm_bus, gfl_bus, sm_bus)
    stab_alg = save_max_nonzero_eig!(sim_alg, [])[1]

    # Save parameter values and results
    push!(df, [sim_alg.status, load_scale, kq, kpv, kiv, kpc, kic, stab_alg, gfm_eta, gfl_eta, sm_eta, gfm_eta+gfl_eta+sm_eta, inv_share, gfm_share])

    # Dynamic line 
    sys_dyn = get_system(file_name) # to be safe, rebuild these each time 
    gfm = get_component(Generator, sys_dyn, gfm_bus)
    gfm.base_power = 100;
    gfl = get_component(Generator, sys_dyn, gfl_bus)
    gfl.base_power = 100;

    p = p1_9bus;
    p.load_scale = load_scale;
    # GFM 
    gfm.dynamic_injector.outer_control.reactive_power_control.kq = kq
    gfm.dynamic_injector.inner_control.kpc = kpc
    gfm.dynamic_injector.inner_control.kpv = kpv
    gfm.dynamic_injector.inner_control.kic = kic
    gfm.dynamic_injector.inner_control.kiv = kiv
    gfm.active_power = gfm_p # ??? check 
    # GFL  
    gfl.dynamic_injector.inner_control.kpc = kpc
    gfl.dynamic_injector.inner_control.kic = kic 
    gfl.active_power = gfl_p  # should this also be scaled by the load factor ?

    # Scale loads according to load scale 
    for l in get_components(PSY.StandardLoad, sys_dyn)
        transform_load_to_constant_impedance(l)
        l.impedance_active_power = l.impedance_active_power * p.load_scale 
        l.impedance_reactive_power = l.impedance_reactive_power * p.load_scale 
    end

    sys_dyn = ETL.build_new_impedance_model!(sys_dyn, p, true, alg_line_name)
    sim_dyn = build_sim(sys_dyn, tspan, perturbation, true, p);
    gfm_eta, gfl_eta, sm_eta = calculate_etas(sim_dyn.sys, gfm_bus, gfl_bus, sm_bus)

    stab_dyn = save_max_nonzero_eig!(sim_dyn, [])[1]
    push!(df_dyn, [sim_dyn.status, load_scale, kq, kpv, kiv, kpc, kic, stab_dyn, gfm_eta, gfl_eta, sm_eta,gfm_eta+gfl_eta+sm_eta, inv_share, gfm_share])

end

# df
# df_dyn 

# Remove nans from dfs 
df = filter(:stable => x -> !any(f -> f(x), (ismissing, isnan)), df)
df_dyn = filter(:stable => x -> !any(f -> f(x), (ismissing, isnan)), df_dyn)


# KQ 
w = 0.01;
offset = w/2;
@df df boxplot(:kq.-offset, :stable, bar_width=w, label="Algebraic", legend=:bottom, xlabel="Kq", ylabel=L"R(\lambda)")
@df df_dyn boxplot!(:kq.+offset, :stable, bar_width=w, label="Dyn", legend=:bottom)

# LOADING
w = 0.05;
offset = w/2;
@df df boxplot(:load_scale.-offset, :stable, bar_width=w, label="Algebraic", legend=:bottom, xlabel="Load scale", ylabel=L"R(\lambda)")
@df df_dyn boxplot!(:load_scale.+offset, :stable, bar_width=w, label="Dyn", legend=:bottom)

# GFM SHARE 
w = 0.01;
offset = w/2;
@df df boxplot(:gfm_eta.-offset, :stable, bar_width=w, label="Algebraic", legend=:bottom, xlabel="gfm share", ylabel=L"R(\lambda)")
@df df_dyn boxplot!(:gfm_eta.+offset, :stable, bar_width=w, label="Dyn", legend=:bottom)

# GFL SHARE 
w = 0.01;
offset = w/2;
@df df boxplot(:gfl_eta.-offset, :stable, bar_width=w, label="Algebraic", legend=:bottom, xlabel="gfl share", ylabel=L"R(\lambda)")
@df df_dyn boxplot!(:gfl_eta.+offset, :stable, bar_width=w, label="Dyn", legend=:bottom)

# KPV 
w = 0.0005;
offset = w/2;
@df df boxplot(:kpv.-offset, :stable, bar_width=w, label="Algebraic", legend=:bottom, xlabel="kpv (inner loop GFM)", ylabel=L"R(\lambda)")
@df df_dyn boxplot!(:kpv.+offset, :stable, bar_width=w, label="Dyn", legend=:bottom)

#KIV 
w = 0.2;
offset = w/2;
@df df boxplot(:kiv.-offset, :stable, bar_width=w, label="Algebraic", legend=:bottom, xlabel="kiv (inner loop GFM)", ylabel=L"R(\lambda)")
@df df_dyn boxplot!(:kiv.+offset, :stable, bar_width=w, label="Dyn", legend=:bottom)

# KPC 
w = 0.01;
offset = w/2;
@df df boxplot(:kpc.-offset, :stable, bar_width=w, label="Algebraic", legend=:bottom, xlabel="kpc (inner loop GFM,GFL)", ylabel=L"R(\lambda)")
@df df_dyn boxplot!(:kpc.+offset, :stable, bar_width=w, label="Dyn", legend=:bottom)

# KIC 
w = 0.3;
offset = w/2;
@df df boxplot(:kic.-offset, :stable, bar_width=w, label="Algebraic", legend=:bottom, xlabel="kic (inner loop GFM,GFL)", ylabel=L"R(\lambda)")
@df df_dyn boxplot!(:kic.+offset, :stable, bar_width=w, label="Dyn", legend=:bottom)
