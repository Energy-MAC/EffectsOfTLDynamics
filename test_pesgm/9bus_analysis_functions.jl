using PowerSimulationsDynamics
using PowerSystems
using Parameters
using EffectsOfTLDynamics
using DataFrames 
using CSV 
using StatsPlots
const ETL = EffectsOfTLDynamics

include("../test/eig_analysis_funcs.jl")
include("generate_samples.jl")

function generate_single_case_2(samples_cut, inv_share, load_scale, gfl_share, line_scale, file_name, gfm_bus, gfl_bus, sm_bus, exp_folder, case)
    # Generates parameter samples, selects N samples, and builds dataframe of results for 2 different line models. Saves results and returns dfs.

    case_name = case*"_"*string(load_scale);
    
    alg, stab_alg, unst_alg, alg_nans = get_small_signal_results_from_samples(file_name, samples_cut, inv_share, load_scale, gfl_share, line_scale, gfm_bus, gfl_bus, sm_bus, false, false, false)
    
    dyn, stab_dyn, unst_dyn, dyn_nans = get_small_signal_results_from_samples(file_name, samples_cut, inv_share, load_scale, gfl_share, line_scale, gfm_bus, gfl_bus, sm_bus, true, false, false)

    # SAVE DATA 
    CSV.write(exp_folder*"/alg_"*case_name*".csv", alg)
    CSV.write(exp_folder*"/dyn_"*case_name*".csv", dyn)

    return alg, dyn, case_name 
end

function generate_single_case_3(samples_cut, inv_share, load_scale, gfl_share, line_scale, file_name, gfm_bus, gfl_bus, sm_bus, exp_folder, case)
    # Takes a sample of parameters, and builds dataframe of results for 3 different line models. Saves results and returns dfs.

    case_name = case*"_"*string(load_scale);
    
    alg = get_small_signal_results_from_samples(file_name, samples_cut, inv_share, load_scale, gfl_share, line_scale, gfm_bus, gfl_bus, sm_bus, false, false, false)
    
    dyn = get_small_signal_results_from_samples(file_name, samples_cut, inv_share, load_scale, gfl_share, line_scale, gfm_bus, gfl_bus, sm_bus, true, false, false)
    
    mssb = get_small_signal_results_from_samples(file_name, samples_cut, inv_share, load_scale, gfl_share, line_scale, gfm_bus, gfl_bus, sm_bus, true, true, false)

    # SAVE DATA 
    CSV.write(exp_folder*"/alg_"*case_name*".csv", alg)
    CSV.write(exp_folder*"/dyn_"*case_name*".csv", dyn)
    CSV.write(exp_folder*"/mssb_"*case_name*".csv", mssb)

    #return alg, dyn, mssb, case_name
end

function generate_single_case_4(samples_cut, inv_share, load_scale, gfl_share, line_scale, file_name, gfm_bus, gfl_bus, sm_bus, exp_folder, case)
    # Takes a sample of parameters, and builds dataframe of results for 3 different line models. Saves results and returns dfs.

    case_name = case*"_"*string(load_scale);
    
    alg, stab_alg, unst_alg, alg_nans = get_small_signal_results_from_samples(file_name, samples_cut, inv_share, load_scale, gfl_share, line_scale, gfm_bus, gfl_bus, sm_bus, false, false, false)
    
    dyn, stab_dyn, unst_dyn, dyn_nans = get_small_signal_results_from_samples(file_name, samples_cut, inv_share, load_scale, gfl_share, line_scale, gfm_bus, gfl_bus, sm_bus, true, false, false)
    
    mssb, stab_mssb, unst_mssb, mssb_nans = get_small_signal_results_from_samples(file_name, samples_cut, inv_share, load_scale, gfl_share, line_scale, gfm_bus, gfl_bus, sm_bus, true, true, false)

    msmb, stab_msmb, unst_msmb, msmb_nans = get_small_signal_results_from_samples(file_name, samples_cut, inv_share, load_scale, gfl_share, line_scale, gfm_bus, gfl_bus, sm_bus, true, true, true)

    # SAVE DATA 
    CSV.write(exp_folder*"/alg_"*case_name*".csv", alg)
    CSV.write(exp_folder*"/dyn_"*case_name*".csv", dyn)
    CSV.write(exp_folder*"/mssb_"*case_name*".csv", mssb)
    CSV.write(exp_folder*"/msmb_"*case_name*".csv", msmb)

    return alg, dyn, mssb, msmb, case_name
end

function make_case_boxplot_2(alg, dyn, legend_loc, case_name, exp_folder)
    # No MSSB 
    w = 0.9;
    labels = ["statpi", "dynpi"];

    @df alg boxplot(ones(size(alg)[1]), :stable, bar_width=w, label=labels[1], legend=legend_loc, ylabel=L"\Re \ (\lambda)", xticks=[])
    @df dyn boxplot!(2*ones(size(alg)[1]), :stable, bar_width=w, label=labels[2])
    #plot!([0.5,2.5], [0,0], linestyle=:dash, color=:black, primary=false)

    # labelfontsize=25;
    # ticksize = 16;
    # legendsize = 16;
    # plot!(xguidefontsize=labelfontsize, yguidefontsize=labelfontsize, tickfontsize=ticksize,legendfontsize=legendsize, legendtitlesize=legendsize)
    
    savefig(exp_folder*"/boxplot_"*case_name*".png")
    savefig(exp_folder*"/boxplot_"*case_name*".svg")

end

function height_adj(bp, target_height)
    yl = ylims(bp)
    mid = (yl[2]+yl[1])/2
    new_max = mid+target_height/2
    new_min = mid-target_height/2
    ylims!(new_min,new_max)
end

function make_case_boxplot_3(alg, dyn, mssb, legend_toggle, case_name, exp_folder, yaxis_height)
    w = 0.9;
    labels = [L"\mathrm{statpi}", L"\mathrm{dynpi}", L"\mathrm{MSSB}"];

    labelfontsize=16;
    ticksize = 14;
    legendsize = 16;
    plot_font = "Computer Modern"
    default(fontfamily=plot_font,xguidefontsize=labelfontsize, yguidefontsize=labelfontsize, tickfontsize=ticksize,legendfontsize=legendsize,framestyle=:box, widen=true)
    
    bp = @df alg boxplot(ones(size(alg)[1]), :stable, bar_width=w, label=labels[1], legend=legend_toggle, ylabel=L"\Re \ (\lambda)", xticks=[])
    @df dyn boxplot!(2*ones(size(alg)[1]), :stable, bar_width=w, label=labels[2])
    @df mssb boxplot!(3*ones(size(alg)[1]), :stable, bar_width=w, label=labels[3])
    height_adj(bp, yaxis_height)

    savefig(exp_folder*"/boxplot_"*case_name*".png")
    #savefig(exp_folder*"/boxplot_"*case_name*".svg")

end

function make_case_boxplot_4(alg, dyn, mssb, msmb, legend_loc, case_name, exp_folder)
    w = 0.9;
    labels = ["statpi", "dynpi", "mssb", "msmb"];
    
    @df alg boxplot(ones(size(alg)[1]), :stable, bar_width=w, label=labels[1], legend=legend_loc, ylabel=L"\Re \ (\lambda)", xticks=[], )
    @df dyn boxplot!(2*ones(size(alg)[1]), :stable, bar_width=w, label=labels[2])
    @df mssb boxplot!(3*ones(size(alg)[1]), :stable, bar_width=w, label=labels[3])
    @df msmb boxplot!(4*ones(size(alg)[1]), :stable, bar_width=w, label=labels[4])

    savefig(exp_folder*"/boxplot_"*case_name*".png")
    savefig(exp_folder*"/boxplot_"*case_name*".svg")

end

function make_gain_boxplots(sel, labels, legend_key, case_name, exp_folder)
    # Produces boxplots for a fixed operating condition (fixed load scale, gfm share, gfl share)
    if exp_folder[end] != "/";
        exp_folder = exp_folder*"/"
    end

    num_line_models = length(sel);
    if num_line_models == 2;
        get_o = get_2o;
    elseif num_line_models == 3;
        get_o = get_3o;
    else
        get_o = get_4o;
    end

    labelfontsize=13;
    ticksize = 11;
    legendsize = 10;
    plot_font = "Computer Modern"
    default(fontfamily=plot_font,xguidefontsize=labelfontsize, yguidefontsize=labelfontsize, tickfontsize=ticksize,legendfontsize=legendsize,framestyle=:box, widen=true)

    labelfontsize=13;
    ticksize = 11;
    legendsize = 10;

    # KQ 
    w = 0.008
    plot()
    for i = 1:num_line_models;
        @df sel[i] boxplot!(:kq.+get_o(w,i), :stable, bar_width=w, label=labels[i],  xlabel=L"k_q", ylabel=L"\Re \ (\lambda)", legend=false)
    end
    savefig(exp_folder*"kq_case_"*case_name*".png")
    savefig(exp_folder*"kq_case_"*case_name*".svg")

    # KPV 
    w = 0.004
    plot()
    for i = 1:num_line_models;
        @df sel[i] boxplot!(:kpv.+get_o(w,i), :stable, bar_width=w, label=labels[i], xlabel=L"k_{pv}", ylabel=L"\Re \ (\lambda)", legend=false)
    end

    savefig(exp_folder*"kpv_case_"*case_name*".png")
    #savefig(exp_folder*"kpv_case_"*case_name*".svg")

    # KIV 
    w = 20.0
    plot()
    for i = 1:num_line_models;
        @df sel[i] boxplot!(:kiv.+get_o(w,i), :stable, bar_width=w, label=labels[i], xlabel=L"k_{iv}", ylabel=L"\Re \ (\lambda)", legend=false)
    end
    savefig(exp_folder*"kiv_case_"*case_name*".png")
    #savefig(exp_folder*"kiv_case_"*case_name*".svg")

    # KPC 
    w = 0.03;
    plot()
    for i = 1:num_line_models;
        @df sel[i] boxplot!(:kpc.+get_o(w,i), :stable, bar_width=w, label=labels[i], xlabel=L"k_{pc}", ylabel=L"\Re \ (\lambda)", legend=false)
    end
    savefig(exp_folder*"kpc_case_"*case_name*".png")
    #savefig(exp_folder*"kpc_case_"*case_name*".svg")

    # KIC 
    w = 0.2;
    plot()
    for i = 1:num_line_models;
        @df sel[i] boxplot!(:kic.+get_o(w,i), :stable, bar_width=w, label=labels[i], xlabel=L"k_{ic}", ylabel=L"\Re \ (\lambda)", legend=false)
    end
    savefig(exp_folder*"kic_case_"*case_name*".png")
    #savefig(exp_folder*"kic_case_"*case_name*".svg")

    # kd
    w = 15
    plot()
    for i = 1:num_line_models;
        @df sel[i] boxplot!(:kd.+get_o(w,i), :stable, bar_width=w, label=labels[i], legend=legend_key, xlabel=L"k_{d}", ylabel=L"\Re \ (\lambda)")
    end
    savefig(exp_folder*"kd_case_"*case_name*".png")
    #savefig(exp_folder*"kd_case_"*case_name*".svg")

    # Ta
    w = 0.07
    plot()
    for i = 1:num_line_models;
        @df sel[i] boxplot!(:Ta.+get_o(w,i), :stable, bar_width=w, label=labels[i], xlabel=L"T_a", ylabel=L"\Re \ (\lambda)", legend=false)
    end
    savefig(exp_folder*"Ta_case_"*case_name*".png")
    #savefig(exp_folder*"Ta_case_"*case_name*".svg")

end

function get_2o(w,n)
    offset2 = [-0.5*w, 0.5*w]
    return offset2[n]
end

function get_3o(w,n)
    offset3 = [-w, 0, w]
    return offset3[n]
end

function get_4o(w,n)
    offset4 = [w*1.5, w/2, -w/2, -w*1.5]
    return offset4[n]
end

function load_min_from_folder(save_folder)
    alg = DataFrame(CSV.File(save_folder*"/alg.csv"))
    dyn = DataFrame(CSV.File(save_folder*"/dyn.csv"))

    return alg, dyn 
end

function load_full_from_folder(save_folder)
    alg = DataFrame(CSV.File(save_folder*"/alg.csv"))
    dyn = DataFrame(CSV.File(save_folder*"/dyn.csv"))
    mssb = DataFrame(CSV.File(save_folder*"/mssb.csv"))
    return alg, dyn, mssb 
end

function get_small_signal_results_from_samples(file_name, df, inv_share, load_scale, gfl_share, line_scale, gfm_bus, gfl_bus, sm_bus, dyn_line, multiseg, fd)

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
        gfm_share = 1.0 - gfl_share;

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
    # stable_df = filter(:stable => x -> x<0, df)
    # unstable_df = filter(:stable => x -> x>0, df)
    
    return df

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

gfm_bus = "generator-3-1";
gfl_bus = "generator-2-1";
sm_bus = "generator-1-1";