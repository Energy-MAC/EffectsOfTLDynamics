cd(@__DIR__)
using PowerSystems
using PowerSimulationsDynamics
using Plots
using Revise
using EffectsOfTLDynamics
using LaTeXStrings
using DataFrames
using StatsPlots
using Random 
const ETL = EffectsOfTLDynamics;
const PSY = PowerSystems;
const PSID = PowerSimulationsDynamics;

include("9bus_defaults.jl")
include("../test/eig_analysis_funcs.jl")
include("generate_samples.jl")

function generate_single_case_full(N, inv_share, load_scale, gfl_share, line_scale, file_name, gfm_bus, gfl_bus, sm_bus)
    # Generates parameter samples, selects N samples, and builds dataframe of results for 3 different line models. Saves results and returns dfs.

    ls_str = string.([N,inv_share,load_scale,gfl_share,line_scale])
    case_name = "case"*join(ls_str, "_")

    samples = generate_single_case_params(inv_share, load_scale, gfl_share, line_scale)
    samples_shuffled = shuffle(MersenneTwister(123), samples);
    samples_cut = samples_shuffled[1:N,:]
    
    alg, stab_alg, unst_alg, alg_nans = get_small_signal_results_from_samples(file_name, samples_cut, gfm_bus, gfl_bus, sm_bus, false, false, false)
    
    dyn, stab_dyn, unst_dyn, dyn_nans = get_small_signal_results_from_samples(file_name, samples_cut, gfm_bus, gfl_bus, sm_bus, true, false, false)
    
    mssb, stab_mssb, unst_mssb, mssb_nans = get_small_signal_results_from_samples(file_name, samples_cut, gfm_bus, gfl_bus, sm_bus, true, true, false)

    # SAVE DATA 
    mkdir(case_name)
    CSV.write(case_name*"/alg.csv", alg)
    CSV.write(case_name*"/dyn.csv", dyn)
    CSV.write(case_name*"/mssb.csv", dyn)

    return alg, dyn, mssb, case_name

end

function generate_single_case_min(N, inv_share, load_scale, gfl_share, line_scale, file_name, gfm_bus, gfl_bus, sm_bus)
    # Generates parameter samples, selects N samples, and builds dataframe of results for 2 different line models. Saves results and returns dfs.
    ls_str = string.([N,inv_share,load_scale,gfl_share,line_scale])
    case_name = "case_"*join(ls_str, "_")

    samples = generate_single_case_params(inv_share, load_scale, gfl_share, line_scale)
    samples_shuffled = shuffle(MersenneTwister(123), samples);
    samples_cut = samples_shuffled[1:N,:]
    
    alg, stab_alg, unst_alg, alg_nans = get_small_signal_results_from_samples(file_name, samples_cut, gfm_bus, gfl_bus, sm_bus, false, false, false)
    
    dyn, stab_dyn, unst_dyn, dyn_nans = get_small_signal_results_from_samples(file_name, samples_cut, gfm_bus, gfl_bus, sm_bus, true, false, false)

    # SAVE DATA 
    mkdir(case_name)
    CSV.write(case_name*"/alg.csv", alg)
    CSV.write(case_name*"/dyn.csv", dyn)

    return alg, dyn, case_name
end


function make_case_boxplot_min(alg, dyn, legend_loc, save_folder)
    # No MSSB 
    w = 0.9;
    labels = ["statpi", "dynpi"];
    
    @df alg boxplot(ones(size(alg)[1]), :stable, bar_width=w, label=labels[1], legend=legend_loc, ylabel=L"R(\lambda)", xticks=[])
    @df dyn boxplot!(2*ones(size(alg)[1]), :stable, bar_width=w, label=labels[2])
    plot!([0.5,2.5], [0,0], linestyle=:dash, color=:black, primary=false)
    
    savefig(save_folder*"/boxplot.png")
    savefig(save_folder*"/boxplot.svg")

end

function make_case_boxplot_full(alg, dyn, mssb, legend_loc, save_folder)
    # No MSSB 
    w = 0.9;
    labels = ["statpi", "dynpi", "mssb"];
    
    @df alg boxplot(ones(size(alg)[1]), :stable, bar_width=w, label=labels[1], legend=legend_loc, ylabel=L"R(\lambda)", xticks=[])
    @df dyn boxplot!(2*ones(size(alg)[1]), :stable, bar_width=w, label=labels[2])
    @df mssb boxplot!(3*ones(size(alg)[1]), :stable, bar_width=w, label=labels[3])
    plot!([0.5,3.5], [0,0], linestyle=:dash, color=:black, primary=false)
    
    
    savefig(save_folder*"/boxplot.png")
    savefig(save_folder*"/boxplot.svg")

end


function make_gain_boxplots(sel, labels, legend_loc, save_folder)
    # Produces boxplots for a fixed operating condition (fixed load scale, gfm share, gfl share)
    if save_folder[end] != "/";
        save_folder = save_folder*"/"
    end

    num_line_models = length(sel);
    if num_line_models == 2;
        get_o = get_2o;
    elseif num_line_models == 3;
        get_o = get_3o;
    else
        get_o = get_4o;
    end

    # KQ 
    w = 0.005
    plot()
    for i = 1:num_line_models;
        @df sel[i] boxplot!(:kq.+get_o(w,i), :stable, bar_width=w, label=labels[i], legend=legend_loc, xlabel=L"k_q", ylabel=L"R(\lambda)")
    end
    savefig(save_folder*"kq_case.png")
    savefig(save_folder*"kq_case.svg")

    # KPV 
    w = 0.002
    plot()
    for i = 1:num_line_models;
        @df sel[i] boxplot!(:kpv.+get_o(w,i), :stable, bar_width=w, label=labels[i], legend=legend_loc, xlabel=L"k_{pv}", ylabel=L"R(\lambda)")
    end

    savefig(save_folder*"kpv_case.png")
    savefig(save_folder*"kpv_case.svg")

    # KIV 
    w = 15.0
    plot()
    for i = 1:num_line_models;
        @df sel[i] boxplot!(:kiv.+get_o(w,i), :stable, bar_width=w, label=labels[i], legend=legend_loc, xlabel=L"k_{iv}", ylabel=L"R(\lambda)")
    end

    savefig(save_folder*"kiv_case.png")
    savefig(save_folder*"kiv_case.svg")

    # KPC 
    w = 0.01;
    plot()
    for i = 1:num_line_models;
        @df sel[i] boxplot!(:kpc.+get_o(w,i), :stable, bar_width=w, label=labels[i], legend=legend_loc, xlabel=L"k_{pc}", ylabel=L"R(\lambda)")
    end
    
    savefig(save_folder*"kpc_case.png")
    savefig(save_folder*"kpc_case.svg")

    # KIC 
    w = 0.2;
    plot()
    for i = 1:num_line_models;
        @df sel[i] boxplot!(:kic.+get_o(w,i), :stable, bar_width=w, label=labels[i], legend=legend_loc, xlabel=L"k_{ic}", ylabel=L"R(\lambda)")
    end

    savefig(save_folder*"kic_case.png")
    savefig(save_folder*"kic_case.svg")

    # kd
    w = 10
    plot()
    for i = 1:num_line_models;
        @df sel[i] boxplot!(:kd.+get_o(w,i), :stable, bar_width=w, label=labels[i], legend=legend_loc, xlabel=L"k_{d}", ylabel=L"R(\lambda)")
    end

    savefig(save_folder*"kd_case.png")
    savefig(save_folder*"kd_case.svg")

    # Ta
    w = 0.05
    plot()
    for i = 1:num_line_models;
        @df sel[i] boxplot!(:Ta.+get_o(w,i), :stable, bar_width=w, label=labels[i], legend=legend_loc, xlabel=L"T_a", ylabel=L"R(\lambda)")
    end

    savefig(save_folder*"Ta_case.png")
    savefig(save_folder*"Ta_case.svg")

end

function get_2o(w,n)
    offset2 = [-0.5*w, 0.5*w]
    return offset2[n]
end

function get_3o(w,n)
    offset3 = [w, 0, -w]
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