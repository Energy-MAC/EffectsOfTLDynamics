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

file_name = "../data/json_data/9bus_VSM3_SM1_GFL2.json"

# defaults 
gfm_bus = "generator-3-1";
gfl_bus = "generator-2-1";
sm_bus = "generator-1-1";

N = 50;

save_folder = "Nov3_extra_loading_share_combinations/";
labels = ["statpi", "dynpi", "mssb"]
w = 0.9
legend_loc = :left

samples = generate_single_case_params(0.4, 0.5, 0.2, 3.0); # inv_share, load_scale, gfl_share, line_scale 
sample_name = "low_load_low_gfl"
samples_shuffled = shuffle(MersenneTwister(123), samples);
samples_cut = samples_shuffled[1:N,:]

alg, stab_alg, unst_alg, alg_nans = get_small_signal_results_from_samples(file_name, samples_cut, gfm_bus, gfl_bus, sm_bus, false, false, false)

dyn, stab_dyn, unst_dyn, dyn_nans = get_small_signal_results_from_samples(file_name, samples_cut, gfm_bus, gfl_bus, sm_bus, true, false, false)

mssb, stab_mssb, unst_mssb, mssb_nans = get_small_signal_results_from_samples(file_name, samples_cut, gfm_bus, gfl_bus, sm_bus, true, true, false)
 
@df alg boxplot(ones(size(alg)[1]), :stable, bar_width=w, label=labels[1], legend=legend_loc, ylabel=L"R(\lambda)", xticks=[])
@df dyn boxplot!(2*ones(size(alg)[1]), :stable, bar_width=w, label=labels[2])
@df mssb boxplot!(3*ones(size(alg)[1]), :stable, bar_width=w, label=labels[3])

plot!([0.5,3.5], [0,0], linestyle=:dash, color=:black, primary=false)

savefig(save_folder*sample_name*"_boxplot.png")
savefig(save_folder*sample_name*"_boxplot.svg")

sample_name = "nom_load_low_gfl";
samples = generate_single_case_params(0.4, 1.0, 0.2, 3.0)
samples_shuffled = shuffle(MersenneTwister(123), samples);
samples_cut = samples_shuffled[1:N,:]

alg, stab_alg, unst_alg, alg_nans = get_small_signal_results_from_samples(file_name, samples_cut, gfm_bus, gfl_bus, sm_bus, false, false, false)

dyn, stab_dyn, unst_dyn, dyn_nans = get_small_signal_results_from_samples(file_name, samples_cut, gfm_bus, gfl_bus, sm_bus, true, false, false)

mssb, stab_mssb, unst_mssb, mssb_nans = get_small_signal_results_from_samples(file_name, samples_cut, gfm_bus, gfl_bus, sm_bus, true, true, false)

@df alg boxplot(ones(size(alg)[1]), :stable, bar_width=w, label=labels[1], legend=legend_loc, ylabel=L"R(\lambda)", xticks=[])
@df dyn boxplot!(2*ones(size(alg)[1]), :stable, bar_width=w, label=labels[2])
@df mssb boxplot!(3*ones(size(alg)[1]), :stable, bar_width=w, label=labels[3])

plot!([0.5,3.5], [0,0], linestyle=:dash, color=:black, primary=false)

savefig(save_folder*sample_name*"_boxplot.png")
savefig(save_folder*sample_name*"_boxplot.svg")

sample_name = "low_load_high_gfl";
samples = generate_single_case_params(0.4, 0.5, 0.8, 3.0)
samples_shuffled = shuffle(MersenneTwister(123), samples);
samples_cut = samples_shuffled[1:N,:]

alg, stab_alg, unst_alg, alg_nans = get_small_signal_results_from_samples(file_name, samples_cut, gfm_bus, gfl_bus, sm_bus, false, false, false)

dyn, stab_dyn, unst_dyn, dyn_nans = get_small_signal_results_from_samples(file_name, samples_cut, gfm_bus, gfl_bus, sm_bus, true, false, false)

mssb, stab_mssb, unst_mssb, mssb_nans = get_small_signal_results_from_samples(file_name, samples_cut, gfm_bus, gfl_bus, sm_bus, true, true, false)

@df alg boxplot(ones(size(alg)[1]), :stable, bar_width=w, label=labels[1], legend=legend_loc, ylabel=L"R(\lambda)", xticks=[])
@df dyn boxplot!(2*ones(size(alg)[1]), :stable, bar_width=w, label=labels[2])
@df mssb boxplot!(3*ones(size(alg)[1]), :stable, bar_width=w, label=labels[3])

plot!([0.5,3.5], [0,0], linestyle=:dash, color=:black, primary=false)

savefig(save_folder*sample_name*"_boxplot.png")
savefig(save_folder*sample_name*"_boxplot.svg")


sample_name = "nom_load_high_gfl";
samples = generate_single_case_params(0.4, 1.0, 0.8, 3.0)
samples_shuffled = shuffle(MersenneTwister(123), samples);
samples_cut = samples_shuffled[1:N,:]

alg, stab_alg, unst_alg, alg_nans = get_small_signal_results_from_samples(file_name, samples_cut, gfm_bus, gfl_bus, sm_bus, false, false, false)

dyn, stab_dyn, unst_dyn, dyn_nans = get_small_signal_results_from_samples(file_name, samples_cut, gfm_bus, gfl_bus, sm_bus, true, false, false)

mssb, stab_mssb, unst_mssb, mssb_nans = get_small_signal_results_from_samples(file_name, samples_cut, gfm_bus, gfl_bus, sm_bus, true, true, false)

@df alg boxplot(ones(size(alg)[1]), :stable, bar_width=w, label=labels[1], legend=legend_loc, ylabel=L"R(\lambda)", xticks=[])
@df dyn boxplot!(2*ones(size(alg)[1]), :stable, bar_width=w, label=labels[2])
@df mssb boxplot!(3*ones(size(alg)[1]), :stable, bar_width=w, label=labels[3])
plot!([0.5,3.5], [0,0], linestyle=:dash, color=:black, primary=false)

savefig(save_folder*sample_name*"_boxplot.png")
savefig(save_folder*sample_name*"_boxplot.svg")


function make_boxplot(N, inv_share, load_scale, gfl_share, line_scale, file_name, gfm_bus, gfl_bus, sm_bus, sample_name, save_folder)

    w = 0.9;
    labels = ["statpi", "dynpi", "mssb"];

    samples = generate_single_case_params(inv_share, load_scale, gfl_share, line_scale)
    samples_shuffled = shuffle(MersenneTwister(123), samples);
    samples_cut = samples_shuffled[1:N,:]
    
    alg, stab_alg, unst_alg, alg_nans = get_small_signal_results_from_samples(file_name, samples_cut, gfm_bus, gfl_bus, sm_bus, false, false, false)
    
    dyn, stab_dyn, unst_dyn, dyn_nans = get_small_signal_results_from_samples(file_name, samples_cut, gfm_bus, gfl_bus, sm_bus, true, false, false)
    
    mssb, stab_mssb, unst_mssb, mssb_nans = get_small_signal_results_from_samples(file_name, samples_cut, gfm_bus, gfl_bus, sm_bus, true, true, false)
    
    @df alg boxplot(ones(size(alg)[1]), :stable, bar_width=w, label=labels[1], legend=legend_loc, ylabel=L"R(\lambda)", xticks=[])
    @df dyn boxplot!(2*ones(size(alg)[1]), :stable, bar_width=w, label=labels[2])
    @df mssb boxplot!(3*ones(size(alg)[1]), :stable, bar_width=w, label=labels[3])
    plot!([0.5,3.5], [0,0], linestyle=:dash, color=:black, primary=false)
    
    savefig(save_folder*sample_name*"_boxplot.png")
    savefig(save_folder*sample_name*"_boxplot.svg")

end


function make_minimum_boxplot(N, inv_share, load_scale, gfl_share, line_scale, file_name, gfm_bus, gfl_bus, sm_bus, sample_name, save_folder)
    # No MSSB 

    w = 0.9;
    labels = ["statpi", "dynpi", "mssb"];

    samples = generate_single_case_params(inv_share, load_scale, gfl_share, line_scale)
    samples_shuffled = shuffle(MersenneTwister(123), samples);
    samples_cut = samples_shuffled[1:N,:]
    
    alg, stab_alg, unst_alg, alg_nans = get_small_signal_results_from_samples(file_name, samples_cut, gfm_bus, gfl_bus, sm_bus, false, false, false)
    
    dyn, stab_dyn, unst_dyn, dyn_nans = get_small_signal_results_from_samples(file_name, samples_cut, gfm_bus, gfl_bus, sm_bus, true, false, false)
    
    @df alg boxplot(ones(size(alg)[1]), :stable, bar_width=w, label=labels[1], legend=legend_loc, ylabel=L"R(\lambda)", xticks=[])
    @df dyn boxplot!(2*ones(size(alg)[1]), :stable, bar_width=w, label=labels[2])
    plot!([0.5,3.5], [0,0], linestyle=:dash, color=:black, primary=false)
    
    savefig(save_folder*sample_name*"_boxplot.png")
    savefig(save_folder*sample_name*"_boxplot.svg")

end

save_folder = "Nov3_boxplots/"
for load_scale = 0.5:0.1:1.0;
    make_minimum_boxplot(50, 0.6, load_scale, 0.4, 1.0, file_name, gfm_bus, gfl_bus, sm_bus, "test"*string(load_scale), save_folder)
end
