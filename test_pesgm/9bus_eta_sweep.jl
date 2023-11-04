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

N = 100;

save_folder = "Nov2_conditioned_params_updated_figs/";
marker = "_";

generate_conditioned_nrel_param_csv()

samples = DataFrame(CSV.File("conditioned_parameters_nrel_darco.csv"));
samples_shuffled = shuffle(MersenneTwister(123), samples);
samples_cut = samples_shuffled[1:N,:]

alg, stab_alg, unst_alg, alg_nans = get_small_signal_results_from_samples(file_name, samples_cut, gfm_bus, gfl_bus, sm_bus, false, false, false)

dyn, stab_dyn, unst_dyn, dyn_nans = get_small_signal_results_from_samples(file_name, samples_cut, gfm_bus, gfl_bus, sm_bus, true, false, false)

mssb, stab_mssb, unst_mssb, mssb_nans = get_small_signal_results_from_samples(file_name, samples_cut, gfm_bus, gfl_bus, sm_bus, true, true, false)

#msmb, stab_msmb, unst_msmb = get_small_signal_results_from_samples(file_name, samples_cut, gfm_bus, gfl_bus, sm_bus, true, true, true);

# for i = 1:100; 
#     if (dyn_nans[i,:stable] > 0) & (mssb_nans[i, :stable] < 0)
#         print(alg_nans[i,:])
#         print(dyn_nans[i,:])
#     elseif (dyn_nans[i, :stable] < 0) & (mssb_nans[i, :stable] > 0)
#         print(alg_nans[i,:])
#         print(dyn_nans[i,:])
#     end
# end

# Plot settings 
legend_loc = :outertopright

#all = [alg, dyn, mssb, msmb]
#stab = [stab_alg, stab_dyn, stab_mssb, stab_msmb]
#unstable = [unst_alg, unst_dyn, unst_mssb, unst_msmb]
labels = ["statpi", "dynpi", "mssb", "msmb"]

sel = [alg, dyn, mssb]
#sel = [stab_alg, stab_dyn, stab_mssb]
#sel = all # select which data to plot


# KQ 
w = 0.005
 
@df sel[1] boxplot(:kq.+get_o(w,1), :stable, bar_width=w, label=labels[1], legend=legend_loc, xlabel=L"k_q", ylabel=L"R(\lambda)")
@df sel[2] boxplot!(:kq.+get_o(w,2), :stable, bar_width=w, label=labels[2])
@df sel[3] boxplot!(:kq.+get_o(w,3), :stable, bar_width=w, label=labels[3])
#@df sel[4] boxplot!(:kq.+get_o(w,4), :stable, bar_width=w, label=labels[4])

savefig(save_folder*marker*"kq.png")
savefig(save_folder*marker*"kq.svg")


# LOADING 
w = 0.02

@df sel[1] boxplot(:load_scale.+get_o(w,1), :stable, bar_width=w, label=labels[1], legend=legend_loc, xlabel=L"load\ scale", ylabel=L"R(\lambda)")
@df sel[2] boxplot!(:load_scale.+get_o(w,2), :stable, bar_width=w, label=labels[2])
@df sel[3] boxplot!(:load_scale.+get_o(w,3), :stable, bar_width=w, label=labels[3])
#@df remove_nans_from_df(sel[4]) boxplot!(:load_scale.+get_o(w,4), :stable, bar_width=w, label=labels[4])

savefig(save_folder*marker*"load_scale.png")
savefig(save_folder*marker*"load_scale.svg")


# KPV 
w = 0.002

@df sel[1] boxplot(:kpv.+get_o(w,1), :stable, bar_width=w, label=labels[1], legend=legend_loc, xlabel=L"k_{pv}", ylabel=L"R(\lambda)")
@df sel[2] boxplot!(:kpv.+get_o(w,2), :stable, bar_width=w, label=labels[2])
@df sel[3] boxplot!(:kpv.+get_o(w,3), :stable, bar_width=w, label=labels[3])
#@df remove_nans_from_df(sel[4]) boxplot!(:kpv.+get_o(w,4), :stable, bar_width=w, label=labels[4])

savefig(save_folder*marker*"kpv.png")
savefig(save_folder*marker*"kpv.svg")

# KIV 
w = 15.0

@df sel[1] boxplot(:kiv.+get_o(w,1), :stable, bar_width=w, label=labels[1], legend=legend_loc, xlabel=L"k_{iv}", ylabel=L"R(\lambda)")
@df sel[2] boxplot!(:kiv.+get_o(w,2), :stable, bar_width=w, label=labels[2])
@df sel[3] boxplot!(:kiv.+get_o(w,3), :stable, bar_width=w, label=labels[3])
#@df remove_nans_from_df(sel[4]) boxplot!(:kiv.+get_o(w,4), :stable, bar_width=w, label=labels[4])

savefig(save_folder*marker*"kiv.png")
savefig(save_folder*marker*"kiv.svg")

# KPC 
w = 0.01;

@df sel[1] boxplot(:kpc.+get_o(w,1), :stable, bar_width=w, label=labels[1], legend=legend_loc, xlabel=L"k_{pc}", ylabel=L"R(\lambda)")
@df sel[2] boxplot!(:kpc.+get_o(w,2), :stable, bar_width=w, label=labels[2])
@df sel[3] boxplot!(:kpc.+get_o(w,3), :stable, bar_width=w, label=labels[3])
#@df remove_nans_from_df(sel[4]) boxplot!(:kpc.+get_o(w,4), :stable, bar_width=w, label=labels[4])

savefig(save_folder*marker*"kpc.png")
savefig(save_folder*marker*"kpc.svg")


# KIC 
w = 0.2;

@df sel[1] boxplot(:kic.+get_o(w,1), :stable, bar_width=w, label=labels[1], legend=legend_loc, xlabel=L"k_{ic}", ylabel=L"R(\lambda)")
@df sel[2] boxplot!(:kic.+get_o(w,2), :stable, bar_width=w, label=labels[2])
@df sel[3] boxplot!(:kic.+get_o(w,3), :stable, bar_width=w, label=labels[3])
#@df remove_nans_from_df(sel[4]) boxplot!(:kic.+get_o(w,4), :stable, bar_width=w, label=labels[4])

savefig(save_folder*marker*"kic.png")
savefig(save_folder*marker*"kic.svg")

# GFM eta  
w = 0.02

@df bin_data(sel[1],:gfm_eta) boxplot(:gfm_eta.+get_o(w,1), :stable, bar_width=w, label=labels[1], legend=legend_loc, xlabel=L"\eta_{gfm}", ylabel=L"R(\lambda)")
@df bin_data(sel[2], :gfm_eta) boxplot!(:gfm_eta.+get_o(w,2), :stable, bar_width=w, label=labels[2])
@df bin_data(sel[3], :gfm_eta) boxplot!(:gfm_eta.+get_o(w,3), :stable, bar_width=w, label=labels[3])
#@df remove_nans_from_df(bin_data(sel[4], :gfm_eta)) boxplot!(:gfm_eta.+get_o(w,4), :stable, bar_width=w, label=labels[4])

savefig(save_folder*marker*"gfm_eta.png")
savefig(save_folder*marker*"gfm_eta.svg")


# GFL eta  
w = 0.02

@df bin_data(sel[1],:gfl_eta) boxplot(:gfl_eta.+get_o(w,1), :stable, bar_width=w, label=labels[1], legend=legend_loc, xlabel=L"\eta_{gfl}", ylabel=L"R(\lambda)")
@df bin_data(sel[2], :gfl_eta) boxplot!(:gfl_eta.+get_o(w,2), :stable, bar_width=w, label=labels[2])
@df bin_data(sel[3], :gfl_eta) boxplot!(:gfl_eta.+get_o(w,3), :stable, bar_width=w, label=labels[3])
#@df remove_nans_from_df(bin_data(sel[4], :gfl_eta)) boxplot!(:gfl_eta.+get_o(w,4), :stable, bar_width=w, label=labels[4])

savefig(save_folder*marker*"gfl_eta.png")
savefig(save_folder*marker*"gfl_eta.svg")

# kd
w = 10
 
@df sel[1] boxplot(:kd.+get_o(w,1), :stable, bar_width=w, label=labels[1], legend=legend_loc, xlabel=L"k_d", ylabel=L"R(\lambda)")
@df sel[2] boxplot!(:kd.+get_o(w,2), :stable, bar_width=w, label=labels[2])
@df sel[3] boxplot!(:kd.+get_o(w,3), :stable, bar_width=w, label=labels[3])
#@df remove_nans_from_df(sel[4]) boxplot!(:kd.+get_o(w,4), :stable, bar_width=w, label=labels[4])

savefig(save_folder*marker*"kd.png")
savefig(save_folder*marker*"kd.svg")

# Ta
w = 0.05

@df sel[1] boxplot(:Ta.+get_o(w,1), :stable, bar_width=w, label=labels[1], legend=legend_loc, xlabel=L"T_a", ylabel=L"R(\lambda)")
@df sel[2] boxplot!(:Ta.+get_o(w,2), :stable, bar_width=w, label=labels[2])
@df sel[3] boxplot!(:Ta.+get_o(w,3), :stable, bar_width=w, label=labels[3])
#@df sel[4] boxplot!(:Ta.+get_o(w,4), :stable, bar_width=w, label=labels[4])

savefig(save_folder*marker*"Ta.png")
savefig(save_folder*marker*"Ta.svg")


# SAVE DATA 
CSV.write(save_folder*"alg.csv", alg)
CSV.write(save_folder*"dyn.csv", dyn)
CSV.write(save_folder*"mssb.csv", mssb)