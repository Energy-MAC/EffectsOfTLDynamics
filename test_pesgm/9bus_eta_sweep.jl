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

file_name = "../data/json_data/9bus_VSM3_SM1_GFL2.json"

# defaults 
gfm_bus = "generator-3-1";
gfl_bus = "generator-2-1";
sm_bus = "generator-1-1";

N = 300;

samples = DataFrame(CSV.File("nrel_parameters.csv"))
samples_shuffled = shuffle(MersenneTwister(123), samples)
samples_cut = samples_shuffled[1:N,:]

@time begin
alg, stab_alg, unst_alg = get_small_signal_results_from_samples(file_name, samples_cut, gfm_bus, gfl_bus, sm_bus, false, false, false)
end

@time begin
dyn, stab_dyn, unst_dyn = get_small_signal_results_from_samples(file_name, samples_cut, gfm_bus, gfl_bus, sm_bus, true, false, false)
end

@time begin
mssb, stab_mssb, unst_mssb = get_small_signal_results_from_samples(file_name, samples_cut, gfm_bus, gfl_bus, sm_bus, true, true, false)
end 

msmb, stab_msmb, unst_msmb = get_small_signal_results_from_samples(file_name, samples_cut, gfm_bus, gfl_bus, sm_bus, true, true, true)

#mssb_df = filter(:stable => x -> !any(f -> f(x), (ismissing, isnan)), mssb_df)

# Plot settings 
legend_loc = :outertopright

all = [alg, dyn, mssb, msmb]
stab = [stab_alg, stab_dyn, stab_mssb, stab_msmb]
unstable = [unst_alg, unst_dyn, unst_mssb, unst_msmb]
labels = ["statpi", "dynpi", "mssb", "msmb"]

sel = all # select which data to plot

# KQ 
w = 0.005
offsets = [w*1.5, w/2, -w/2, -w*1.5]
 
@df remove_nans_from_df(sel[1]) boxplot(:kq.+offsets[1], :stable, bar_width=w, label=labels[1], legend=legend_loc, xlabel="Kq", ylabel=L"R(\lambda)")
@df remove_nans_from_df(sel[2]) boxplot!(:kq.+offsets[2], :stable, bar_width=w, label=labels[2])
@df remove_nans_from_df(sel[3]) boxplot!(:kq.+offsets[3], :stable, bar_width=w, label=labels[3])
@df remove_nans_from_df(sel[4]) boxplot!(:kq.+offsets[4], :stable, bar_width=w, label=labels[4])

# LOADING 
w = 0.02
offsets = [w*1.5, w/2, -w/2, -w*1.5]
 

@df remove_nans_from_df(sel[1]) boxplot(:load_scale.+offsets[1], :stable, bar_width=w, label=labels[1], legend=legend_loc, xlabel="load_scale", ylabel=L"R(\lambda)")
@df remove_nans_from_df(sel[2]) boxplot!(:load_scale.+offsets[2], :stable, bar_width=w, label=labels[2])
@df remove_nans_from_df(sel[3]) boxplot!(:load_scale.+offsets[3], :stable, bar_width=w, label=labels[3])
@df remove_nans_from_df(sel[4]) boxplot!(:load_scale.+offsets[4], :stable, bar_width=w, label=labels[4])


# KPV 
w = 0.0002
offsets = [w*1.5, w/2, -w/2, -w*1.5]

@df remove_nans_from_df(sel[1]) boxplot(:kpv.+offsets[1], :stable, bar_width=w, label=labels[1], legend=legend_loc, xlabel="kpv", ylabel=L"R(\lambda)")
@df remove_nans_from_df(sel[2]) boxplot!(:kpv.+offsets[2], :stable, bar_width=w, label=labels[2])
@df remove_nans_from_df(sel[3]) boxplot!(:kpv.+offsets[3], :stable, bar_width=w, label=labels[3])
@df remove_nans_from_df(sel[4]) boxplot!(:kpv.+offsets[4], :stable, bar_width=w, label=labels[4])

# KIV 
w = 15.0
offsets = [w*1.5, w/2, -w/2, -w*1.5]

@df remove_nans_from_df(sel[1]) boxplot(:kiv.+offsets[1], :stable, bar_width=w, label=labels[1], legend=legend_loc, xlabel="kiv", ylabel=L"R(\lambda)")
@df remove_nans_from_df(sel[2]) boxplot!(:kiv.+offsets[2], :stable, bar_width=w, label=labels[2])
@df remove_nans_from_df(sel[3]) boxplot!(:kiv.+offsets[3], :stable, bar_width=w, label=labels[3])
@df remove_nans_from_df(sel[4]) boxplot!(:kiv.+offsets[4], :stable, bar_width=w, label=labels[4])

# KPC 
w = 0.01;
offsets = [w*1.5, w/2, -w/2, -w*1.5]

@df remove_nans_from_df(sel[1]) boxplot(:kpc.+offsets[1], :stable, bar_width=w, label=labels[1], legend=legend_loc, xlabel="kpc", ylabel=L"R(\lambda)")
@df remove_nans_from_df(sel[2]) boxplot!(:kpc.+offsets[2], :stable, bar_width=w, label=labels[2])
@df remove_nans_from_df(sel[3]) boxplot!(:kpc.+offsets[3], :stable, bar_width=w, label=labels[3])
@df remove_nans_from_df(sel[4]) boxplot!(:kpc.+offsets[4], :stable, bar_width=w, label=labels[4])

# KIC 
w = 0.2;
offsets = [w*1.5, w/2, -w/2, -w*1.5]

@df remove_nans_from_df(sel[1]) boxplot(:kic.+offsets[1], :stable, bar_width=w, label=labels[1], legend=legend_loc, xlabel="kic", ylabel=L"R(\lambda)")
@df remove_nans_from_df(sel[2]) boxplot!(:kic.+offsets[2], :stable, bar_width=w, label=labels[2])
@df remove_nans_from_df(sel[3]) boxplot!(:kic.+offsets[3], :stable, bar_width=w, label=labels[3])
@df remove_nans_from_df(sel[4]) boxplot!(:kic.+offsets[4], :stable, bar_width=w, label=labels[4])



# GFM eta  
w = 0.005
offsets = [w*1.5, w/2, -w/2, -w*1.5]

@df remove_nans_from_df(bin_data(sel[1],:gfm_eta)) boxplot(:gfm_eta.+offsets[1], :stable, bar_width=w, label=labels[1], legend=legend_loc, xlabel="gfm_eta", ylabel=L"R(\lambda)")
@df remove_nans_from_df(bin_data(sel[2], :gfm_eta)) boxplot!(:gfm_eta.+offsets[2], :stable, bar_width=w, label=labels[2])
@df remove_nans_from_df(bin_data(sel[3], :gfm_eta)) boxplot!(:gfm_eta.+offsets[3], :stable, bar_width=w, label=labels[3])
@df remove_nans_from_df(bin_data(sel[4], :gfm_eta)) boxplot!(:gfm_eta.+offsets[4], :stable, bar_width=w, label=labels[4])

# GFL eta  
w = 0.01
offsets = [w*1.5, w/2, -w/2, -w*1.5]

@df remove_nans_from_df(bin_data(sel[1],:gfl_eta)) boxplot(:gfl_eta.+offsets[1], :stable, bar_width=w, label=labels[1], legend=legend_loc, xlabel="gfm_eta", ylabel=L"R(\lambda)")
@df remove_nans_from_df(bin_data(sel[2], :gfl_eta)) boxplot!(:gfl_eta.+offsets[2], :stable, bar_width=w, label=labels[2])
@df remove_nans_from_df(bin_data(sel[3], :gfl_eta)) boxplot!(:gfl_eta.+offsets[3], :stable, bar_width=w, label=labels[3])
@df remove_nans_from_df(bin_data(sel[4], :gfl_eta)) boxplot!(:gfl_eta.+offsets[4], :stable, bar_width=w, label=labels[4])

