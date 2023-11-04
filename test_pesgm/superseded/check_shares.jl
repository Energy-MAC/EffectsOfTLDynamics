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


alg = DataFrame(CSV.File("Nov1_conditioned_params_1000_samples/alg.csv"))

#stab_alg = filter(:stable => x -> x < 0, alg)
alg_low_load = filter(:load_scale => x -> x==0.5, alg)
# want inv_share=0.4, gfm_share=0.2
x = filter(:inv_share => x -> x==0.4, alg_low_load)
maximum(x[!,:stable])
minimum(x[!,:stable])

# y = filter(:gfm_share => x -> x==0.2, x)
# minimum(stab_alg_low_load[!,:gfl_eta])

# stab_alg_low_load_low_gfl = filter(:gfl_eta => x -> x <0.1, stab_alg_low_load)
# stab_alg_low_load_high_gfl = filter(:gfl_eta => x -> x >0.1, stab_alg_low_load)


alg_nom_load = filter(:load_scale => x -> x ==1.1, alg)
x = filter(:inv_share => x -> x==0.4, alg_nom_load)
maximum(x[!,:stable])
minimum(x[!,:stable])

minimum(stab_alg_nom_load[!,:gfl_eta])

stab_alg_low_load_low_gfl = filter(:gfl_eta => x -> x <0.1, stab_alg_nom_load)
stab_alg_low_load_high_gfl = filter(:gfl_eta => x -> x >0.1, stab_alg_nom_load)

