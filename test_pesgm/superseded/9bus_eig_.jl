cd(@__DIR__)
using PowerSystems
using PowerSimulationsDynamics
using Plots
using Revise
using EffectsOfTLDynamics
using LaTeXStrings
using DataFrames

const ETL = EffectsOfTLDynamics;
const PSY = PowerSystems;
const PSID = PowerSimulationsDynamics;

include("defaults.jl")
include("../test/eig_analysis_funcs.jl")

file_name = "../data/json_data/9bus_VSM_SM_GFL_.json"

p1 = p1_9bus;
p1.load_scale = 1.0;
p1.line_scale = 1.0;
tspan = (0.0, p1.sim_params.t_max)
perturbation = choose_disturbance(sys, p1.perturbation, p1)
alg_line_name = p1.perturbation_params.branch_trip_params.line_to_trip

sim = ETL.build_9bus_sim_from_file(file_name, false, false, p1)

summary_eigenvalues(small_signal_analysis(sim))

summary_participation_factors(small_signal_analysis(sim))[!,[:Name, :λ_47]]

stability = save_max_nonzero_eig!(sim, [])[1]

g3 = get_component(Generator, sim.sys, "generator-3-1")
g3.dynamic_injector.outer_control.states