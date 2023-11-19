cd(@__DIR__)
using PowerSystems
using PowerSimulationsDynamics
using Plots
using Revise
using EffectsOfTLDynamics

const ETL = EffectsOfTLDynamics;
const PSY = PowerSystems;
const PSID = PowerSimulationsDynamics;

include("9bus_analysis_functions.jl")

# defaults 
file_name = "../data/json_data/9bus_VSM3_SM1_GFL2.json"
gfm_bus = "generator-3-1";
gfl_bus = "generator-2-1";
sm_bus = "generator-1-1";

# Parameters that stay constant across imulations
line_scale = 1.0;
inv_share = 0.4;
N = 200; # Number of samples of parameters 
legend_loc = :outerright

# Define parameters 
load_scale = 0.5;
gfl_share = 0.2;

alg, dyn, case_name = generate_single_case_min(N, inv_share, load_scale, gfl_share, line_scale, file_name, gfm_bus, gfl_bus, sm_bus)

make_case_boxplot_min(alg, dyn, legend_loc, case_name)
make_gain_boxplots([alg, dyn], ["statpi", "dynpi"], legend_loc, case_name)

# Define parameters 
load_scale = 1.0;
gfl_share = 0.2;

alg, dyn, case_name = generate_single_case_min(N, inv_share, load_scale, gfl_share, line_scale, file_name, gfm_bus, gfl_bus, sm_bus)

make_case_boxplot_min(alg, dyn, legend_loc, case_name)
make_gain_boxplots([alg, dyn], ["statpi", "dynpi"], legend_loc, case_name)

# Define parameters 
load_scale = 0.5;
gfl_share = 0.8;

alg, dyn, case_name = generate_single_case_min(N, inv_share, load_scale, gfl_share, line_scale, file_name, gfm_bus, gfl_bus, sm_bus)

make_case_boxplot_min(alg, dyn, legend_loc, case_name)
make_gain_boxplots([alg, dyn], ["statpi", "dynpi"], legend_loc, case_name)


# Define parameters 
load_scale = 1.0;
gfl_share = 0.8;

alg, dyn, case_name = generate_single_case_min(N, inv_share, load_scale, gfl_share, line_scale, file_name, gfm_bus, gfl_bus, sm_bus)

make_case_boxplot_min(alg, dyn, legend_loc, case_name)
make_gain_boxplots([alg, dyn], ["statpi", "dynpi"], legend_loc, case_name)






