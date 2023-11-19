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

N = 1000; # Number of samples of parameters 

samples = generate_single_case_gain_params()
samples_cut = get_random_sample(samples, N)

exp_folder = "Nov6_1000_samples_lseg_10" # Make folder for saving exp results 
mkdir(exp_folder)

# Save sample in folder 
CSV.write(exp_folder*"/samples.csv", samples_cut)

# Define DF of experiment parameters 
exp_df = DataFrame(CSV.File("exp_df.csv"))

for i = 1:size(exp_df)[1];
    line_scale = exp_df[i,:line_scale]
    inv_share = exp_df[i,:inv_share]
    load_scale = exp_df[i,:load_scale]
    gfl_share = exp_df[i,:gfl_share]
    case = exp_df[i,:case_name]

    generate_single_case_3(samples_cut, inv_share, load_scale, gfl_share, line_scale, file_name, gfm_bus, gfl_bus, sm_bus, exp_folder, case)

end

# Plotting 
case1_scale = 0.1;
case2_scale = 0.02;
case3_scale = 0.4;
case4_scale = 0.2;

yscales = [case1_scale, case1_scale, case2_scale, case2_scale*10, case3_scale, case3_scale, case4_scale, case4_scale]
legend_key = [:bottomright, false, false, false, false, false, false, false]

include("9bus_analysis_functions.jl")
# Make boxplots 
for i = 1:size(exp_df)[1];
    load_scale = exp_df[i,:load_scale]
    case = exp_df[i,:case_name]
    case_name = case*"_"*string(load_scale)
    suffix = case*"_"*string(load_scale)*".csv"
    alg = DataFrame(CSV.File(exp_folder*"/alg_"*suffix))
    dyn = DataFrame(CSV.File(exp_folder*"/dyn_"*suffix))
    mssb = DataFrame(CSV.File(exp_folder*"/mssb_"*suffix))

    make_case_boxplot_3(alg, dyn, mssb, legend_key[i], case_name, exp_folder, yscales[i])
end

# choose case to make gain plots for 
i = 3
load_scale = exp_df[i,:load_scale]
case = exp_df[i,:case_name]
case_name = case*"_"*string(load_scale)
suffix = case*"_"*string(load_scale)*".csv"
alg = DataFrame(CSV.File(exp_folder*"/alg_"*suffix))
dyn = DataFrame(CSV.File(exp_folder*"/dyn_"*suffix))
mssb = DataFrame(CSV.File(exp_folder*"/mssb_"*suffix))


include("9bus_analysis_functions.jl")
make_gain_boxplots([alg, dyn, mssb], [L"\mathrm{statpi}", L"\mathrm{dynpi}", L"\mathrm{MSSB}"], :inside, case_name, exp_folder)
