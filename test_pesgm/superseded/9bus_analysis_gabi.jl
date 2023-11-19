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

N = 2; # Number of samples of parameters 
exp_folder = "gabi";

samples = generate_single_case_gain_params()
samples_cut = get_random_sample(samples, N)

# Save sample of parameters 
CSV.write(exp_folder*"/samples.csv", samples_cut)

exp_df = DataFrame(CSV.File("gabi/exp_df.csv"))

for i = 1:size(exp_df)[1];
    line_scale = exp_df[i,:line_scale]
    inv_share = exp_df[i,:inv_share]
    load_scale = exp_df[i,:load_scale]
    gfl_share = exp_df[i,:gfl_share]
    case = exp_df[i,:case_name]

    alg, dyn, mssb, msmb, case_name = generate_single_case_full_msmb(samples_cut, inv_share, load_scale, gfl_share, line_scale, file_name, gfm_bus, gfl_bus, sm_bus, exp_folder, case)

end

### WORK IN PROGRESS - DON'T NEED TO RUN THIS YET ### 

# # load and plot saved data 
# for i = 1:size(exp_df)[1];
#     load_scale = exp_df[i,:load_scale]
#     case = exp_df[i,:case_name]
#     alg = DataFrame(CSV.File(exp_folder*"/alg_"*case_name*".csv"))
#     dyn = DataFrame(CSV.File(exp_folder*"/dyn_"*case_name*".csv"))
#     mssb = DataFrame(CSV.File(exp_folder*"/mssb"*case_name*".csv"))
#     msmb = DataFrame(CSV.File(exp_folder*"/msmb_"*case_name*".csv"))

#     # make_case_boxplot_full_msmb(alg, dyn, mssb, msmb, legend_loc, case_name, exp_folder)
#     # make_gain_boxplots([alg, dyn, mssb, msmb], ["statpi", "dynpi", "mssb", "msmb"], legend_loc, case_name, exp_folder)

# end

i = 1;
load_scale = exp_df[i,:load_scale]
case = exp_df[i,:case_name]
suffix = string(case)*"_"*string(load_scale)*".csv"
alg = DataFrame(CSV.File(exp_folder*"/alg_"*suffix))
dyn = DataFrame(CSV.File(exp_folder*"/dyn_"*suffix))
mssb = DataFrame(CSV.File(exp_folder*"/mssb_"*suffix))
msmb = DataFrame(CSV.File(exp_folder*"/msmb_"*suffix))


make_case_boxplot_full_msmb(alg, dyn, mssb, msmb, legend_loc, case_name, exp_folder)