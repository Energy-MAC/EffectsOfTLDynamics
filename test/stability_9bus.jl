cd(@__DIR__)
using PowerSystems
using PowerSimulationsDynamics
using Plots
using Revise
using LaTeXStrings
using EffectsOfTLDynamics

using CSV
using DataFrames

const ETL = EffectsOfTLDynamics;
const PSY = PowerSystems;
const PSID = PowerSimulationsDynamics;

include("eig_analysis_funcs.jl")

line_dict = default_9_bus_line_dict

impedance_csv = "../data/cable_data/dommel_data.csv"
capacitance_csv = "../data/cable_data/dommel_data_C.csv"

t_fault = 0.25
perturbation_type = "BranchTrip"
perturbation_params = get_default_perturbation(t_fault, perturbation_type)
perturbation_params.branch_trip_params = BTParam("Bus 7-Bus 5-i_1")

### Define simulation parameters
sim_p = SimParams(
    abstol = 1e-13,
    reltol = 1e-10,
    maxiters = Int(1e10),
    dtmax = 1e-4,
    solver = "Rodas4",
    t_max = 0.1,
)

### Extract line data from files for M=1 and M=3
M = 1
factor_z = 1.0; 
factor_y = 1.0;
z_km_1, y_km_1, Z_c_abs_1, z_km_ω_1, z_km_ω_5_to_1_1, Z_c_5_to_1_abs_1 = get_line_parameters(impedance_csv, capacitance_csv, M, factor_z, factor_y)

M = 3
z_km_3, y_km_3, Z_c_abs_3, z_km_ω_3, z_km_ω_5_to_1_3, Z_c_5_to_1_abs_3 = get_line_parameters(impedance_csv, capacitance_csv, M, factor_z, factor_y)


l_seg = 10.0; # km 

p1 = ETL.ExpParams(
    nothing, 
    1, 
    100, # doesn't matter for 9 bus 
    l_seg, 
    Z_c_abs_1, 
    z_km_1,
    y_km_1,
    z_km_ω_1,
    z_km_ω_5_to_1_1,
    Z_c_5_to_1_abs_1,
    line_dict,
    sim_p, 
    perturbation_type, 
    perturbation_params,
    0.0, # p load 
    0.0, # q load 
    1.0, # line scale 
    1.0 # load scale 
);

p3 = ETL.ExpParams(
    nothing, 
    3, # M
    100, # doesn't matter for 9 bus 
    l_seg, 
    Z_c_abs_3, 
    z_km_3,
    y_km_3,
    z_km_ω_3,
    z_km_ω_5_to_1_3,
    Z_c_5_to_1_abs_3,
    line_dict,
    sim_p, 
    perturbation_type, 
    perturbation_params,
    0.0, # p load 
    0.0, # q load 
    1.0, # line scale 
    1.0 # load scale 
);


file_name = "../data/json_data/9bus_slackless.json"; # choose system  

load_scale = 1.0;
line_scales = collect(1.0:7.0)

#line_scales = [0.2,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0,6.0, 7.0] 
alg, dyn, mssb, msmb = generate_9bus_line_sweep_data(file_name, p1, p3, load_scale, line_scales, l_seg)
plot_9bus_line_sweep(alg, dyn, mssb, msmb, line_scales, (800,600), 2, L"Load\ scale "*string(load_scale), :top)

#line_to_scale = "Bus 7-Bus 8-i_1";
line_to_scale = "Bus 5-Bus 4-i_1"; # line from inverter ?
alg, dyn, mssb, msmb = generate_9bus_individual_line_sweep_data(file_name, p1, p3, load_scale, line_scales, l_seg, line_to_scale)

plot_9bus_line_sweep(alg, dyn, mssb, msmb, line_scales, (800,600), 2, "Load scale "*string(load_scale), :top)



# Scaling all lines 
load_scale = 1.2;
line_scales = collect(1.0:5.0)
alg, dyn, mssb, msmb = generate_9bus_line_sweep_data(file_name, p1, p3, load_scale, line_scales, l_seg)

plot_9bus_line_sweep(alg, dyn, mssb, msmb, line_scales, (800,600), 2, "Load scale "*string(load_scale), :top)


#savefig(foldername*"/load_20_shortest_line_scaled.png")


# # Save all these results in a jld2 file in a new folder 
foldername = "../figures/Ruth/ninebus"

savepath = foldername*"/9bus_line_sweep.jld2"
@save savepath alg dyn mssb msmb load_scale line_scales p1 p3 l_seg 
