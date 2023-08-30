cd(@__DIR__)
using PowerSystems
using PowerSimulationsDynamics
using Plots
using Revise

using CSV
using DataFrames

const PSY = PowerSystems;
const PSID = PowerSimulationsDynamics;

include("../src/ExperimentStructs.jl")
include("../src/experiment_funcs.jl")

file_name = "../data/json_data/inv_v_machine.json"; # choose system 
line_dict = default_2_bus_line_dict
l = 100; # Line length (km)
line_dict["BUS 1-BUS 2-i_1"] = l;

# Define line parameters 
# Kundur params 
Z_c = 380.0; # Ω
r_km = 0.05; # Ω/km
x_km = 0.488; # Ω/km
g_km = 0; # S/km
b_km = 3.371e-6; # S/km
z_km = r_km + im*x_km;
y_km = im*b_km;

M = 1; # Number of parallel branches 
N = 1; # Number of segments 

# Define perturbation 
perturbation_type = "CRC"
t_fault = 0.25
perturbation_params = get_default_perturbation(t_fault, perturbation_type)

p = ExpParams(
    N, 
    M, 
    l, 
    Z_c, 
    z_km,
    y_km,
    z_km, 
    line_dict,
    sim_p,
    perturbation_type, 
    perturbation_params)

sys = System(joinpath(pwd(), file_name));
sys_ms = build_seg_model!(sys, p)

dist = choose_disturbance(sys_ms, perturbation_type, p)


sim = PSID.Simulation(
        MassMatrixModel, #Type of model used
        sys_ms, #system
        pwd(), #folder to output results
        (0.0, 20.0), #time span
        dist, #Type of perturbation
        all_lines_dynamic = true
    )

execute_sim!(sim, p);
results = results_sim(sim);

#v1 = get_state_series(results, ("generator-102-1", :vr_filter));
v1d = get_state_series(results, ("V_1", :R));
v1q = get_state_series(results, ("V_1", :I));

v2d = get_state_series(results, ("V_2", :R));
v2q = get_state_series(results, ("V_2", :I));

plot(grid=true)
#plot!(v1, xlabel = "time", ylabel = "vr p.u.", label = "vr")
plot!(v1d, xlabel="time", ylabel="p.u.", label="v1d, N="*string(N))
plot!(v2d, xlabel="time", label="v2d, N="*string(N))

xlims!(0.249,0.253)
ylims!(1.035, 1.051)
xlims!(0,2.0)
ylims!(0.95, 1.05)

xlims!(0.249,0.252)
ylims!(1.02, 1.05)


#title!("N="*string(N)*", L="*string(l)*"km")
title!("L="*string(l)*"km")
ylabel!("voltage p.u.")
xlabel!("time")
