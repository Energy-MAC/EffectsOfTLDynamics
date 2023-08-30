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
l = 500; # Line length (km)
line_dict["BUS 1-BUS 2-i_1"] = l;


# Define line parameters 
Z_c = 380.0; # Ω
r_km = 0.05; # Ω/km
x_km = 0.488; # Ω/km
g_km = 0; # S/km
b_km = 3.371e-6; # S/km
z_km = (r_km + im*x_km)*3;
y_km = im*b_km;


M = 1; # Number of parallel branches 

perturbation_type = "CRC"
t_fault = 0.25
perturbation_params = get_default_perturbation(t_fault, perturbation_type)

Nrange = [1,50,100];

stb = [];
max_λ = []
second_λ = []
plot()
for N in Nrange;
    
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
           (0.0, 2.0), #time span
           dist, #Type of perturbation
           all_lines_dynamic = true
       )

    ss = small_signal_analysis(sim)
    push!(stb, ss)
    push!(max_λ, maximum(real(ss.eigenvalues)))
    n_eigs = length(ss.eigenvalues)
    push!(second_λ, ss.eigenvalues[n_eigs-2])
    
    # # Plot eigenvalues 
    display(plot!(real(ss.eigenvalues), imag(ss.eigenvalues), seriestype=:scatter, label="N="*string(N)))

end

max_λ
second_λ


xlabel!("Real")
ylabel!("Imag")
xlims!((-.36234, -0.3623))
xlims!((-4.28445, -4.2844))
title!("L="*string(l)*"km")

xlims!(-1,0.1)
