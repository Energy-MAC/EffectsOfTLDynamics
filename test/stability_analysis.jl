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

impedance_csv = "../data/cable_data/impedance_data.csv"
capacitance_csv = "../data/cable_data/C_per_km.csv"


# Define line parameters 
Z_c = 380.0; # Ω
r_km = 0.05; # Ω/km
x_km = 0.488; # Ω/km
g_km = 0; # S/km
b_km = 3.371e-6; # S/km
z_km = r_km + im*x_km;
y_km = im*b_km;

N = 1;
M = 1; # Number of parallel branches 

perturbation_type = "CRC"
t_fault = 0.25
perturbation_params = get_default_perturbation(t_fault, perturbation_type)


#z_km, y_km, Z_c, z_km_ω = get_line_parameters(impedance_csv, capacitance_csv, M)

p = ExpParams(N, M, l, abs(Z_c), z_km,y_km, z_km, line_dict,sim_p, perturbation_type,perturbation_params)
sys = System(joinpath(pwd(), file_name));
dist = choose_disturbance(sys, perturbation_type, p)
sys_ms = build_seg_model!(sys, p)

sim = PSID.Simulation(
        MassMatrixModel, #Type of model used
        sys_ms, #system
        pwd(), #folder to output results
        (0.0, 2.0), #time span
        dist, #Type of perturbation
        all_lines_dynamic = true
    )

#frange = [0.1, 0.5, 2, 5, 10, 100]
frange = [1]
plot()
max_λ = [];
second_λ = [];
for f in frange;
    z_km_new = z_km;
    y_km_new = y_km;

    p = ExpParams(N, M, l, Z_c, z_km_new,y_km_new,z_km_new, line_dict,sim_p, perturbation_type,perturbation_params)

    sys = System(joinpath(pwd(), file_name));

    dist = choose_disturbance(sys, perturbation_type, p)

    sys_ms = build_seg_model!(sys, p)

    sim = PSID.Simulation(
            MassMatrixModel, #Type of model used
            sys_ms, #system
            pwd(), #folder to output results
            (0.0, 2.0), #time span
            dist, #Type of perturbation
            all_lines_dynamic = true
        )

    if sim.status != PSID.BUILD_FAILED
        ss = small_signal_analysis(sim)

        # # Plot eigenvalues 
        display(plot!(real(ss.eigenvalues), imag(ss.eigenvalues), seriestype=:scatter, label="f="*string(f)))
        push!(max_λ, maximum(real(ss.eigenvalues)))
        n_eigs = length(ss.eigenvalues)
        push!(second_λ, ss.eigenvalues[n_eigs-2])
    end

end

max_λ
second_λ

xlabel!("Real")
ylabel!("Imag")

xlims!(-20,-17)