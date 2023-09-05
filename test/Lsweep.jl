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

### Load relevant line data
impedance_csv = "../data/cable_data/impedance_data.csv"
capacitance_csv = "../data/cable_data/C_per_km.csv"


# Define line parameters 
Z_c = 380.0; # Ω
r_km = 0.05; # Ω/km
x_km = 0.488; # Ω/km
g_km = 0; # S/km
b_km = 3.371e-6; # S/km

f = 2;

z_km = (r_km + im*x_km);
z_km_ω = z_km;
y_km = im*b_km*f;

#M = 1; # Number of parallel branches 
#z_km, y_km, Z_c, z_km_ω = get_line_parameters(impedance_csv, capacitance_csv, M)

perturbation_type = "CRC"
t_fault = 0.25
perturbation_params = get_default_perturbation(t_fault, perturbation_type)

# lseg_opt = 60;
# # find closest N 
# Nup = ceil(l/lseg_opt);

Nrange = [1,2,3,5];
#Nrange = 1:10;
#lRange = [900,920,930,940,950];
#lRange = 730:5:755;


#lRange = 900:925;
#lRange = 600:25:800;
#lRange = 400:20:600;
#lRange = 495:3:515;

# normal inv_machine 
#lRange = 850:5:960;

# double Z, inv_v_machine
#lRange = 600:10:770;

# double Y, inv_machine
lRange = 450:5:550;

# lRange = 1050:1:1100;

results = zeros(length(lRange),length(Nrange));
dyn_results = [];
l_idx = 1;
for l = lRange;
    n_idx = 1;
    stb = [];
    max_λ = []
    second_λ = []
    line_dict["BUS 1-BUS 2-i_1"] = l;
    p = ExpParams(
        N, 
        M, 
        l, 
        Z_c, 
        z_km,
        y_km,
        z_km_ω, 
        line_dict,
        sim_p,
        perturbation_type, 
        perturbation_params)

    N = nothing
    sys = System(joinpath(pwd(), file_name));
    sys_dyn = build_new_impedance_model!(sys,p,true,"")
    sim_dyn = PSID.Simulation(
        MassMatrixModel, #Type of model used
        sys_dyn, #system
        pwd(), #folder to output results
        (0.0, 2.0), #time span
        dist, #Type of perturbation
        all_lines_dynamic = true
    )

    if sim_dyn.status != PSID.BUILD_FAILED;
        ss_dyn = small_signal_analysis(sim_dyn)
        push!(dyn_results, maximum(real(ss_dyn.eigenvalues)))
    end

    for N in Nrange;
        
        p = ExpParams(
        N, 
        M, 
        l, 
        Z_c, 
        z_km,
        y_km,
        z_km_ω, 
        line_dict,
        sim_p,
        perturbation_type, 
        perturbation_params)

        sys = System(joinpath(pwd(), file_name));
        sys_ms = build_seg_model!(sys, p, true, "")
        dist = choose_disturbance(sys_ms, perturbation_type, p)

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
            push!(stb, ss)
            push!(max_λ, maximum(real(ss.eigenvalues)))
            n_eigs = length(ss.eigenvalues)
            push!(second_λ, ss.eigenvalues[n_eigs-2])
            
            # # # Plot eigenvalues 
            # display(plot!(real(ss.eigenvalues), imag(ss.eigenvalues), seriestype=:scatter, label="N="*string(N)))
            #push!(results, [l, N, maximum(real(ss.eigenvalues))])
            results[l_idx, n_idx] = maximum(real(ss.eigenvalues));
        else
            results[l_idx, n_idx] = -1; # to indicate sim that hasn't run properly 
            
        end
        n_idx += 1;  
    end
    l_idx += 1
end

# Multiseg results 

h1 = heatmap(lRange, Nrange, results', clim=(0,0.1),colorbar_title="Largest real λ", size=(800,600))
# Dynamic pi results
h2 = heatmap(lRange, [1], dyn_results', clim=(0,0.1), colorbar_title="Largest real λ",)

plot(h1,h2, layout=@layout[a;b])

xlabel!("Line length (km)")
ylabel!("N")
title!("Stability heatmap")
# max_λ
# second_λ
# stb


xlabel!("Real")
ylabel!("Imag")
xlims!((-.36234, -0.3623))
xlims!((-4.28445, -4.2844))
title!("L="*string(l)*"km")

xlims!(-1,0.1)
