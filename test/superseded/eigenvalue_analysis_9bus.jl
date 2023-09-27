cd(@__DIR__)
using PowerSystems
using PowerSimulationsDynamics
using Plots
using Revise
using LaTeXStrings
using EffectsOfTLDynamics

const ETL = EffectsOfTLDynamics
using CSV
using DataFrames

const PSY = PowerSystems;
const PSID = PowerSimulationsDynamics;

function save_max_nonzero_eig!(sim, output);
    if sim.status != PSID.BUILD_FAILED;
        ss = small_signal_analysis(sim);
        if maximum(real(ss.eigenvalues)) == 0.0;
            # Choose second largest eig 
            push!(output, real(ss.eigenvalues[end-1]));
        else
            push!(output, maximum(real(ss.eigenvalues)));
        end
    else
        push!(output, NaN)
    end
end

file_name = "../data/json_data/9bus_slackless.json"; # choose system 
line_dict = default_9_bus_line_dict

impedance_csv = "../data/cable_data/dommel_data.csv"
capacitance_csv = "../data/cable_data/dommel_data_C.csv"

# Define perturbation - not important for eig analysis 
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
factor_z = 1.0 # to get it closer to kundur 
factor_y = 1.0
z_km_1, y_km_1, Z_c_abs_1, z_km_ω_1, z_km_ω_5_to_1_1, Z_c_5_to_1_abs_1 = get_line_parameters(impedance_csv, capacitance_csv, M, factor_z, factor_y)

M = 3
factor_z = 1.0 # to get it closer to kundur 
factor_y = 1.0
z_km_3, y_km_3, Z_c_abs_3, z_km_ω_3, z_km_ω_5_to_1_3, Z_c_5_to_1_abs_3 = get_line_parameters(impedance_csv, capacitance_csv, M, factor_z, factor_y)

p1 = ETL.ExpParams(
    nothing, 
    1, 
    100,
    50, 
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
    0.0,
    0.0,
    1.0,
    1.0
);

p3 = ETL.ExpParams(
    nothing, 
    3, #
    100,
    50, 
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
    0.0,
    0.0,
    1.0,
    1.0
);

load_scales = [2];
line_scales = [[0.1:0.1:2]];

alg_lims = [];
dyn_lims = [];
mssb_lims = [];
msmb_lims = [];

idx = 1;
mssb_results = [];
msmb_results = [];
dyn_results = [];
alg_results = [];

for load_scale in load_scales;
    # mssb_results = [];
    # msmb_results = [];
    # dyn_results = [];
    # alg_results = [];
    p1.load_scale = load_scale
    p3.load_scale = load_scale
    current_lscale = line_scales[idx][1];
    for line_scale = current_lscale;

        p1.line_scale = line_scale
        p3.line_scale = line_scale

        sim_alg = ETL.build_9bus_sim_from_file(file_name, false, false, p1);
        save_max_nonzero_eig!(sim_alg, alg_results)

        # Dynamic pi  
        sim_dyn = ETL.build_9bus_sim_from_file(file_name, true, false, p1);
        save_max_nonzero_eig!(sim_dyn, dyn_results)

        # MSSB
        sim_ms = ETL.build_9bus_sim_from_file(file_name, true, true, p1);
        save_max_nonzero_eig!(sim_ms, mssb_results)

        # Do MSMB 
        sim_msmb = ETL.build_9bus_sim_from_file(file_name, true, true, p3);
        save_max_nonzero_eig!(sim_msmb, msmb_results)

    end

    # Find first length where we lose stability - NaN indicates that we don't encounter stability limit 
    if findfirst(alg_results.>0) !== nothing;
        push!(alg_lims, current_lscale[findfirst(alg_results.>0)])
    else
        push!(alg_lims, NaN)
    end
    if findfirst(dyn_results.>0) !== nothing;
        push!(dyn_lims, current_lscale[findfirst(dyn_results.>0)])
    else
        push!(dyn_lims, NaN)
    end
    if findfirst(mssb_results.>0) !== nothing;
        push!(mssb_lims, current_lscale[findfirst(mssb_results.>0)])
    else
        push!(mssb_lims, NaN)
    end
    if findfirst(msmb_results.>0) !== nothing;
        push!(msmb_lims, current_lscale[findfirst(msmb_results.>0)])
    else
        push!(msmb_lims, NaN)
    end
    idx += 1
end


plot(line_scales[1], alg_results,xlabel="Line scale", ylabel="Max real λ != 0", label="Algebraic", seriestype=:line, size=(800,600))
plot!(line_scales[1], dyn_results, seriestype=:line, label="Dynpi")
plot!(line_scales[1], mssb_results, seriestype=:line, linewidth=2, label="MSSB lseg="*string(p1.l_seg))
plot!(line_scales[1], msmb_results, seriestype=:line, linestyle=:dashdot, linewidth=4, label="MSMB lseg="*string(p1.l_seg))
title!("Load scale="*string(load_scales[1]))

#savefig("../figures/Ruth/nine_bus/line_length_v_eig_3.png")




plot(load_scales, alg_lims,xlabel="Load scale", ylabel="Line length @ stability boundary", label="Algebraic", seriestype=:line, size=(800,600))
plot!(load_scales, dyn_lims, seriestype=:line, label="Dynpi")
plot!(load_scales, mssb_lims, seriestype=:line, linestyle=:dashdot, label="MSSB N="*string(N))
plot!(load_scales, msmb_lims, seriestype=:line, linestyle=:dash, label="MSMB N="*string(N))

title!("Base p="*string(p_load)*", q="*string(q_load))
#





### Eigenvalue plots #####
# Select params
p_load = 1.0;
q_load = 0.0;

load_scale = 1.0
line_scale = 1.0

l_seg = 20;
l = 300; # not doing anything 

p1 = ExpParams(
    nothing, 
    1, 
    l,
    l_seg, 
    Z_c_abs, 
    z_km_1,
    y_km_1,
    z_km_ω_1,
    z_km_ω_5_to_1_1,
    Z_c_5_to_1_abs_1,
    line_dict,
    sim_p, 
    perturbation_type, 
    perturbation_params,
    p_load,
    q_load,
    line_scale,
    load_scale
)

p3 = ExpParams(
    nothing, 
    3, #
    l,
    l_seg, 
    Z_c_abs, 
    z_km_3,
    y_km_3,
    z_km_ω_3,
    z_km_ω_5_to_1_3,
    Z_c_5_to_1_abs_3,
    line_dict,
    sim_p, 
    perturbation_type, 
    perturbation_params,
    p_load,
    q_load,
    line_scale,
    load_scale
);

# Algebraic 
sim_alg = build_9bus_sim_from_file(file_name, false, false, p1)
eigs_alg = small_signal_analysis(sim_alg).eigenvalues;
sim_dyn = build_9bus_sim_from_file(file_name, true, false, p1)
eigs_dyn = small_signal_analysis(sim_dyn).eigenvalues;
sim_mssb = build_9bus_sim_from_file(file_name, true, true, p1)
eigs_mssb = small_signal_analysis(sim_mssb).eigenvalues;
sim_msmb = build_9bus_sim_from_file(file_name, true, true, p3)
eigs_msmb = small_signal_analysis(sim_msmb).eigenvalues;


plot(real(eigs_alg), imag(eigs_alg), seriestype=:scatter, xlabel="Real", ylabel="Imag", label="Algebraic",legend = :outertopright)
plot!(real(eigs_dyn), imag(eigs_dyn), seriestype=:scatter,label="Dynpi")
plot!(real(eigs_mssb), imag(eigs_mssb), seriestype=:scatter,label="MSSB")
plot!(real(eigs_msmb), imag(eigs_msmb), seriestype=:scatter,label="MSMB")

#xlims!(-1000,0.1)
#ylims!(-4,4)

title!("System eigs, p="*string(p_load*load_scale)*", q="*string(q_load*load_scale))

#savefig("../figures/Ruth/dommel_params/eig_comparison_1.png")

# Analyse the participation in the least stable modes?
