cd(@__DIR__)
using PowerSystems
using PowerSimulationsDynamics
using Plots
using Revise
using LaTeXStrings
#using EffectsOfTLDynamics

include("../src/ExperimentStructs.jl")
include("../src/experiment_funcs.jl")

using CSV
using DataFrames

const PSY = PowerSystems;
const PSID = PowerSimulationsDynamics;

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

p_load = 0.0; # do these matter for 9bus? no.
q_load = 0.0;

load_scale = 1.0

Nrange = [10];

line_scales = [0.3, 0.8, 1, 1.5];

mssb_results = zeros(length(line_scales),length(Nrange));
msmb_results = zeros(length(line_scales),length(Nrange));
dyn_results = [];
alg_results = [];
l_idx = 1;
l = 1000; # does this matter for 9bus?
for line_scale = line_scales;
    
    n_idx = 1;
    l_seg = l;
    
    p = ExpParams(
        nothing, 
        1, 
        l,
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
        p_load,
        q_load,
        line_scale,
        load_scale
    );

    sim_alg = build_9bus_sim_from_file(file_name, false, false, p);

    if sim_alg.status != PSID.BUILD_FAILED;
        ss_alg = small_signal_analysis(sim_alg);
        if maximum(real(ss_alg.eigenvalues)) == 0.0;
            # Choose second largest eig 
            push!(alg_results, real(ss_alg.eigenvalues[end-1]));
        else
            push!(alg_results, maximum(real(ss_alg.eigenvalues)));
        end
    else
        push!(alg_results, NaN)
    end

    # Do dynamic pi  
    sim_dyn = build_9bus_sim_from_file(file_name, true, false, p);

    if sim_dyn.status != PSID.BUILD_FAILED;
        ss_dyn = small_signal_analysis(sim_dyn);
        if maximum(real(ss_dyn.eigenvalues)) == 0.0;
            # Choose second largest eig 
            push!(dyn_results, real(ss_dyn.eigenvalues[end-1]));
        else
            push!(dyn_results, maximum(real(ss_dyn.eigenvalues)));
        end
    else
        push!(dyn_results, NaN)
    end

    
    for N in Nrange;
        l_seg = l/N;

        # Do MSSB 

        p = ExpParams(
            nothing, 
            1, # M=1
            l,
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
            p_load,
            q_load,
            line_scale,
            load_scale
        );

        sim_ms = build_9bus_sim_from_file(file_name, true, true, p);

        if sim_ms.status != PSID.BUILD_FAILED
            ss_ms = small_signal_analysis(sim_ms);
            if maximum(real(ss_ms.eigenvalues)) == 0.0;
                # Choose second largest eig 
                mssb_results[l_idx, n_idx] = real(ss_ms.eigenvalues[end-1]);
            else
                mssb_results[l_idx, n_idx] = maximum(real(ss_ms.eigenvalues));
            end
        else
            mssb_results[l_idx, n_idx] = NaN; # to indicate sim that hasn't run properly 
        end

        # Do MSMB 
        p = ExpParams(
            nothing, 
            3, #
            l,
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
            p_load,
            q_load,
            line_scale,
            load_scale
        );

        sim_msmb = build_9bus_sim_from_file(file_name, true, true, p);

        if sim_msmb.status != PSID.BUILD_FAILED
            ss_msmb = small_signal_analysis(sim_msmb);
            if maximum(real(ss_msmb.eigenvalues)) == 0.0;
                # Choose second largest eig 
                msmb_results[l_idx, n_idx] = real(ss_msmb.eigenvalues[end-1]);
            else
                msmb_results[l_idx, n_idx] = maximum(real(ss_msmb.eigenvalues));
            end
        else
            msmb_results[l_idx, n_idx] = NaN; # to indicate sim that hasn't run properly 
        end
        n_idx += 1;  
    end
    l_idx += 1
end

#### Stability margin plots 
plot(line_scales, alg_results, label="Algebraic",linestyle=:dash, legend = :outertopright, size=(800,600),xlabel="Line scale", ylabel="Largest real λ, != 0")
plot!(line_scales, dyn_results, label="Dynpi",linestyle=:dash)
for k in 1:length(Nrange);
    display(plot!(line_scales, mssb_results[:,k],label="MSSB: N="*string(Nrange[k]),linestyle=:dashdotdot, linewidth=2))
    display(plot!(line_scales, msmb_results[:,k],label="MSMB: N="*string(Nrange[k]), linestyle=:dashdot,linewidth=2))
end

function get_max(itr)
    return maximum([x for x in itr if !isnan(x)])
end

plt_ub = maximum([get_max(alg_results), get_max(dyn_results), get_max(mssb_results), get_max(msmb_results)])

plot!(line_scales, zeros(length(line_scales)), fillrange = ones(length(lRange))*plt_ub, fillalpha = 0.1, linealpha=0, c = 1, label="Unstable")

title!("p="*string(p_load*load_scale)*", q="*string(q_load*load_scale))

savefig("../figures/Ruth/nine_bus/stability_margin_Nsweep_1.png")


###### Heatmaps
clims = (0,0.2)

# MSSB
h1 = heatmap(lRange, Nrange, mssb_results', clim=clims,colorbar_title="Largest real λ", size=(800,600), ylabel="MSSB - N", title="p="*string(p_load*load_scale)*"\nq="*string(q_load*load_scale));

# MSMB
h2 = heatmap(lRange, Nrange, msmb_results', clim=clims,colorbar_title="Largest real λ",ylabel="MSMB - N");

# Dynpi 
h3 = heatmap(lRange, [1], dyn_results', clim=clims, colorbar_title="Largest real λ", ylabel="Dynpi");

# Algebraic 
h4 = heatmap(lRange, [1], alg_results', clim=clims,colorbar_title="Largest real λ", ylabel="Algebraic", xlabel="Line length (km)");

plot(h1,h2,h3,h4, layout=@layout[a;b;c;d], size=(1000,800))

#savefig("../figures/Ruth/dommel_params/loading1_heatmap.png")


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

savefig("../figures/Ruth/dommel_params/eig_comparison_1.png")

# Analyse the participation in the least stable modes?
