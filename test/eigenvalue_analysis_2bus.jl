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

file_name = "../data/json_data/inv_v_machine.json"; # choose system 
line_dict = default_2_bus_line_dict

impedance_csv = "../data/cable_data/dommel_data.csv"
capacitance_csv = "../data/cable_data/dommel_data_C.csv"

perturbation_type = "CRC"
t_fault = 0.25
perturbation_params = get_default_perturbation(t_fault, perturbation_type)

### Define simulation parameters
sim_p = SimParams(
    abstol = 1e-13,
    reltol = 1e-10,
    maxiters = Int(1e10),
    dtmax = 1e-4,
    solver = "Rodas4",
    t_max = 0.1,
)

### Extract line data from files
M = 3

factor_z = 1.0 # to get it closer to kundur 
factor_y = 1.0
z_km, y_km, Z_c_abs, z_km_ω, z_km_ω_5_to_1, Z_c_5_to_1_abs = get_line_parameters(impedance_csv, capacitance_csv, M, factor_z, factor_y)


# # Calculate SIL  
# Z_c_total = sqrt(z_km_ω_5_to_1/y_km[1]);
# Vnom = 230; # KV 
# SIL = (Vnom^2)/Z_c_total

p_load = 0.5;
q_load = 0.5;

load_scale = 1.0
line_scale = 1.0

lseg_max = [10, 50, 100];

# p=0.5, q=0.5
#lRange = 475:5:580; 
#lRange = 150:10:250;
#lRange = 220:10:300;
#lRange = 100:50:500;
lRange = 900:20:1100;

mssb_results = zeros(length(lRange),length(lseg_max));
msmb_results = zeros(length(lRange),length(lseg_max));
dyn_results = [];
alg_results = [];
l_idx = 1;
for l = lRange;
    line_dict["BUS 1-BUS 2-i_1"] = l;
    line_dict["BUS 1-BUS 2-i_1_static"] = l;
    
    n_idx = 1;
    l_seg = l;
    
    p = ExpParams(
        nothing, 
        M, 
        l,
        l_seg, 
        Z_c_abs, 
        z_km,
        y_km,
        z_km_ω,
        z_km_ω_5_to_1,
        Z_c_5_to_1_abs,
        line_dict,
        sim_p, 
        perturbation_type, 
        perturbation_params,
        p_load,
        q_load,
        line_scale,
        load_scale
    )

    sim_alg = build_2bus_sim_from_file(file_name, false, false, p)

    if sim_alg.status != PSID.BUILD_FAILED;
        ss_alg = small_signal_analysis(sim_alg)
        if maximum(real(ss_alg.eigenvalues)) == 0.0;
            # Choose second largest eig 
            push!(alg_results, real(ss_alg.eigenvalues[end-1]))
        else
            push!(alg_results, maximum(real(ss_alg.eigenvalues)))
        end
    else
        push!(alg_results, NaN)
    end

    # Do dynamic pi  
    sim_dyn = build_2bus_sim_from_file(file_name, true, false, p)

    if sim_dyn.status != PSID.BUILD_FAILED;
        ss_dyn = small_signal_analysis(sim_dyn)
        if maximum(real(ss_dyn.eigenvalues)) == 0.0;
            # Choose second largest eig 
            push!(dyn_results, real(ss_dyn.eigenvalues[end-1]))
        else
            push!(dyn_results, maximum(real(ss_dyn.eigenvalues)))
        end
    else
        push!(dyn_results, NaN)
    end

    

    for l_seg in lseg_max;

        # Do MSSB 

        p = ExpParams(
            nothing, 
            1, # M=1
            l,
            l_seg, 
            Z_c_abs, 
            z_km,
            y_km,
            z_km_ω,
            z_km_ω_5_to_1,
            Z_c_5_to_1_abs,
            line_dict,
            sim_p, 
            perturbation_type, 
            perturbation_params,
            p_load,
            q_load,
            line_scale,
            load_scale
        )

        sim_ms = build_2bus_sim_from_file(file_name, true, true, p)

        if sim_ms.status != PSID.BUILD_FAILED
            ss_ms = small_signal_analysis(sim_ms)
            if maximum(real(ss_ms.eigenvalues)) == 0.0;
                # Choose second largest eig 
                mssb_results[l_idx, n_idx] = real(ss_ms.eigenvalues[end-1])
            else
                mssb_results[l_idx, n_idx] = maximum(real(ss_ms.eigenvalues))
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
            Z_c_abs, 
            z_km,
            y_km,
            z_km_ω,
            z_km_ω_5_to_1,
            Z_c_5_to_1_abs,
            line_dict,
            sim_p, 
            perturbation_type, 
            perturbation_params,
            p_load,
            q_load,
            line_scale,
            load_scale
        )

        sim_msmb = build_2bus_sim_from_file(file_name, true, true, p)

        if sim_msmb.status != PSID.BUILD_FAILED
            ss_msmb = small_signal_analysis(sim_msmb)
            if maximum(real(ss_msmb.eigenvalues)) == 0.0;
                # Choose second largest eig 
                msmb_results[l_idx, n_idx] = real(ss_msmb.eigenvalues[end-1])
            else
                msmb_results[l_idx, n_idx] = maximum(real(ss_msmb.eigenvalues))
            end
        else
            msmb_results[l_idx, n_idx] = NaN; # to indicate sim that hasn't run properly 
        end
        n_idx += 1;  
    end
    l_idx += 1
end


save_name = "loading p="*string(p_load)*" q="*string(q_load)

plot(lRange, alg_results, label="Algebraic",linestyle=:dash, legend = :outertopright, size=(800,600),xlabel="Line length (km)", ylabel="Largest real λ, != 0")
plot!(lRange, dyn_results, label="Dynpi",linestyle=:dash)
for k in 1:length(lseg_max);
    display(plot!(lRange,mssb_results[:,k],label="MSSB: lseg_max="*string(lseg_max[k])))
    display(plot!(lRange,msmb_results[:,k],label="MSMB: lseg_max="*string(lseg_max[k])))
end

function get_max(itr)
    return maximum([x for x in itr if !isnan(x)])
end

plt_ub = maximum([get_max(alg_results), get_max(dyn_results), get_max(mssb_results), get_max(msmb_results)])

plot!(lRange, zeros(length(lRange)), fillrange = ones(length(lRange))*plt_ub, fillalpha = 0.1, linealpha=0, c = 1, label="Unstable")
   
#title!("Kundur params, loading: p="*string(p_load)*", q="*string(q_load))

title!("p="*string(p_load)*"\nq="*string(q_load))

#savefig("../figures/Ruth/dommel_params/loading1_stability_margin.png")

# Heatmaps
clims = (0,0.2)

# MSSB
h1 = heatmap(lRange, lseg_max, mssb_results', clim=clims,colorbar_title="Largest real λ", size=(800,600), ylabel="MSSB - max l_seg", title="p="*string(p_load)*"\nq="*string(q_load))

# MSMB
h2 = heatmap(lRange, lseg_max, msmb_results', clim=clims,colorbar_title="Largest real λ",ylabel="MSMB - max l_seg")

# Dynpi 
h3 = heatmap(lRange, [1], dyn_results', clim=clims, colorbar_title="Largest real λ", ylabel="Dynpi")

# Algebraic 
h4 = heatmap(lRange, [1], alg_results', clim=clims,colorbar_title="Largest real λ", ylabel="Algebraic", xlabel="Line length (km)")

plot(h1,h2,h3,h4, layout=@layout[a;b;c;d], size=(1000,800))

savefig("../figures/Ruth/dommel_params/loading1_heatmap.png")
