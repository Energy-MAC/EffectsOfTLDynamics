cd(@__DIR__)
using PowerSystems
using PowerSimulationsDynamics
using Plots
using Revise
#using EffectsOfTLDynamics

include("../src/ExperimentStructs.jl")
include("../src/experiment_funcs.jl")

using CSV
using DataFrames

const PSY = PowerSystems;
const PSID = PowerSimulationsDynamics;

file_name = "../data/json_data/inv_v_machine.json"; # choose system 
line_dict = default_2_bus_line_dict

impedance_csv = "../data/cable_data/impedance_data.csv"
capacitance_csv = "../data/cable_data/C_per_km.csv"

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
M = 5

factor_z = 1.0 # to get it closer to kundur 
factor_y = 1.0
z_km, y_km, Z_c_abs, z_km_ω, z_km_ω_5_to_1, Z_c_5_to_1_abs = get_line_parameters(impedance_csv, capacitance_csv, M, factor_z, factor_y)

# # Replace with kundur
# Z_c_abs = 380 # Ω
# r_km = 0.05 # Ω/km
# x_km = 0.488 # Ω/km
# g_km = 0.0 # S/km
# b_km = 3.371e-6 # S/km
# z_km = r_km + im*x_km;
# y_km = im*b_km;
# z_km_ω = z_km;
# z_km_ω_5_to_1 = z_km;
# Z_c_5_to_1_abs = Z_c_abs;

p_load = 1.0;
q_load = 0.0;

Nrange = [1,2,5,10,20];

# p=0.5, q=0.5
lRange = 400:5:500;
lRange = 150:10:250;

# p=q=0
#lRange = 400:10:800;


# Kundur 
#lRange = 880:5:940; # P=1
#lRange = 1000:10:1100; # P=1


results = zeros(length(lRange),length(Nrange));
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
        q_load
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

    # Do MS 

    for N in Nrange;
        l_seg = l/N

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
            q_load
        )

        sim_ms = build_2bus_sim_from_file(file_name, true, true, p)

        if sim_ms.status != PSID.BUILD_FAILED
            ss_ms = small_signal_analysis(sim_ms)
            if maximum(real(ss_ms.eigenvalues)) == 0.0;
                # Choose second largest eig 
                results[l_idx, n_idx] = real(ss_ms.eigenvalues[end-1])
            else
                results[l_idx, n_idx] = maximum(real(ss_ms.eigenvalues))
            end
        else
            results[l_idx, n_idx] = NaN; # to indicate sim that hasn't run properly 
        end
        n_idx += 1;  
    end
    l_idx += 1
end

# title_str = "Stability heatmap, loading: p="*string(p_load)*", q="*string(q_load);
# minc = -.1;
# maxc = .1;
# h1 = heatmap(lRange, Nrange, results', clim=(minc, maxc),colorbar_title="Largest real λ", ylabel="N", title=title_str)

# # Dynamic pi results
# h2 = heatmap(lRange, [1], dyn_results', clim=(minc, maxc), colorbar_title="Largest real λ", ylabel="Dyn")

# #algebraic/static pi results
# h3 = heatmap(lRange, [1], alg_results', clim=(minc, maxc), colorbar_title="Largest real λ", ylabel="Alg", xlabel="Line length (km)");

# plot(h1,h2, h3, layout=@layout[a;b;c], size=(800,600))


plot(lRange, alg_results, label="Alg",linestyle=:dash)
plot!(lRange, dyn_results, label="Dyn",linestyle=:dash)
for k in 1:length(Nrange);
    display(plot!(lRange,results[:,k],label="N="*string(Nrange[k])*",M="*string(M)))
end

function get_max(itr)
    return maximum([x for x in itr if !isnan(x)])
end

plt_ub = maximum([get_max(alg_results), get_max(dyn_results), get_max(results)])

plot!(lRange, zeros(length(lRange)), fillrange = ones(length(lRange))*plt_ub, fillalpha = 0.1, linealpha=0, c = 1, label="Unstable")
   
xlabel!("Line length (km)")
ylabel!("Largest real λ, != 0")
#title!("Kundur params, loading: p="*string(p_load)*", q="*string(q_load))

title!("Fitted params, loading: p="*string(p_load)*", q="*string(q_load))

#savefig("../figures/fitted_2bus_M5_p3.png")

#savefig("../figures/kundur_2bus_M1_p2.png")