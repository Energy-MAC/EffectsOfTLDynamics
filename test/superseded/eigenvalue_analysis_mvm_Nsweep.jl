cd(@__DIR__)
using PowerSystems
using PowerSimulationsDynamics
using Plots
using Revise
using LaTeXStrings
using EffectsOfTLDynamics

const ETL = EffectsOfTLDynamics
const PSID = PowerSimulationsDynamics
const PSY = PowerSystems

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

file_name = "../data/json_data/mach_v_mach.json"; # choose system 
line_dict = default_2_bus_line_dict

impedance_csv = "../data/cable_data/dommel_data.csv"
capacitance_csv = "../data/cable_data/dommel_data_C.csv"

perturbation_type = "CRC"
t_fault = 0.25
perturbation_params = get_default_perturbation(t_fault, perturbation_type)
perturbation_params.t_fault = t_fault; 
perturbation_params.crc_params = CRCParam(DynamicGenerator, "generator-102-1", :V_ref, 0.95)


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


# Calculate SIL  
# SIL = VLL^2/ZC where ZC is in ohms, and VLL is in kV 
# Typical values for SIL - 14 - 140 MW for voltages up to 230 kV, and 
# 250 kV - 300-400 MW
# 500 kV - 850-1000 MW 

Vnom = 230; # kV 
# Zc in units of ohms (p.u?)
SIL = (Vnom^2)/Z_c_5_to_1_abs_3 # P in p.u., MW 

p_load = 1.0;  # base for load is ? 100?
q_load = 0.25;

load_scale = 1.0
line_scale = 1.0

l_seg = 10; # 
l = 100;

p1 = ETL.ExpParams(
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

p3 = ETL.ExpParams(
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

N = 20;

load_scales = [1,1.5, 2];
lRange = [[1075:2:1120],[1050:2:1150], [950:2:1020]];

load_scales = [0];
lRange = [[200:100:1200]];

load_bus = "BUS 2";

alg_lims = [];
dyn_lims = [];
mssb_lims = [];
msmb_lims = [];
mssb_results = [];
msmb_results = [];
dyn_results = [];
alg_results = [];
idx = 1;
for load_scale in load_scales;
    # mssb_results = [];
    # msmb_results = [];
    # dyn_results = [];
    # alg_results = [];
    p1.load_scale = load_scale
    p3.load_scale = load_scale
    current_lRange = lRange[idx][1];
    for l = current_lRange;
        line_dict["BUS 1-BUS 2-i_1"] = l;
        line_dict["BUS 1-BUS 2-i_1_static"] = l;

        p1.l = l; # update params 
        p1.load_scale = load_scale

        sim_alg = ETL.build_2bus_sim_from_file(file_name, false, false, p1, load_bus);
        save_max_nonzero_eig!(sim_alg, alg_results)

        # Dynamic pi  
        sim_dyn = ETL.build_2bus_sim_from_file(file_name, true, false, p1, load_bus);
        save_max_nonzero_eig!(sim_dyn, dyn_results)

        # MSSB
        # Define MSSB l_seg
        p1.l_seg = l/N; 
        sim_ms = ETL.build_2bus_sim_from_file(file_name, true, true, p1, load_bus);
        save_max_nonzero_eig!(sim_ms, mssb_results)

        # Do MSMB 
        p3.l_seg = l/N; 
        sim_msmb = ETL.build_2bus_sim_from_file(file_name, true, true, p3, load_bus);
        save_max_nonzero_eig!(sim_msmb, msmb_results)

    end

    display(alg_results)
    display(mssb_results)

    # Find first length where we lose stability - NaN indicates that we don't encounter stability limit 
    if findfirst(alg_results.>0) !== nothing;
        push!(alg_lims, current_lRange[findfirst(alg_results.>0)])
    else
        push!(alg_lims, NaN)
    end
    if findfirst(dyn_results.>0) !== nothing;
        push!(dyn_lims, current_lRange[findfirst(dyn_results.>0)])
    else
        push!(dyn_lims, NaN)
    end
    if findfirst(mssb_results.>0) !== nothing;
        push!(mssb_lims, current_lRange[findfirst(mssb_results.>0)])
    else
        push!(mssb_lims, NaN)
    end
    if findfirst(msmb_results.>0) !== nothing;
        push!(msmb_lims, current_lRange[findfirst(msmb_results.>0)])
    else
        push!(msmb_lims, NaN)
    end
    idx += 1
end

alg_lims
dyn_lims
mssb_lims
msmb_lims

plot(lRange[1], alg_results,xlabel="Line length (km)", ylabel="Max real λ != 0", label="Algebraic", seriestype=:line, size=(800,600))
plot!(lRange[1], dyn_results, seriestype=:line, label="Dynpi")
plot!(lRange[1], mssb_results, seriestype=:line, linewidth=2, label="MSSB N="*string(N))
plot!(lRange[1], msmb_results, seriestype=:line, linestyle=:dashdot, linewidth=4, label="MSMB N="*string(N))
title!("p="*string(p_load*load_scales[1])*", q="*string(q_load*load_scales[1]))

savefig("../figures/Ruth/mach_v_mach/line_length_v_eig_no_load.png")



plot(load_scales, alg_lims,xlabel="Load scale", ylabel="Line length @ stability boundary", label="Algebraic", seriestype=:line, size=(800,600))
plot!(load_scales, dyn_lims, seriestype=:line, label="Dynpi")
plot!(load_scales, mssb_lims, seriestype=:line, linestyle=:dashdot, label="MSSB N="*string(N))
plot!(load_scales, msmb_lims, seriestype=:line, linestyle=:dash, label="MSMB N="*string(N))

title!("Base p="*string(p_load)*", q="*string(q_load))


savefig("../figures/Ruth/mach_v_mach/loading_v_line_lim.png")


#### Stability margin plots 
plot(lRange, alg_results, label="Algebraic",linestyle=:dash, legend = :outertopright, size=(800,600),xlabel="Line length (km)", ylabel="Largest real λ, != 0")
plot!(lRange, dyn_results, label="Dynpi",linestyle=:dash)
for k in 1:length(Nrange);
    display(plot!(lRange,mssb_results[:,k],label="MSSB: N="*string(Nrange[k]),linestyle=:dashdotdot, linewidth=2))
    display(plot!(lRange,msmb_results[:,k],label="MSMB: N="*string(Nrange[k]), linestyle=:dashdot,linewidth=2))
end

function get_max(itr)
    return maximum([x for x in itr if !isnan(x)])
end

plt_ub = maximum([get_max(alg_results), get_max(dyn_results), get_max(mssb_results), get_max(msmb_results)])

plot!(lRange, zeros(length(lRange)), fillrange = ones(length(lRange))*plt_ub, fillalpha = 0.1, linealpha=0, c = 1, label="Unstable")

title!("p="*string(p_load*load_scale)*", q="*string(q_load*load_scale))

savefig("../figures/Ruth/inv_v_mach/stability_margin_Nsweep_large_range_2.png")


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
l = 300;
line_dict["BUS 1-BUS 2-i_1"] = l;
line_dict["BUS 1-BUS 2-i_1_static"] = l;

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
sim_alg = build_2bus_sim_from_file(file_name, false, false, p1)
eigs_alg = small_signal_analysis(sim_alg).eigenvalues;
sim_dyn = build_2bus_sim_from_file(file_name, true, false, p1)
eigs_dyn = small_signal_analysis(sim_dyn).eigenvalues;
sim_mssb = build_2bus_sim_from_file(file_name, true, true, p1)
eigs_mssb = small_signal_analysis(sim_mssb).eigenvalues;
sim_msmb = build_2bus_sim_from_file(file_name, true, true, p3)
eigs_msmb = small_signal_analysis(sim_msmb).eigenvalues;


plot(real(eigs_alg), imag(eigs_alg), seriestype=:scatter, xlabel="Real", ylabel="Imag", label="Algebraic",legend = :outertopright)
plot!(real(eigs_dyn), imag(eigs_dyn), seriestype=:scatter,label="Dynpi")
plot!(real(eigs_mssb), imag(eigs_mssb), seriestype=:scatter,label="MSSB")
plot!(real(eigs_msmb), imag(eigs_msmb), seriestype=:scatter,label="MSMB")

xlims!(-1000,0.1)
#ylims!(-4,4)

title!("System eigs, p="*string(p_load*load_scale)*", q="*string(q_load*load_scale))

#savefig("../figures/Ruth/dommel_params/eig_comparison_1.png")

# Analyse the participation in the least stable modes?