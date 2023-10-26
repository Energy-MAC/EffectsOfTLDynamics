cd(@__DIR__)
using PowerSystems
using PowerSimulationsDynamics
using Plots
using Revise
using LaTeXStrings
using EffectsOfTLDynamics

const ETL = EffectsOfTLDynamics;
const PSY = PowerSystems;
const PSID = PowerSimulationsDynamics;

include("defaults.jl")
include("../test/eig_analysis_funcs.jl")

##------------------##
## Systems ##

# inv_v_machine
# twobus_2inv
# mach_v_mach
# droop_v_machine
# dvoc_v_machine
# GFL_v_machine
# VSM_v_GFL
# VSM_v_dVOC
# VSM_V_droop


p1 = dommel_M1 # default dommel parameters, with no line or load scaling, SIL load. 

file_name = "../data/json_data/gorno.json"; 
update_ref_bus = false;

inv_bus = "generator-102-1"; # bus with generator we want to adjust share of 
non_inv_bus = "generator-101-1"; # other bus 

load_bus = "BUS 2"; # this does make a small difference to results 

gen_scale_range = 0.05:0.1:3.0

# Algebraic 
etas = [];
eigs = [];
for gen_scale = gen_scale_range;
    sim = build_sim_from_file(file_name, false, false, p1, load_bus, 0.5, gen_scale, update_ref_bus);
    inv = get_component(Generator, sim.sys, inv_bus);
    sm = get_component(Generator, sim.sys, non_inv_bus);
    eta = inv.active_power/(inv.active_power + sm.active_power);
    push!(etas, eta)
    save_max_nonzero_eig!(sim, eigs);
end

plot(etas, eigs, xlabel=L"\eta\ (inv)", ylabel=L"\mathcal{R}(\lambda)", label="Algebraic")

# Dynamic 
etas_d = [];
eigs_d = [];
for gen_scale = gen_scale_range;
    sim_d = build_sim_from_file(file_name, true, false, p1, load_bus, 0.5, gen_scale,update_ref_bus);
    inv_d = get_component(Generator, sim_d.sys, inv_bus);
    sm_d = get_component(Generator, sim_d.sys, non_inv_bus);
    eta_d = inv_d.active_power/(inv_d.active_power + sm_d.active_power);
    push!(etas_d, eta_d)
    save_max_nonzero_eig!(sim_d, eigs_d);
end

plot!(etas_d, eigs_d, label="Dynamic")


# # MSSB 
# etas_ms = [];
# eigs_ms = [];
# for gen_scale = gen_scale_range;
#     sim_ms = build_sim_from_file(file_name, true, true, p1, load_bus, 0.5, gen_scale, update_ref_bus);
#     inv_ms = get_component(Generator, sim_ms.sys,"generator-102-1");
#     sm_ms = get_component(Generator, sim_ms.sys, "generator-101-1");
#     eta_ms = inv_ms.active_power/(inv_ms.active_power + sm_ms.active_power);
#     push!(etas_ms, eta_ms)
#     save_max_nonzero_eig!(sim_ms, eigs_ms);
# end

# plot!(etas_ms, eigs_ms, label="MSSB")

# p3 = dommel_M3;

# # MSMB
# etas_mb = [];
# eigs_mb = [];
# for gen_scale = gen_scale_range;
#     sim_mb = build_sim_from_file(file_name, true, true, p3, load_bus, 0.5, gen_scale, update_ref_bus);
#     inv_mb = get_component(Generator, sim_mb.sys,"generator-102-1");
#     sm_mb = get_component(Generator, sim_mb.sys, "generator-101-1");
#     eta_mb = inv_mb.active_power/(inv_mb.active_power + sm_mb.active_power);
#     push!(etas_mb, eta_mb)
#     save_max_nonzero_eig!(sim_mb, eigs_mb);
# end

# plot!(etas_mb, eigs_mb, label="MSMB")



system_name = split(file_name, "/")[end][1:end-5]
title!(system_name*", load at "*load_bus)

savefig("fig_results/initial_figs/"*system_name*".png")

