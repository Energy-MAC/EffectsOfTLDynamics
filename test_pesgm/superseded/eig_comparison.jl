cd(@__DIR__)
using PowerSystems
using PowerSimulationsDynamics
using Plots
using Revise
using LaTeXStrings
using EffectsOfTLDynamics

using CSV
using DataFrames

const ETL = EffectsOfTLDynamics;
const PSY = PowerSystems;
const PSID = PowerSimulationsDynamics;

include("../test/eig_analysis_funcs.jl")
include("defaults.jl")


## System options ##

# inv_v_machine
# twobus_2inv
# mach_v_mach
# droop_v_GFL
# droop_v_machine
# dVOC_v_GFL
# dvoc_v_machine
# GFL_v_machine
# VSM_v_GFL

file_name = "../data/json_data/droop_v_machine.json"; 
load_bus = "BUS 1"; # this does make a small difference to results 

p1 = dommel_M1

inv_share = 0.5;
rating_scale = 1.0;
sim = build_sim_from_file(file_name, false, false, p1, load_bus, inv_share, rating_scale, false);
inv = get_component(Generator, sim.sys,"generator-102-1");
sm = get_component(Generator, sim.sys, "generator-101-1");
eta = inv.active_power/(inv.active_power + sm.active_power)
# plot eigs
ss_alg = small_signal_analysis(sim)
alg_eigs = ss_alg.eigenvalues
plot(real(alg_eigs), imag(alg_eigs), label="Algebraic", seriestype=:scatter, ms=6)
xlabel!(L"\mathrm{Re}\ [\lambda]")
ylabel!(L"\mathrm{Im}\ [\lambda]")

sim_d = build_sim_from_file(file_name, true, false, p1, load_bus, frac, rating_scale, false);
inv_d = get_component(Generator, sim_d.sys,"generator-102-1");
sm_d = get_component(Generator, sim_d.sys, "generator-101-1");
eta_d = inv_d.active_power/(inv_d.active_power + sm_d.active_power);
ss_dyn = small_signal_analysis(sim_d)
dyn_eigs = ss_dyn.eigenvalues
plot!(real(dyn_eigs), imag(dyn_eigs), label="Dynamic", seriestype=:scatter, ms=4)

sim_ms = build_sim_from_file(file_name, true, true, p1, load_bus, frac, rating_scale, false);
inv_ms = get_component(Generator, sim_ms.sys,"generator-102-1");
sm_ms = get_component(Generator, sim_ms.sys, "generator-101-1");
eta_ms = inv_ms.active_power/(inv_ms.active_power + sm_ms.active_power);
ss_ms = small_signal_analysis(sim_ms)
ms_eigs = ss_ms.eigenvalues
plot!(real(ms_eigs), imag(ms_eigs), label="MSSB", seriestype=:scatter, ms=2)

xlims!(-5000, 100)

system_name = split(file_name, "/")[end][1:end-5]
title!(system_name*", load @ "*load_bus)

savefig("fig_results/eig_figs/"*system_name*".png")


