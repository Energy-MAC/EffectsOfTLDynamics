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

file_name = "../data/json_data/inv_v_machine.json"; 
update_ref_bus = false;

inv_bus = "generator-102-1"; # bus with generator we want to adjust share of 
non_inv_bus = "generator-101-1"; # other bus 

load_bus = "BUS 2"; # this does make a small difference to results 

gen_scale_range = 0.05:0.1:3.0

# Algebraic 
etas = [];
eigs = [];
gen_scale = 1.0;

sys = System(joinpath(pwd(), file_name));

vsm = get_component(Generator, sys, "generator-102-1")



perturbation = ETL.choose_disturbance(sys, p1.perturbation, p1);
# Add load
load = StandardLoad(
    name = "load1",
    available = true,
    bus = get_component(Bus, sys, load_bus),
    base_power = 100.0,
    constant_active_power = 0.0,
    constant_reactive_power = 0.0,
    impedance_active_power = p1.p_load*p1.load_scale,
    impedance_reactive_power = p1.q_load*p1.load_scale,
    current_active_power = 0.0,
    current_reactive_power = 0.0,
    max_constant_active_power = 0.0,
    max_constant_reactive_power = 0.0,
    max_impedance_active_power = p1.p_load*p1.load_scale,
    max_impedance_reactive_power = p1.q_load*p1.load_scale,
    max_current_active_power = 0.0,
    max_current_reactive_power = 0.0,
);
add_component!(sys, load);

for g in get_components(PSY.Generator, sys)
    set_rating!(g, g.rating * gen_scale) # adjust the ratings of all generators
    set_active_power!(g, g.rating * 0.5); # adjust the share of all generators 
end
sys = ETL.build_new_impedance_model!(sys, p1, false, "");

sim = build_sim(sys, tspan, perturbation, false, p1);

inv = get_component(Generator, sim.sys, inv_bus);
sm = get_component(Generator, sim.sys, non_inv_bus);
eta = inv.active_power/(inv.active_power + sm.active_power);
push!(etas, eta)
save_max_nonzero_eig!(sim, eigs);







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


