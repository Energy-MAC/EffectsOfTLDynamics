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

file_name = "../data/json_data/VSM_v_GFL.json"; 
load_bus = "BUS 1"; # this does make a small difference to results 
p = dommel_M1

sys = System(joinpath(pwd(), file_name));

#vsm = get_component(Generator, sys, "generator-102-1")
#gfl = get_component(Generator, sys, "generator-101-1")

# Update which bus is the ref bus 
#b1 = get_component(Bus, sys, "BUS 1")
#b1.bustype = BusTypes.PV

#b2 = get_component(Bus, sys, "BUS 2")
#b2.bustype = BusTypes.REF

# Simulation time span
tspan = (0.0, p.sim_params.t_max);
perturbation = choose_disturbance(sys, p.perturbation, p);

# Add load
load = StandardLoad(
    name = "load1",
    available = true,
    bus = get_component(Bus, sys, load_bus),
    base_power = 100.0,
    constant_active_power = 0.0,
    constant_reactive_power = 0.0,
    impedance_active_power = p.p_load*p.load_scale,
    impedance_reactive_power = p.q_load*p.load_scale,
    current_active_power = 0.0,
    current_reactive_power = 0.0,
    max_constant_active_power = 0.0,
    max_constant_reactive_power = 0.0,
    max_impedance_active_power = p.p_load*p.load_scale,
    max_impedance_reactive_power = p.q_load*p.load_scale,
    max_current_active_power = 0.0,
    max_current_reactive_power = 0.0,
);
add_component!(sys, load);

# Try build algebraic 
sys = ETL.build_new_impedance_model!(sys, p1, false, "");
sim = build_sim(sys, tspan, perturbation, false, p1);
