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


function build_sim_from_file(file_name::String, dyn_lines::Bool, multi_segment::Bool, p::ExpParams, load_bus::String, inv_share, rating_scale)

    # build system
    sys = System(joinpath(pwd(), file_name));

    # Simulation time span
    tspan = (0.0, p.sim_params.t_max);
    perturbation = ETL.choose_disturbance(sys, p.perturbation, p);

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

    for g in get_components(PSY.Generator, sys)
        set_rating!(g, g.rating * rating_scale) # adjust the ratings of all generators
        set_active_power!(g, g.rating * inv_share); # adjust the share of all generators 
        # if g.bus.name == inv_bus # at the inverter bus 
        #     set_active_power!(g, g.rating * inv_share); # adjust P setpoint to be correct fraction of rating 
        # end
        
    end
    
    # build segments model
    if (multi_segment == true)
        sys = ETL.build_seg_model!(sys, p, dyn_lines, "");
    else
        sys = ETL.build_new_impedance_model!(sys, p, dyn_lines, "");
    end

    sim = build_sim(sys, tspan, perturbation, dyn_lines, p);
    return sim 

end

include("eig_analysis_funcs.jl")

line_dict = default_2_bus_line_dict
foldername = "../figures/Ruth/current/twobus";

impedance_csv = "../data/cable_data/dommel_data.csv"
capacitance_csv = "../data/cable_data/dommel_data_C.csv"

# Perturbation type not relevant for stability analysis 
perturbation_type = "CRC"
t_fault = 0.25
perturbation_params = get_default_perturbation(t_fault, perturbation_type)

### Extract line data from files for M=1 and M=3
M = 1
factor_z = 1.0 # to get it closer to kundur 
factor_y = 1.0
z_km_1, y_km_1, Z_c_abs_1, z_km_ω_1, z_km_ω_5_to_1_1, Z_c_5_to_1_abs_1 = get_line_parameters(impedance_csv, capacitance_csv, M, factor_z, factor_y)

M = 3
z_km_3, y_km_3, Z_c_abs_3, z_km_ω_3, z_km_ω_5_to_1_3, Z_c_5_to_1_abs_3 = get_line_parameters(impedance_csv, capacitance_csv, M, factor_z, factor_y)

# Calculate SIL 
V_nom = 230 # kV
Z_o = sqrt(z_km_ω_5_to_1_1/y_km_1)
SIL = V_nom ^2/Z_o
p_load = real(SIL)/100
q_load = imag(SIL)/100

l_seg = 10.0; # km 

### Define simulation parameters
sim_p = SimParams(
    abstol = 1e-13,
    reltol = 1e-10,
    maxiters = Int(1e10),
    dtmax = 1e-4,
    solver = "Rodas4",
    t_max = 0.1,
)

p1 = ETL.ExpParams(
    nothing, # N 
    1, # M
    100, # base line length 
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
    1.0, # line_scale
    1.0 # load_scale
);

p3 = ETL.ExpParams(
    nothing, # N
    3, # M
    100, # base line length 
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
    1.0, # line_scale
    1.0 # load_scale
);

##------------------##

# inv_v_machine
# twobus_2inv
# mach_v_mach

file_name = "../data/json_data/inv_v_machine.json"; 
load_bus = "BUS 1"; # this does make a small difference to inverter penetration 

p1.line_scale = 5.0;
p3.line_scale = p1.line_scale;
p1.load_scale = 1.0;
p3.load_scale = p1.load_scale;

gen_scale_range = 0.1:0.1:3.0

sim = build_sim_from_file(file_name, false, false, p1, load_bus, 0.5, 1.0);
inv = get_component(Generator, sim.sys,"generator-102-1");
sm = get_component(Generator, sim.sys, "generator-101-1");

# Algebraic 
etas = [];
eigs = [];
for gen_scale = gen_scale_range;
    sim = build_sim_from_file(file_name, false, false, p1, load_bus, 0.5, gen_scale);
    inv = get_component(Generator, sim.sys,"generator-102-1");
    sm = get_component(Generator, sim.sys, "generator-101-1");
    eta = inv.active_power/(inv.active_power + sm.active_power);
    push!(etas, eta)
    save_max_nonzero_eig!(sim, eigs);
end

plot(etas, eigs, xlabel=L"\eta\ (inv)", ylabel=L"Max\ \lambda", label="Algebraic")

# Dynamic 
etas_d = [];
eigs_d = [];
for gen_scale = gen_scale_range;
    sim_d = build_sim_from_file(file_name, true, false, p1, load_bus, 0.5, gen_scale);
    inv_d = get_component(Generator, sim_d.sys,"generator-102-1");
    sm_d = get_component(Generator, sim_d.sys, "generator-101-1");
    eta_d = inv_d.active_power/(inv_d.active_power + sm_d.active_power);
    push!(etas_d, eta_d)
    save_max_nonzero_eig!(sim_d, eigs_d);
end

ylabel = L"\mathcal{R}(\lambda)"
plot!(etas_d, eigs_d, xlabel=L"\eta", ylabel=ylabel, label="Dynamic")

# MSSB 
etas_ms = [];
eigs_ms = [];
for gen_scale = gen_scale_range;
    sim_ms = build_sim_from_file(file_name, true, true, p1, load_bus, 0.5, gen_scale);
    inv_ms = get_component(Generator, sim_ms.sys,"generator-102-1");
    sm_ms = get_component(Generator, sim_ms.sys, "generator-101-1");
    eta_ms = inv_ms.active_power/(inv_ms.active_power + sm_ms.active_power);
    push!(etas_ms, eta_ms)
    save_max_nonzero_eig!(sim_ms, eigs_ms);
end

plot!(etas_ms, eigs_ms, xlabel=L"\eta", ylabel=L"Max\ \lambda", label="MSSB")











# # Check setpoint of modified system

# inv = get_component(Generator, sim.sys,"generator-102-1");
# inv.active_power
# inv.reactive_power # this changes based on P (bc of powerflow)
# inv.rating
# sm = get_component(Generator, sim.sys, "generator-101-1");
# sm.active_power
# sm.reactive_power
# sm.rating
# eta = inv.active_power/(inv.active_power + sm.active_power)


# sim = build_sim_from_file(file_name, false, false, p1, load_bus, "BUS 2", 0.5, 2.0);
# inv = get_component(Generator, sim.sys,"generator-102-1");
# inv.active_power
# inv.reactive_power # this changes based on P (bc of powerflow)
# inv.rating
# sm = get_component(Generator, sim.sys, "generator-101-1");
# sm.active_power
# sm.reactive_power
# sm.rating
# eta = inv.active_power/(inv.active_power + sm.active_power)




# #ss_alg = small_signal_analysis(sim_alg) 
# output = [];
# save_max_nonzero_eig!(sim_alg, output);