cd(@__DIR__)
using PowerSystems
using PowerSimulationsDynamics
using Plots
using Revise
using EffectsOfTLDynamics
using LaTeXStrings
using DataFrames
using StatsPlots

const ETL = EffectsOfTLDynamics;
const PSY = PowerSystems;
const PSID = PowerSimulationsDynamics;

# define buses 

gfm_bus = "generator-3-1";
gfl_bus = "generator-2-1";
sm_bus = "generator-1-1";
file_name = "../data/json_data/9bus_VSM3_SM1_GFL2.json"
include("../test/eig_analysis_funcs.jl")


include("9bus_defaults.jl")
sys = get_system(file_name)
p1 = p1_9bus;
p3 = p3_9bus;

gfm = get_component(Generator, sys, gfm_bus);
gfl = get_component(Generator, sys, gfl_bus);
sm = get_component(Generator, sys, sm_bus);

#gfm_p = inv_share*gfm_share; # Create p setpoints 
#gfl_p = inv_share*(1.0-gfm_share); #
load_scale = 1.0;

gfm.rating
gfm.base_power
gfm.active_power # is this wrt the gfm's base? so really 0.85*100 = 85 MW ?
gfm.active_power*gfm.base_power
gfm.reactive_power 
gfm.dynamic_injector.filter.rf
#gfm.dynamic_injector.filter.rf = 0.03
#gfm.active_power = 0.4 
gfm.active_power = 0.0


# GFL 
gfl.rating
gfl.base_power
gfl.active_power

gfl.reactive_power
gfl.active_power*gfl.base_power
#gfl.dynamic_injector.filter.cf = 0.2
gfl.active_power = 0.0
gfl.active_power*gfl.base_power

#gfl.dynamic_injector.filter.rf
#gfl.dynamic_injector.filter.rf = 0.03

# sm.rating
# sm.base_power
# sm.active_power
# sm.active_power*sm.base_power
# sm.reactive_power
#sm.active_power = 100.0


p1.load_scale = load_scale
p3.load_scale = load_scale
# Scale loads according to load scale 
for l in get_components(PSY.StandardLoad, sys)
    transform_load_to_constant_impedance(l)
    l.impedance_active_power = l.impedance_active_power * p1.load_scale 
    l.impedance_reactive_power = l.impedance_reactive_power * p1.load_scale 
end

dyn_lines = false
multi_seg = false
fd = false

if multi_seg == true;
    if fd == true
        sys = ETL.build_seg_model!(sys, p3, dyn_lines, "")
    else
        sys = ETL.build_seg_model!(sys, p1, dyn_lines, "")
    end
else
    sys = ETL.build_new_impedance_model!(sys, p1, dyn_lines, "")
end
 
perturbation = ETL.choose_disturbance(sys, p1.perturbation, p1)
tspan = (0.0, 0.25)
if fd == true;
    sim = build_sim(sys, tspan, perturbation, dyn_lines, p3)
else
    sim = build_sim(sys, tspan, perturbation, dyn_lines, p1)
end
save_max_nonzero_eig!(sim, [])[1]
ss = small_signal_analysis(sim)
summary_eigenvalues(ss)



gfm_eta, gfl_eta, sm_eta = calculate_etas(sim.sys, gfm_bus, gfl_bus, sm_bus)

inv_eta = gfm_eta + gfl_eta




gfm_after = get_component(Generator, sim.sys, gfm_bus)
gfm_after.active_power*gfm_after.base_power
#gfm_after.reactive_power
gfl_after = get_component(Generator, sim.sys, gfl_bus)
gfl_after.active_power*gfl_after.base_power
#gfl_after.reactive_power
sm_after = get_component(Generator, sim.sys, sm_bus)
sm_after.active_power*sm_after.base_power
#sm_after.reactive_power

gfm_after.active_power + gfl_after.active_power + sm_after.active_power
gfm_after.reactive_power + gfl_after.reactive_power + sm_after.reactive_power

