cd(@__DIR__)
using PowerSystems
using PowerSimulationsDynamics
using Plots
using Revise
using EffectsOfTLDynamics

const ETL = EffectsOfTLDynamics;
const PSY = PowerSystems;
const PSID = PowerSimulationsDynamics;

file_name = "../data/json_data/9bus_allgens.json"
include("../test/eig_analysis_funcs.jl")
include("9bus_defaults.jl")
sm_bus = "generator-3-1"
sys = get_system(file_name)
p1 = p1_9bus;
p1.line_scale = 1.0;
p1.load_scale = 1.0;

total_p_load = 0.0;
for l in get_components(PSY.StandardLoad, sys)
    transform_load_to_constant_impedance(l)
    l.impedance_active_power = l.impedance_active_power * p1.load_scale 
    l.impedance_reactive_power = l.impedance_reactive_power * p1.load_scale 
    total_p_load += l.impedance_active_power;
end

total_p_load = total_p_load*sys.units_settings.base_value

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


p1.perturbation_params.crc_params = CRCParam(DynamicGenerator, sm_bus, :V_ref, 0.95)
perturbation = choose_disturbance(sys, p1.perturbation, p1)

tspan = (0.0, 0.25)
if fd == true;
    sim = build_sim(sys, tspan, perturbation, dyn_lines, p3)
else
    sim = build_sim(sys, tspan, perturbation, dyn_lines, p1)
end
ss = small_signal_analysis(sim)
all_gens = 


res = [];
for load_scale = collect(0.0:0.1:1.0)
    ss = check_SM(load_scale)
    push!(res, real(ss.eigenvalues[end-2:end]))

end
res
