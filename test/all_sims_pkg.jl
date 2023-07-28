cd(@__DIR__)
using PowerSystems
using PowerSimulationsDynamics
using Plots
using EffectsOfTLDynamics

const PSY = PowerSystems;
const PSID = PowerSimulationsDynamics;

file_name = "test_sys.json"
t_max = 2.0

# "CRC"
# "NetworkSwitch"
# "InfBusChange"
dist = "CRC"

Z_c = 380 # Ω
r_km = 0.05 # Ω/km
x_km = 0.488 # Ω/km
g_km = 0 # S/km
b_km = 3.371e-6 # S/km
l = 1000 #km
N = 1
abstol = 1e-13
reltol = 1e-10
maxiters = Int(1e10)
dtmax = 1e-4
sim_params = SimParams(abstol, reltol, maxiters, 1e-4, "Rodas4", t_max)
# TODO: We can't build the perturbation without initializing the system.
perturbation = "CRC"
t_fault = 0.25
perturbation_params = PerturbationParams(0.25, nothing, nothing, nothing, nothing, nothing, nothing)
p = ExpParams(N, l, Z_c, r_km, x_km, g_km, b_km, sim_params, perturbation, perturbation_params)

line_model_1 = "Algebraic"
results_alg, sys = run_experiment(file_name, t_max, dist, line_model_1, p);

line_model_2 = "Dynamic"
results_dyn, dyn_sys = run_experiment(file_name, t_max, dist, line_model_2, p);

vr_alg = get_state_series(results_alg, ("generator-102-1", :vr_filter));
# sim_time, vr_alg = get_voltage_magnitude_series(results_alg, 1)
plot(vr_alg, xlabel = "time", ylabel = "vr p.u.", label = "vr")
vr_dyn = get_state_series(results_dyn, ("generator-102-1", :vr_filter));
plot!(vr_dyn, xlabel = "time", ylabel = "vr p.u.", label = "vr_dyn")

line_model_3 = "Multi-Segment Dynamic"

for N in [1]
    print(N)
    p = ExpParams(N, l, Z_c, r_km, x_km, g_km, b_km, abstol, reltol, maxiters)
    
    results_ms_dyn, seg_sys = nothing, nothing
    results_ms_dyn, seg_sys = run_experiment(file_name, t_max, dist, line_model_3, p)

    vr_ms_dyn = get_state_series(results_ms_dyn, ("generator-102-1", :vr_filter));
    display(plot!(vr_ms_dyn, xlabel = "time", ylabel = "vr p.u.", label = "vr_segs_$(N)"))
end

plot!(xlims=(0, 1))
