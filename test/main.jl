cd(@__DIR__)
using PowerSystems
using PowerSimulationsDynamics
using Plots
using Revise
using EffectsOfTLDynamics

const PSY = PowerSystems;
const PSID = PowerSimulationsDynamics;

sim_p = SimParams(
    abstol = 1e13,
    reltol = 1e10,
    maxiters = Int(1e10),
    dtmax = 1e-4,
    solver = "Rodas4",
    t_max = 2.0,
)

file_name = "OMIB.json"

# "BIC"
# "GenTrip"
# "CRC"
# "LoadChange"
# "LoadTrip"
# "InfBusChange"
perturbation = "CRC"

Z_c = 380 # Ω
r_km = 0.05 # Ω/km
x_km = 0.488 # Ω/km
g_km = 0 # S/km
b_km = 3.371e-6 # S/km
l = 100 #km
N = nothing
t_fault = 0.25
perturbation_params = get_default_perturbation(t_fault, perturbation)
p = ExpParams(N, l, Z_c, r_km, x_km, g_km, b_km, sim_p, perturbation, perturbation_params)

line_model_1 = "Algebraic"
results_alg, sys = run_experiment(file_name, line_model_1, p);

line_model_2 = "Dynamic"
results_dyn, dyn_sys = run_experiment(file_name, line_model_2, p);

vr_alg = get_state_series(results_alg, ("generator-102-1", :vr_filter));
# sim_time, vr_alg = get_voltage_magnitude_series(results_alg, 1)
plot(vr_alg, xlabel = "time", ylabel = "vr p.u.", label = "vr")
vr_dyn = get_state_series(results_dyn, ("generator-102-1", :vr_filter));
plot!(vr_dyn, xlabel = "time", ylabel = "vr p.u.", label = "vr_dyn")

line_model_3 = "Multi-Segment Dynamic"


for n in [1, 2, 3, 4, 5, 10, 15]
    print(n)
    p.N = n
    results_ms_dyn, seg_sys = nothing, nothing
    results_ms_dyn, seg_sys = run_experiment(file_name, line_model_3, p)

    vr_ms_dyn = get_state_series(results_ms_dyn, ("generator-102-1", :vr_filter));
    display(plot!(vr_ms_dyn, xlabel = "time", ylabel = "vr p.u.", label = "vr_segs_$(p.N)"))
end

plot!(xlims=(0.249, 0.255))
plot!(ylims=(0.98, 0.99))
