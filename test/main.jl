cd(@__DIR__)
using PowerSystems
using PowerSimulationsDynamics
using Plots
using Revise
using EffectsOfTLDynamics

using CSV
using DataFrames

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

# "SIIB.json"
# "9bus.json"
# "inv_v_machine.json"
# "twobus_2inv.json"
# "9bus_slackless.json"
file_name = "../data/json_data/SIIB.json"

# "BIC"
# "GenTrip"
# "CRC"
# "LoadChange"
# "LoadTrip"
# "InfBusChange"
perturbation = "InfBusChange"

M = 2
impedance_csv = "../data/cable_data/impedance_data.csv"
capacitance_csv = "../data/cable_data/C_per_km.csv"
r_km, l_km, c_km = get_line_parameters(impedance_csv, capacitance_csv, M)
# Z_c, r_km, x_km, g_km, b_km = get_line_parameters(impedance_csv, capacitance_csv, M)

Z_c = 380 # Ω
r_km = 0.05 # Ω/km
x_km = 0.488 # Ω/km
g_km = 0 # S/km
b_km = 3.371e-6 # S/km
l = 100 #km
N = nothing
t_fault = 0.25
perturbation_params = get_default_perturbation(t_fault, perturbation)
# perturbation_params.crc_params = CRCParam(DynamicInverter, "generator-1-1", :V_ref, 0.95)
p = ExpParams(N, M, l, Z_c, r_km, x_km, g_km, b_km, sim_p, perturbation, perturbation_params)

line_model_1 = "Algebraic"
results_alg, sys = run_experiment(file_name, line_model_1, p);

line_model_2 = "Dynamic"
results_dyn, dyn_sys = run_experiment(file_name, line_model_2, p);

vr_alg = get_state_series(results_alg, ("generator-102-1", :vr_filter));
vr_dyn = get_state_series(results_dyn, ("generator-102-1", :vr_filter));

plot(vr_alg, xlabel = "time", ylabel = "vr p.u.", label = "vr")
plot!(vr_dyn, xlabel = "time", ylabel = "vr p.u.", label = "vr_dyn")

line_model_3 = "Multi-Segment Dynamic"

for n in [2, 3, 4, 5, 10]
    print(n)
    p.N = n
    results_ms_dyn, seg_sys = nothing, nothing
    results_ms_dyn, seg_sys = run_experiment(file_name, line_model_3, p)

    vr_ms_dyn = get_state_series(results_ms_dyn, ("generator-102-1", :vr_filter));
    display(plot!(vr_ms_dyn, xlabel = "time", ylabel = "vr p.u.", label = "vr_segs_$(p.N)_branch_$(p.M)"))
end

plot!(xlims=(0.24, 0.5))
plot!(ylims=(0.8,1))
plot!(legend = true)
plot!(legend=:bottomright)
