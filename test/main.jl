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
# "OMIB.json"
# "9bus.json"
file_name = "OMIB.json"

# "BIC"
# "GenTrip"
# "CRC"
# "LoadChange"
# "LoadTrip"
# "InfBusChange"
perturbation = "CRC"

M = 1
# Z_c, r_km, x_km, g_km, b_km = get_line_parameters(data_file, M)

Z_c = 380 # Ω
r_km = 0.05 # Ω/km
x_km = 0.488 # Ω/km
g_km = 0 # S/km
b_km = 3.371e-6 # S/km
l = 1000 #km
N = nothing
t_fault = 0.25
perturbation_params = get_default_perturbation(t_fault, perturbation)
p = ExpParams(N, M, l, Z_c, r_km, x_km, g_km, b_km, sim_p, perturbation, perturbation_params)

line_model_1 = "Algebraic"
results_alg, sys = run_experiment(file_name, line_model_1, p);

line_model_2 = "Dynamic"
results_dyn, dyn_sys = run_experiment(file_name, line_model_2, p);

vr_alg = get_state_series(results_alg, ("generator-102-1", :vr_filter));
vr_dyn = get_state_series(results_dyn, ("generator-102-1", :vr_filter));

using PlotlyJS, LaTeXStrings

trace1 = scatter(
    x = vr_alg[1],
    y = vr_alg[2],
    mode = "lines",
)

layout = Layout(
    xaxis_title = L"$\text{time} (s)$",
    yaxis_title = L"$V_f [\text{p.u.}$"
)

plot(trace1)

plot(vr_alg, xlabel = "time", ylabel = "vr p.u.", label = "vr")
plot!(vr_dyn, xlabel = "time", ylabel = "vr p.u.", label = "vr_dyn")

line_model_3 = "Multi-Segment Dynamic"

for n in [10]
    print(n)
    p.N = n
    results_ms_dyn, seg_sys = nothing, nothing
    results_ms_dyn, seg_sys = run_experiment(file_name, line_model_3, p)

    vr_ms_dyn = get_state_series(results_ms_dyn, ("generator-102-1", :vr_filter));
    display(plot!(vr_ms_dyn, xlabel = "time", ylabel = "vr p.u.", label = "vr_segs_$(p.N)_branch_$(p.M)"))
end

plot!(xlims=(0.24, 0.5))
plot!(ylims=(0.8,0.95))
