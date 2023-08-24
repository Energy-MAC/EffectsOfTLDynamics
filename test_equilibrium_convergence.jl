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

include("../src/ExperimentStructs.jl")
include("../src/experiment_funcs.jl")

### Choose test case
# "SIIB.json"
# "9bus.json"
# "inv_v_machine.json"
# "twobus_2inv.json"
# "9bus_slackless.json"
file_name = "../data/json_data/inv_v_machine.json"

### Load relevant line data
impedance_csv = "../data/cable_data/impedance_data.csv"
capacitance_csv = "../data/cable_data/C_per_km.csv"

### Choose perturbation to be applied
# "BIC"
# "GenTrip"
# "CRC"
# "LoadChange"
# "LoadTrip"
# "InfBusChange"
perturbation = "CRC"

### Define simulation parameters
sim_p = SimParams(
    abstol = 1e13,
    reltol = 1e10,
    maxiters = Int(1e10),
    dtmax = 1e-4,
    solver = "Rodas4",
    t_max = 2.0,
)

### Extract line data from files
M = 1

#r_km, x_km, g_km, b_km, Z_c, r_km_pi, x_km_pi, Z_c_pi = get_line_parameters(impedance_csv, capacitance_csv, M)

#Kundur parameters for testing
Z_c = 380 # Ω
r_km = 0.05 # Ω/km
x_km = 0.488 # Ω/km
g_km = 0.0 # S/km
b_km = 3.371e-6 # S/km

r_km_pi = r_km;
x_km_pi = x_km;
Z_c_pi = Z_c;

r_km = [r_km]
x_km = [x_km]
g_km = [g_km]
b_km = [b_km]

### Define more data
l = 500 #km
N = nothing
t_fault = 0.25

### Get perturbation struct
perturbation_params = get_default_perturbation(t_fault, perturbation)
# perturbation_params.crc_params = CRCParam(DynamicInverter, "generator-1-1", :V_ref, 0.95)

p = ExpParams(
    N, 
    M, 
    l, 
    Z_c, 
    r_km, 
    x_km, 
    g_km, 
    b_km, 
    r_km_pi,
    x_km_pi,
    Z_c_pi,
    sim_p, 
    perturbation, 
    perturbation_params)

line_model_1 = "Algebraic"
results_alg, sys, alg_sim = run_experiment(file_name, line_model_1, p);

line_model_2 = "Dynamic"
results_dyn, sys_dyn, dyn_sim = run_experiment(file_name, line_model_2, p);

vr_alg = get_state_series(results_alg, ("generator-102-1", :vr_filter));
vr_dyn = get_state_series(results_dyn, ("generator-102-1", :vr_filter));

plot(vr_alg, xlabel = "time", ylabel = "vr p.u.", label = "vr")
plot!(vr_dyn, xlabel = "time", ylabel = "vr p.u.", label = "vr_dyn")

line_model_3 = "Multi-Segment Dynamic"

for n in [20]
    print(n)
    p.N = n
    results_ms_dyn, seg_sys = nothing, nothing
    results_ms_dyn, seg_sys = run_experiment(file_name, line_model_3, p);

    vr_ms_dyn = get_state_series(results_ms_dyn, ("generator-102-1", :vr_filter));
    display(plot!(vr_ms_dyn, xlabel = "time", ylabel = "vr p.u.", label = "vr_segs_$(p.N)_branch_$(p.M)"))
end

plot!(xlims=(0.24, 0.275))
plot!(ylims=(0.9814,0.909))
plot!(legend = true)
plot!(legend=:bottomright)

l = first(get_components(Line, sys))

l_dyn = first(get_components(Line, sys_dyn))

n = 1
p.N = n
results_ms_dyn, sys_seg = nothing, nothing
results_ms_dyn, sys_seg = run_experiment(file_name, line_model_3, p);

y_total = 0;

b_total = 0
for l in get_components(Line, sys)
    r = l.r;
    x = l.x;
    z = r + im*x;
    b = l.b.from
    b_total = b_total + b
    y_total = y_total + 1/z;    
end

z_total = 1/y_total;
r_total = real(z_total)
x_total = imag(z_total)
b_total

l.r
l_dyn.r
l_seg.r

l.x
l_dyn.x
l_seg.x

l.b.from
l.b




# Extract l and r from initial states

# Need to get voltage states at t=0, and line current states
ilr = get_state_series(results_dyn, "Il_R")

ilr_dyn = get_state_series(results_dyn, ("BUS 1-BUS 2-i_1", :Il_R))[2][1] # 2 is values, 1 is first timestep
ili_dyn = get_state_series(results_dyn, ("BUS 1-BUS 2-i_1", :Il_I))[2][1]
v1r_dyn = get_state_series(results_dyn, ("V_1", :R))[2][1]
v1i_dyn = get_state_series(results_dyn, ("V_1", :I))[2][1]
v2r_dyn = get_state_series(results_dyn, ("V_2", :R))[2][1]
v2i_dyn = get_state_series(results_dyn, ("V_2", :I))[2][1]

# Solve for r, x
# v1r - v2r - Rilr + omega*L*ili = 0 
# v1i - v2i - Rili - omega*L*ilr = 0

# R, L
A = [ilr_dyn -ili_dyn; ili_dyn ilr_dyn];
b = [v1r_dyn-v2r_dyn; v1i_dyn-v1i_dyn];

x = inv(A)*b
r_line = x[1]
x_line = x[2]

# Get our line
line = get_component(Line, sys_dyn, "BUS 1-BUS 2-i_1")


show_states_initial_value(dyn_sim)