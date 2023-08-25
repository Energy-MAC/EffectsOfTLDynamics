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

### Choose test case
# "SIIB.json"
# "9bus.json"
# "inv_v_machine.json"
# "twobus_2inv.json"
# "9bus_slackless.json"
file_name = "../data/json_data/9bus.json"
# default_2_bus_line_dict - For 2 bus system
# default_9_bus_line_dict - For 9 bus system
line_dict = default_9_bus_line_dict

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
perturbation = "InfBusChange"

### Define simulation parameters
sim_p = SimParams(
    abstol = 1e-13,
    reltol = 1e-10,
    maxiters = Int(1e10),
    dtmax = 1e-4,
    solver = "Rodas4",
    t_max = 20.0,
)

### Extract line data from files
M = 1

z_km, y_km, Z_c_abs, z_km_ω = get_line_parameters(impedance_csv, capacitance_csv, M)
# r_km = [0.0]
# Kundur parameters for testing
# Z_c = 380 # Ω
# r_km = 0.05 # Ω/km
# x_km = 0.488 # Ω/km
# g_km = 0.0 # S/km
# b_km = 3.371e-6 # S/km

### Define more data
l = 500 #, 500, 750, 100 #km
N = nothing
t_fault = 0.25

### Get perturbation struct
perturbation_params = get_default_perturbation(t_fault, perturbation)
# perturbation_params.crc_params = CRCParam(DynamicInverter, "generator-1-1", :V_ref, 0.95)

p = ExpParams(
    N, 
    M, 
    l, 
    Z_c_abs, 
    z_km,
    y_km,
    z_km_ω, 
    line_dict,
    sim_p, 
    perturbation, 
    perturbation_params)

line_model_1 = "Algebraic"
results_alg, sys = run_experiment(file_name, line_model_1, p);

line_model_2 = "Dynamic"
results_dyn, sys_dyn = run_experiment(file_name, line_model_2, p);

vr_alg = get_state_series(results_alg, ("generator-102-1", :vr_filter));
vr_dyn = get_state_series(results_dyn, ("generator-102-1", :vr_filter));

plot(vr_alg, xlabel = "time", ylabel = "vr p.u.", label = "vr")
plot!(vr_dyn, xlabel = "time", ylabel = "vr p.u.", label = "vr_dyn")
plot!(title = "Line length = "*string(p.l)*" km, perturbation = "*perturbation)

line_model_3 = "Multi-Segment Dynamic"

for n in [1,2,3,4,5,10]
    print(n)
    p.N = n
    results_ms_dyn, seg_sys = nothing, nothing
    results_ms_dyn, seg_sys = run_experiment(file_name, line_model_3, p)

    vr_ms_dyn = get_state_series(results_ms_dyn, ("generator-102-1", :vr_filter));
    display(plot!(vr_ms_dyn, xlabel = "time", ylabel = "vr p.u.", label = "vr_segs_$(p.N)_branch_$(p.M)_L"))
end

plot!(xlims=(0.24, 0.5))
plot!(ylims=(0.981,0.983))
plot!(legend = false)
plot!(legend=:bottomright)

# l = first(get_components(Line, sys))

# l_dyn = first(get_components(Line, sys_dyn))

# n = 1
# p.N = n
# results_ms_dyn, sys_seg = nothing, nothing
# results_ms_dyn, sys_seg = run_experiment(file_name, line_model_3, p);

# y_total = 0;

# b_total = 0
# for l in get_components(Line, sys)
#     r = l.r;
#     x = l.x;
#     z = r + im*x;
#     b = l.b.from
#     b_total = b_total + b
#     y_total = y_total + 1/z;    
# end

# z_total = 1/y_total;
# r_total = real(z_total)
# x_total = imag(z_total)
# b_total

# l.r
# l_dyn.r
# l_seg.r

# l.x
# l_dyn.x
# l_seg.x

# l.b.from
# l.b