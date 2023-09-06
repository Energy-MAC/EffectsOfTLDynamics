cd(@__DIR__)
using PowerSystems
using PowerSimulationsDynamics
using Plots
using Revise

using CSV
using DataFrames

const PSY = PowerSystems;
const PSID = PowerSimulationsDynamics;

include("../src/ExperimentStructs.jl")
include("../src/experiment_funcs.jl")

### Choose test case
file_name = "../data/json_data/inv_v_machine.json"

line_dict = default_2_bus_line_dict

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
# "BranchTrip"
perturbation_name = "CRC"

### Define simulation parameters
sim_p = SimParams(
    abstol = 1e-13,
    reltol = 1e-10,
    maxiters = Int(1e10),
    dtmax = 1e-4,
    solver = "Rodas4",
    t_max = 2.0,
)

### Extract line data from files
M = 1

z_km, y_km, Z_c_abs, z_km_ω = get_line_parameters(impedance_csv, capacitance_csv, M)

# Calculate SIL - applies only to lossless lines
# Get L from z_km_ω
l_sil = imag(z_km_ω);
c_sil = imag(y_km)
Z_c_lossless = sqrt(l_sil/c_sil) # dimension is ohms - very very small 
Vref = 230e3;
sil = (Vref^2)/Z_c_lossless

sil_mw = sil/1e6

# SIL calc with resistance.
Z_c = sqrt(z_km_ω/y_km[1])
sil = (Vref^2)/Z_c
sil_mw = sil/1e6

# convert SIL to per unit. Isn't base 100 MW? 100 MVA


# try scale this data down a bit

# z_km = z_km*0.5;
# y_km = y_km*0.5;
# z_km_ω = z_km_ω*0.5;
# Z_c_abs = abs(sqrt(z_km_ω/y_km[1]))

# r_km = [0.0]
# Kundur parameters for testing
Z_c = 380 # Ω
r_km = 0.05 # Ω/km
x_km = 0.488 # Ω/km
g_km = 0.0 # S/km
b_km = 3.371e-6 # S/km
z_km = r_km + im*x_km
y_km = im*b_km
z_km_ω = z_km[1]

sil_kundur_mw = (Vref^2)/(Z_c*1e6)
sil_kundur_lossy_mw = (Vref^2)/((sqrt(z_km_ω/y_km))*1e6)

### Define more data
l = 100 #, 500, 750, 100 #km
# line_dict["BUS 1-BUS 2-i_1"] = l
# line_dict["BUS 1-BUS 2-i_1_static"] = l

N = nothing
t_fault = 0.25

### Get perturbation struct
perturbation_params = get_default_perturbation(t_fault, perturbation_name)
perturbation_params.crc_params = CRCParam(DynamicInverter, "generator-102-1", :V_ref, 0.95)

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
    perturbation_name, 
    perturbation_params,
    nothing,
    nothing
)

p.p_load = 0
p.q_load = 3

# Verify impedance values of raw file vs CSV data
verifying(file_name, M, impedance_csv, capacitance_csv, p)

line_model_1 = "Algebraic"
results_alg, sys, sim = run_experiment(file_name, line_model_1, p);
s = small_signal_analysis(sim)
s.eigenvalues

line_model_2 = "Dynamic"
sys = System(joinpath(pwd(), file_name));
sys_dyn = build_new_impedance_model!(sys, p, true, "")
perturbation = choose_disturbance(sys_dyn, p.perturbation, p)
sim_dyn = build_sim(sys_dyn, (0.0,2.0), perturbation, true, p)
execute_sim!(sim_dyn, p)
results_dyn = results_sim(sim_dyn)

s_dyn = small_signal_analysis(sim_dyn)
s_dyn.eigenvalues

vr_alg = get_state_series(results_alg, ("generator-3-1", :vr_filter));
vr_dyn = get_state_series(results_dyn, ("generator-3-1", :vr_filter));

plot(vr_alg, xlabel = "time", ylabel = "vr p.u.", label = "vr")
plot!(vr_dyn, xlabel = "time", ylabel = "vr p.u.", label = "vr_dyn")
plot!(title = "Line length = "*string(p.l)*" km, perturbation = "*perturbation)

line_model_3 = "Multi-Segment Dynamic"
p.N = 10

results_ms_dyn, seg_sim, seg_sys, s_seg = nothing, nothing, nothing, nothing
results_ms_dyn, seg_sys, seg_sim = run_experiment(file_name, line_model_3, p);
s_seg = small_signal_analysis(seg_sim)
s_seg.eigenvalues

vr_ms_dyn = get_state_series(results_ms_dyn, ("generator-3-1", :vr_filter));
plot(vr_ms_dyn)


plot!(xlims=(0.25, 0.3))
# plot!(ylims=(0.981,0.983))
# plot!(legend = false)
# plot!(legend=:bottomright)

l = get_component(Line, sys, "BUS 1-BUS 2-i_1")
p_branch = get_activepower_branch_flow(results_alg, "BUS 1-BUS 2-i_1", :from)
p_branch_dyn = get_activepower_branch_flow(results_dyn, "BUS 1-BUS 2-i_1", :from)
p_branch_seg = get_activepower_branch_flow(results_ms_dyn, "BUS 1-BUS 2-i_1_segment_1", :from)

plot(p_branch)
plot!(p_branch_dyn)
plot!(ylim = (-1, -0.5))

q_branch = get_reactivepower_branch_flow(results_alg, "BUS 1-BUS 2-i_1", :from)
q_branch_dyn = get_reactivepower_branch_flow(results_dyn, "BUS 1-BUS 2-i_1", :from)
plot(q_branch)
plot!(q_branch_dyn)