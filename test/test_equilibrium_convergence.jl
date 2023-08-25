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
    abstol = 1e-13,
    reltol = 1e-10,
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
#r_km = 0.0;
x_km = 0.488 # Ω/km
g_km = 0.0 # S/km
b_km = 3.371e-6 # S/km

# Calculate Z_c 
# z_km = r_km + im*x_km;
# y_km = im*b_km;
# Z_c = abs(sqrt(z_km[1]/y_km[1]))

r_km_pi = r_km;
x_km_pi = x_km;
Z_c_pi = Z_c;

r_km = [r_km]
x_km = [x_km]
g_km = [g_km]
b_km = [b_km]

### Define more data
l = 600 #km
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
sys = System(joinpath(pwd(), file_name));
tspan = (0.0, p.sim_params.t_max)
sys = build_new_impedance_model!(sys, p)

dist = choose_disturbance(sys, p.perturbation, p)

model = MassMatrixModel
sim = PSID.Simulation(
    model, #Type of model used
    sys, #system
    pwd(), #folder to output results
    tspan, #time span
    dist, #Type of perturbation
    all_lines_dynamic = false
)

#sim = build_sim(sys, tspan, perturbation, false, p)
#show_states_initial_value(sim)

V1_r = sim.x0_init[1]
V1_i = sim.x0_init[3]
V2_r = sim.x0_init[2]
V2_i = sim.x0_init[4]
# ir = sim.x0_init[5]
# ii = sim.x0_init[6]

line_model_3 = "Multi-Segment Dynamic"

v2rs = [];

for n in 1:50
    p.N = n
    sys = System(joinpath(pwd(), file_name));
    sys = build_seg_model_Lshape!(sys, p);
    sim = build_sim(sys, tspan, dist, true, p);
    V2_r = sim.x0_init[2];
    ir = sim.x0_init[5];
    push!(v2rs, V2_r);
end

plot(v2rs)

# Check line parameters
p.N = 2
sys = System(joinpath(pwd(), file_name));
sys = build_seg_model_Lshape!(sys, p);
line = first(get_components(Line, sys))


# Need to get voltage states at t=0, and line current states

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