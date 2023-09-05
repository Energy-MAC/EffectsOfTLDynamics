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
# default_2_bus_line_dict - For 2 bus system
# default_9_bus_line_dict - For 9 bus system
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
perturbation = "BranchTrip"

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
Z_c = 380 # Ω
r_km = 0.05 # Ω/km
x_km = 0.488 # Ω/km
g_km = 0.0 # S/km
b_km = 3.371e-6 # S/km
z_km = r_km + im*x_km
y_km = im*b_km
z_km_ω = z_km[1]

### Define more data
l = 100 #, 500, 750, 100 #km
# line_dict["BUS 1-BUS 2-i_1"] = l
# line_dict["BUS 1-BUS 2-i_1_static"] = l

N = nothing
t_fault = 0.25

### Get perturbation struct
perturbation_params = get_default_perturbation(t_fault, perturbation)
#perturbation_params.crc_params = CRCParam(DynamicInverter, "generator-1-1", :V_ref, 0.95)

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
    perturbation_params,
    nothing,
    nothing
)

# Verify impedance values of raw file vs CSV data
verifying(file_name, M, impedance_csv, capacitance_csv, p)

line_model_1 = "Algebraic"
results_alg, sim = run_experiment(file_name, line_model_1, p);
sys = sim.sys
s = small_signal_analysis(sim)

line_model_2 = "Dynamic"
results_dyn, sim_dyn = run_experiment(file_name, line_model_2, p);
sys_dyn = sim_dyn.sys
s_dyn = small_signal_analysis(sim_dyn)

vr_alg = get_state_series(results_alg, ("generator-102-1", :vr_filter));
vr_dyn = get_state_series(results_dyn, ("generator-102-1", :vr_filter));

plot(vr_alg, xlabel = "time", ylabel = "vr p.u.", label = "vr")
plot!(vr_dyn, xlabel = "time", ylabel = "vr p.u.", label = "vr_dyn")
plot!(title = "Line length = "*string(p.l)*" km, perturbation = "*perturbation)

line_model_3 = "Multi-Segment Dynamic"
results_ms_dyn, seg_sim, seg_sys, s_seg = nothing, nothing, nothing, nothing
for n in [5]
    for l in [100, 300, 500]
        # initialize storage for simulation results
        vr_ms_dyns = Vector{Vector{Float64}}(undef, 3)
        count = 1
        for (p_load, q_load) in [(0.5, 0.5), (0.75, 0.25), (1.0, 0.0)]
            p.N = n
            p.l = l
            p.p_load = p_load
            p.q_load = q_load
            results_ms_dyn, seg_sim, seg_sys, s_seg = nothing, nothing, nothing, nothing
            results_ms_dyn, seg_sim = run_experiment(file_name, line_model_3, p)
            seg_sys = seg_sim.sys
            s_seg = small_signal_analysis(seg_sim)

            vr_ms_dyn = get_state_series(results_ms_dyn, ("generator-102-1", :vr_filter));
            vr_ms_dyns[count] = vr_ms_dyn
            count += 1
        end
        # Plot the results for all 3 load values
        plot(vr_ms_dyns[1], xlabel = "time", ylabel = "vr p.u.", label = "vr_segs_$(p.N)_branch_$(p.M)")
        plot!(vr_ms_dyns[2], xlabel = "time", ylabel = "vr p.u.", label = "vr_segs_$(p.N)_branch_$(p.M)")
        plot!(vr_ms_dyns[3], xlabel = "time", ylabel = "vr p.u.", label = "vr_segs_$(p.N)_branch_$(p.M)")
        plot!(title = "Line length = "*string(p.l)*" km, perturbation = "*perturbation)
        # Display plot
        display(plot!())
    end
    print(n)
end

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