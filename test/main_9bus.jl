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
file_name = "../data/json_data/9bus_slackless.json"
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

factor_z = 1.0
factor_y = 1.0
z_km, y_km, Z_c_abs, z_km_ω, z_km_ω_5_to_1, Z_c_5_to_1_abs = get_line_parameters(impedance_csv, capacitance_csv, M, factor_z, factor_y)

# Kundur parameters for testing
# Z_c = 380 # Ω
# r_km = 0.05 # Ω/km
# x_km = 0.488 # Ω/km
# g_km = 0.0 # S/km
# b_km = 3.371e-6 # S/km

### Define more data
l = 100 #km
# line_dict["BUS 1-BUS 2-i_1"] = l
# line_dict["BUS 1-BUS 2-i_1_static"] = l

N = nothing
t_fault = 0.25

### Get perturbation struct
# Longest line
perturbation_params = get_default_perturbation(t_fault, perturbation)
# perturbation_params.branch_trip_params = BTParam("Bus 9-Bus 6-i_1")
# Heaviest loaded line
perturbation_params.branch_trip_params = BTParam("Bus 7-Bus 5-i_1")

p_load = 0.5
q_load = 0.5
l_seg = 50 #km

p = ExpParams(
    N, 
    M, 
    l,
    l_seg, 
    Z_c_abs, 
    z_km,
    y_km,
    z_km_ω,
    z_km_ω_5_to_1,
    Z_c_5_to_1_abs,
    line_dict,
    sim_p, 
    perturbation, 
    perturbation_params,
    p_load,
    q_load
)

# Verify impedance values of raw file vs CSV data
# verifying(file_name, M, impedance_csv, capacitance_csv, p)

line_model_1 = "Algebraic"
line_model_2 = "Dynamic"
line_model_3 = "Multi-Segment Dynamic"

results_alg, sim = run_experiment(file_name, line_model_1, p);
sys = sim.sys
s = small_signal_analysis(sim)

# Find heaviest loaded line
# sol = solve_powerflow(ACPowerFlow(), sys)
# sol["flow_results"]

results_dyn, sim_dyn = run_experiment(file_name, line_model_2, p);
sys_dyn = sim_dyn.sys
s_dyn = small_signal_analysis(sim_dyn)

results_ms_dyn, sim_ms_dyn = run_experiment(file_name, line_model_3, p);
sys_ms_dyn = sim_ms_dyn.sys
s_ms_dyn = small_signal_analysis(sim_ms_dyn)   

M = 5
p.M = M
z_km, y_km, Z_c_abs, z_km_ω, z_km_ω_5_to_1, Z_c_5_to_1_abs = get_line_parameters(impedance_csv, capacitance_csv, M, factor_z, factor_y)
p.z_km = z_km;
p.y_km = y_km;
p.Z_c_abs = Z_c_abs;
p.z_km_ω = z_km_ω;
p.z_km_ω_5_to_1 = z_km_ω_5_to_1;
p.Z_c_5_to_1_abs = Z_c_5_to_1_abs;

results_ms_mb_dyn, sim_ms_mb = run_experiment(file_name, line_model_3, p);
sys_ms_mb = sim_ms_mb.sys
s_ms_mb = small_signal_analysis(sim_ms_mb)   

vr_alg = get_voltage_magnitude_series(results_alg, 6);
vr_dyn = get_voltage_magnitude_series(results_dyn, 6);
#ω_alg = get_state_series(results_alg, ("generator-101-1", :ω))
#ω_dyn = get_state_series(results_dyn, ("generator-101-1", :ω))

plot(vr_alg, xlabel = "time", ylabel = "vr p.u.", label = "V1")
plot!(vr_dyn, xlabel = "time", ylabel = "vr p.u.", label = "V1_dyn")
# plot!(title = "Line length = "*string(p.l_dict["BUS 1-BUS 2-i_1"])*" km, perturbation = "*perturbation)
    
vr_ms_dyn = get_voltage_magnitude_series(results_ms_dyn, 6);
plot!(vr_ms_dyn, label = "V1_ms_dyn")
#ω_ms_dyn = get_state_series(results_ms_dyn, ("generator-101-1", :ω))

vr_ms_mb_dyn = get_voltage_magnitude_series(results_ms_mb_dyn, 6);
#ω_ms_mb_dyn = get_state_series(results_ms_mb_dyn, ("generator-101-1", :ω))

plot!(vr_ms_mb_dyn, label = "V1_ms_mb_dyn")
plot!(title = "Line length = "*string(p.l)*" km, perturbation = "*perturbation)

#savefig("../figures/diff M.png")

plot(ω_alg)
plot!(ω_dyn)
plot!(ω_ms_dyn)
plot!(ω_ms_mb_dyn)

results_alg, sim, sys, s, vr_alg = nothing, nothing, nothing, nothing, nothing;
results_dyn, sim_dyn, sys_dyn, s_dyn, vr_dyn = nothing, nothing, nothing, nothing, nothing;
results_ms_dyn, seg_sim, seg_sys, s_seg, vr_ms_dyn = nothing, nothing, nothing, nothing, nothing;
results_ms_b_dyn, sim_ms_mb, sys_ms_mb, s_ms_mb, vr_ms_mb_dyn = nothing, nothing, nothing, nothing, nothing;

"""

line_lengths = [100, 250, 500]
loading_scenarios = [(0.5, 0.5), (0.75, 0.25), (1.0, 0.0)]

plots = []
plt = []

for l in line_lengths
    p.l_dict["BUS 1-BUS 2-i_1"] = l
    p.l_dict["BUS 1-BUS 2-i_1_static"] = l
    for (p_load, q_load) in loading_scenarios
        p.p_load = p_load
        p.q_load = q_load
    
        M = 1
        z_km, y_km, Z_c_abs, z_km_ω, z_km_ω_5_to_1, Z_c_5_to_1_abs = get_line_parameters(impedance_csv, capacitance_csv, M)
        p.M = M
        p.z_km = z_km
    
        results_alg, sim, sys, s, vr_alg = nothing, nothing, nothing, nothing, nothing;
        results_dyn, sim_dyn, sys_dyn, s_dyn, vr_dyn = nothing, nothing, nothing, nothing, nothing;
        results_ms_dyn, seg_sim, seg_sys, s_seg, vr_ms_dyn = nothing, nothing, nothing, nothing, nothing;
        results_ms_mb_dyn, sim_ms_mb, sys_ms_mb, s_ms_mb, vr_ms_mb_dyn = nothing, nothing, nothing, nothing, nothing;

        results_alg, sim = run_experiment(file_name, line_model_1, p);
        sys = sim.sys
        s = small_signal_analysis(sim)
        vr_alg = get_voltage_magnitude_series(results_alg, 102);
        plt = plot(vr_alg, label = "V1_alg")

        results_dyn, sim_dyn = run_experiment(file_name, line_model_2, p);
        sys_dyn = sim_dyn.sys
        s_dyn = small_signal_analysis(sim_dyn)
        vr_dyn = get_voltage_magnitude_series(results_dyn, 102);
        plot!(plt, vr_dyn, label = "V1_dyn")

        results_ms_dyn, sim_ms_dyn = run_experiment(file_name, line_model_3, p);
        sys_ms_dyn = sim_ms_dyn.sys
        s_ms_dyn = small_signal_analysis(sim_ms_dyn)            
        vr_ms_dyn = get_voltage_magnitude_series(results_ms_dyn, 102);
        plot!(plt, vr_ms_dyn, label = "V1_ms_dyn")

        M = 5
        p.M = M
        z_km, y_km, Z_c_abs, z_km_ω, z_km_ω_5_to_1, Z_c_5_to_1_abs = get_line_parameters(impedance_csv, capacitance_csv, M)
        p.z_km = z_km;
        p.y_km = y_km;
        p.Z_c_abs = Z_c_abs;
        p.z_km_ω = z_km_ω;
        p.z_km_ω_5_to_1 = z_km_ω_5_to_1;
        p.Z_c_5_to_1_abs = Z_c_5_to_1_abs;

        results_ms_mb_dyn, sim_ms_mb = run_experiment(file_name, line_model_3, p)
        sys_ms_mb = sim_ms_mb.sys
        s_ms_mb = small_signal_analysis(sim_ms_mb)            
        vr_ms_mb_dyn = get_voltage_magnitude_series(results_ms_mb_dyn, 102);
        plot!(plt, vr_ms_mb_dyn, label = "V1_ms_mb_dyn")
        
        plot!(plt, legend = true, title = "l = $l, p = $p_load, q = $q_load")        
        push!(plots, plt)
    end
end

combined_plot = plot(plots..., layout=(3,3))
plot!(combined_plot, legend = true, layout = (1,3))
# savefig("../figures/2s_zoom.png")

plot!(xlims=(0.249, 0.255))
# plot!(ylims=(0.981,0.983))
# plot!(legend = false)
# plot!(legend=:bottomright)