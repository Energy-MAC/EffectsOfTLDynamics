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
file_name = "../data/json_data/inv_v_machine.json"
# default_2_bus_line_dict - For 2 bus system
# default_9_bus_line_dict - For 9 bus system
line_dict = default_2_bus_line_dict

### Load relevant line data
impedance_csv = "../data/cable_data/dommel_data.csv"
capacitance_csv = "../data/cable_data/dommel_data_C.csv"

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
    t_max = 2.0,
)

### Extract line data from files
M = 1

# Do not change these factors for the 2 bus case
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
line_dict["BUS 1-BUS 2-i_1"] = l
line_dict["BUS 1-BUS 2-i_1_static"] = l

N = nothing
t_fault = 0.25

### Get perturbation struct
perturbation_params = get_default_perturbation(t_fault, perturbation)
# perturbation_params.crc_params = CRCParam(DynamicInverter, "generator-1-1", :V_ref, 0.95)# 
# perturbation_params.branch_trip_params = BTParam("Bus 9-Bus 6-i_1")

p_load = 1.0
q_load = 0.25
l_seg = 50 #km

load_scale = 1.0
line_scale = 1.0

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
    q_load,
    line_scale,
    load_scale
)

# Verify impedance values of raw file vs CSV data
# verifying(file_name, M, impedance_csv, capacitance_csv, p)

line_model_1 = "Algebraic"
line_model_2 = "Dynamic"
line_model_3 = "Multi-Segment Dynamic"

results_alg, sim = run_experiment(file_name, line_model_1, p);
sys = sim.sys
s = small_signal_analysis(sim)

"""
results_dyn, sim_dyn = run_experiment(file_name, line_model_2, p);
sys_dyn = sim_dyn.sys
s_dyn = small_signal_analysis(sim_dyn)

vr_alg = get_voltage_magnitude_series(results_alg, 102);
vr_dyn = get_voltage_magnitude_series(results_dyn, 102);
#ω_alg = get_state_series(results_alg, ("generator-101-1", :ω))
#ω_dyn = get_state_series(results_dyn, ("generator-101-1", :ω))


plot(vr_alg, xlabel = "time", ylabel = "vr p.u.", label = "V1")
plot!(vr_dyn, xlabel = "time", ylabel = "vr p.u.", label = "V1_dyn")
plot!(title = "Line length = "*string(p.l_dict["BUS 1-BUS 2-i_1"])*" km, perturbation = "*perturbation)

results_ms_dyn, sim_ms_dyn = run_experiment(file_name, line_model_3, p);
sys_ms_dyn = sim_ms_dyn.sys
s_ms_dyn = small_signal_analysis(sim_ms_dyn)            
vr_ms_dyn = get_voltage_magnitude_series(results_ms_dyn, 102);
plot!(vr_ms_dyn, label = "V1_ms_dyn")
#ω_ms_dyn = get_state_series(results_ms_dyn, ("generator-101-1", :ω))

M = 5
p.M = M
z_km, y_km, Z_c_abs, z_km_ω, z_km_ω_5_to_1, Z_c_5_to_1_abs = get_line_parameters(impedance_csv, capacitance_csv, M)
p.z_km = z_km;
p.y_km = y_km;
p.Z_c_abs = Z_c_abs;
p.z_km_ω = z_km_ω;
p.z_km_ω_5_to_1 = z_km_ω_5_to_1;
p.Z_c_5_to_1_abs = Z_c_5_to_1_abs;

results_ms_mb_dyn, sim_ms_mb = run_experiment(file_name, line_model_3, p);
sys_ms_mb = sim_ms_mb.sys
s_ms_mb = small_signal_analysis(sim_ms_mb)            
vr_ms_mb_dyn = get_voltage_magnitude_series(results_ms_mb_dyn, 102);
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

# line_lengths = [100, 250, 500]
# loading_scenarios = [(0.5, 0.5), (0.75, 0.25), (1.0, 0.0)]

line_scales = [1.0]#, 1.25, 1.5, 2.0]
load_scales = [1.0]#, 1.25, 1.5, 1.75]

plots = []
plt = []

using LaTeXStrings

for line_scale in line_scales
    p.line_scale = line_scale
# for l in line_lengths
#     p.l_dict["Bus 9-Bus 6-i_1"] = l
#     p.l_dict["Bus 9-Bus 6-i_1_static"] = l
    
    for load_scale in load_scales
    p.load_scale = load_scale

    # for (p_load, q_load) in loading_scenarios
    #    p.p_load = p_load
    #    p.q_load = q_load
    
        M = 1
        z_km, y_km, Z_c_abs, z_km_ω, z_km_ω_5_to_1, Z_c_5_to_1_abs = get_line_parameters(impedance_csv, capacitance_csv, M, factor_z, factor_y)
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
        plt = plot(vr_alg, label = L"$statpi$")

        results_dyn, sim_dyn = run_experiment(file_name, line_model_2, p);
        sys_dyn = sim_dyn.sys
        s_dyn = small_signal_analysis(sim_dyn)
        vr_dyn = get_voltage_magnitude_series(results_dyn, 102);
        plot!(plt, vr_dyn, label = L"$dynpi$")

        results_ms_dyn, sim_ms_dyn = run_experiment(file_name, line_model_3, p);
        sys_ms_dyn = sim_ms_dyn.sys
        s_ms_dyn = small_signal_analysis(sim_ms_dyn)            
        vr_ms_dyn = get_voltage_magnitude_series(results_ms_dyn, 102);
        plot!(plt, vr_ms_dyn, label = L"MSSB")

        M = 3
        p.M = M
        z_km, y_km, Z_c_abs, z_km_ω, z_km_ω_5_to_1, Z_c_5_to_1_abs = get_line_parameters(impedance_csv, capacitance_csv, M, factor_z, factor_y)
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
        plot!(plt, vr_ms_mb_dyn, label = L"MSMB")
        
        plot!(plt, legend = true)        
        push!(plots, plt)
    end
end

combined_plot = plot(plots..., layout=(1,1))
plot!(combined_plot, xlabel = L"$ Time \quad [s]$", title = "")
plot!(combined_plot, ylabel = L" $||V_2|| \quad [p.u.]$")

plot!(combined_plot, xlims = (0.0, 2.0))
plot!(combined_plot, dpi = 300)
Plots.savefig("../figures/Week 2/Weds/ivm_2s.png")

plot!(combined_plot, xlims = (0.249, 0.260))
plot!(combined_plot, ylims = (0.85,1.1))
Plots.savefig("../figures/Week 2/Weds/ivm_2s_zoom.png")

plot!(xlims=(0.249, 0.255))
# plot!(ylims=(0.981,0.983))
# plot!(legend = false)
# plot!(legend=:bottomright)

l = get_component(Line, sys, "BUS 1-BUS 2-i_1")
p_branch = get_activepower_branch_flow(results_alg, "BUS 1-BUS 2-i_1", :from)
p_branch_dyn = get_activepower_branch_flow(results_dyn, "BUS 1-BUS 2-i_1", :from)
p_branch_ms = get_activepower_branch_flow(results_ms_dyn, "BUS 1-BUS 2-i_1_segment_1_branch_1", :from)
p_branch_ms_mb = [get_activepower_branch_flow(results_ms_mb_dyn, "BUS 1-BUS 2-i_1_segment_1_branch_"*string(i), :from) for i in 1:5]

plot(p_branch)
plot!(p_branch_dyn)
plot!(p_branch_ms)
plot!(ylim = (0, 1))

q_branch = get_reactivepower_branch_flow(results_alg, "BUS 1-BUS 2-i_1", :from)
q_branch_dyn = get_reactivepower_branch_flow(results_dyn, "BUS 1-BUS 2-i_1", :from)
plot(q_branch)
plot!(q_branch_dyn)