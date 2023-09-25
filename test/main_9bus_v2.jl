cd(@__DIR__)
using PowerSystems
using PowerSimulationsDynamics
using Plots
using Revise
using EffectsOfTLDynamics

using CSV
using DataFrames
using Dates
using LaTeXStrings

const PSY = PowerSystems;
const PSID = PowerSimulationsDynamics;

### Choose test case
# "SIIB.json"
# "9bus.json"
# "inv_v_machine.json"
# "twobus_2inv.json"
# "9bus_slackless.json"
model = "9bus_slackless"
file_name = "../data/json_data/"*model*".json"
# default_2_bus_line_dict - For 2 bus system
# default_9_bus_line_dict - For 9 bus system
line_dict = default_9_bus_line_dict

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
perturbation_params.branch_trip_params = BTParam("Bus 7-Bus 5-i_1")



V_nom = 230 # kV
Z_o = sqrt(z_km_ω_5_to_1/y_km)
SIL = V_nom ^2/Z_o
p_load = real(SIL)/100
q_load = imag(SIL)/100

l_seg = 50 #km

load_scale = 1.0;
line_scale = 1.0;

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
);

# Verify impedance values of raw file vs CSV data
# verifying(file_name, M, impedance_csv, capacitance_csv, p)

line_model_1 = "Algebraic";
line_model_2 = "Dynamic";
line_model_3 = "Multi-Segment Dynamic";

results_alg, sim = run_experiment(file_name, line_model_1, p);
sys = sim.sys;

now_date = now()
rn = string(now_date) 

main_path = "../results/"*rn*"/"
mkdir(main_path)

partial_path = main_path*"$(p.line_scale)_$(p.load_scale)"
mkdir(partial_path)

folder_path = partial_path*"/statpi"
mkdir(folder_path)

store_bus_voltages(results_alg, sys, folder_path*"/bus_voltages.csv")
store_filter_currents(results_alg, sys, folder_path*"/filter_currents.csv")
store_branch_power_flows(results_alg, sys, folder_path*"/branch_power_flows.csv")
store_generator_speeds(results_alg, sys, folder_path*"/generator_speeds.csv")
store_branch_currents(results_alg, sys, folder_path*"/branch_currents.csv")

using PowerFlows
sol = solve_powerflow(ACPowerFlow(), sys)
sol["bus_results"]
s = small_signal_analysis(sim)

p.load_scale = 10.0
results_alg2, sim2 = run_experiment(file_name, line_model_1, p);
s2 = small_signal_analysis(sim2)
sys2 = sim2.sys;
sol2 = solve_powerflow(ACPowerFlow(), sys2)
sol2["bus_results"]

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

line_scales = collect(2.0)
load_scales = collect(0.5:0.5:3.0)

now_date = now()
rn = string(now_date)
prim_path = "../results/"*model
mkdir(prim_path)
main_path = prim_path*"/"*rn*"/"
mkdir(main_path)

for line_scale in line_scales
    p.line_scale = line_scale

    for load_scale in load_scales
        p.load_scale = load_scale
        
        partial_path = main_path*"$(p.line_scale)_$(p.load_scale)"
        mkdir(partial_path)

        M = 1
        p.M = M
        z_km, y_km, Z_c_abs, z_km_ω, z_km_ω_5_to_1, Z_c_5_to_1_abs = get_line_parameters(impedance_csv, capacitance_csv, M, factor_z, factor_y)
        p.z_km = z_km
    
        folder_path = partial_path*"/statpi"
        mkdir(folder_path)
        results_alg, sim = run_experiment(file_name, line_model_1, p);
        sys = sim.sys
        store_bus_voltages(results_alg, sys, folder_path*"/bus_voltages.csv")
        store_filter_currents(results_alg, sys, folder_path*"/inverter_currents.csv")
        #store_branch_power_flows(results_alg, sys, folder_path*"/branch_power_flows.csv")
        store_generator_speeds(results_alg, sys, folder_path*"/generator_speeds.csv")
        #store_branch_currents(results_alg, sys, folder_path*"/branch_currents.csv")


        folder_path = partial_path*"/dynpi"
        mkdir(folder_path)
        results_dyn, sim_dyn = run_experiment(file_name, line_model_2, p);
        sys_dyn = sim_dyn.sys
        store_bus_voltages(results_dyn, sys_dyn, folder_path*"/bus_voltages.csv")
        store_filter_currents(results_dyn, sys_dyn, folder_path*"/inverter_currents.csv")
        #store_branch_power_flows(results_dyn, sys_dyn, folder_path*"/branch_power_flows.csv")
        store_generator_speeds(results_dyn, sys_dyn, folder_path*"/generator_speeds.csv")
        #store_branch_currents(results_dyn, sys_dyn, folder_path*"/branch_currents.csv")

        folder_path = partial_path*"/MSSB"
        mkdir(folder_path)
        results_ms_dyn, sim_ms_dyn = run_experiment(file_name, line_model_3, p);
        sys_ms_dyn = sim_ms_dyn.sys
        store_bus_voltages(results_ms_dyn, sys_ms_dyn, folder_path*"/bus_voltages.csv")
        store_filter_currents(results_ms_dyn, sys_ms_dyn, folder_path*"/inverter_currents.csv")
        #store_branch_power_flows(results_ms_dyn, sys_ms_dyn, folder_path*"/branch_power_flows.csv")
        store_generator_speeds(results_ms_dyn, sys_ms_dyn, folder_path*"/generator_speeds.csv")
        #store_branch_currents(results_ms_dyn, sys_ms_dyn, folder_path*"/branch_currents.csv")

        M = 3
        p.M = M
        z_km, y_km, Z_c_abs, z_km_ω, z_km_ω_5_to_1, Z_c_5_to_1_abs = get_line_parameters(impedance_csv, capacitance_csv, M, factor_z, factor_y)
        p.z_km = z_km;
        p.y_km = y_km;
        p.Z_c_abs = Z_c_abs;
        p.z_km_ω = z_km_ω;
        p.z_km_ω_5_to_1 = z_km_ω_5_to_1;
        p.Z_c_5_to_1_abs = Z_c_5_to_1_abs;

        folder_path = partial_path*"/MSMB"
        mkdir(folder_path)  
        results_ms_mb_dyn, sim_ms_mb = run_experiment(file_name, line_model_3, p)
        sys_ms_mb = sim_ms_mb.sys
        store_bus_voltages(results_ms_mb_dyn, sys_ms_mb, folder_path*"/bus_voltages.csv")
        store_filter_currents(results_ms_mb_dyn, sys_ms_mb, folder_path*"/inverter_currents.csv")
        #store_branch_power_flows(results_ms_mb_dyn, sys_ms_mb, folder_path*"/branch_power_flows.csv")
        store_generator_speeds(results_ms_mb_dyn, sys_ms_mb, folder_path*"/generator_speeds.csv")
        #store_branch_currents(results_ms_mb_dyn, sys_ms_mb, folder_path*"/branch_currents.csv")

    end
end

# Not necessary?
# results_alg, sim, sys, s, vr_alg = nothing, nothing, nothing, nothing, nothing;
# results_dyn, sim_dyn, sys_dyn, s_dyn, vr_dyn = nothing, nothing, nothing, nothing, nothing;
# results_ms_dyn, seg_sim, seg_sys, s_seg, vr_ms_dyn = nothing, nothing, nothing, nothing, nothing;
# results_ms_mb_dyn, sim_ms_mb, sys_ms_mb, s_ms_mb, vr_ms_mb_dyn = nothing, nothing, nothing, nothing, nothing;

plots = []
plt = []

# Plotting loop
plt = plot(vr_alg, label = L"\mathrm{algebraic}")
plot!(plt, vr_dyn, label = L"\mathrm{dynpi}")
plot!(plt, vr_ms_dyn, label = L"\mathrm{MSSB}")
plot!(plt, vr_ms_mb_dyn, label = L"\mathrm{MSMB}")
plot!(plt, legend = true)        
push!(plots, plt)

combined_plot = plot(plots..., layout=(1,1))
plot!(combined_plot, xlabel = L"$ \mathrm{Time} \quad [s]$", title = "")
plot!(combined_plot, ylabel = L"$||V_2|| \quad \mathrm{[\ p.u.]}$")

plot!(combined_plot, xlims = (0.0, 2.0))
plot!(combined_plot, dpi = 300)
Plots.savefig("../figures/Week 2/Thurs/ivm_2s.png")
plot!(combined_plot, title = L"$\mathrm{Load \ scale} = %$load_scale \ \mathrm{Line \ scale} = %$line_scale$")
annotate!(0.5,1, (L"Label", 8, :red, :top))
plot!(combined_plot, legend_title = L"$\mathrm{Load \ scale} = %$load_scale \ \mathrm{Line \ scale} = %$line_scale$")
plot!(combined_plot, xlims = (0.249, 0.260))
plot!(combined_plot, ylims = (0.85,1.1))
Plots.savefig("../figures/Week 2/Thurs/ivm_2s_zoom.png")

plot!(xlims=(0.249, 0.255))
# plot!(ylims=(0.981,0.983))
# plot!(legend = false)
# plot!(legend=:bottomright)

l = get_component(Line, sys, "BUS 1-BUS 2-i_1")
p_branch = get_activepower_branch_flow(results_alg, "BUS 1-BUS 2-i_1", :from)
p_branch_dyn = get_activepower_branch_flow(results_alg2, "BUS 1-BUS 2-i_1", :from)
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