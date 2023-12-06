cd(@__DIR__)
using PowerSystems
using PowerSimulationsDynamics
using InfrastructureSystems
using Sundials
using Plots
using EffectsOfTLDynamics

const PSY = PowerSystems;
const PSID = PowerSimulationsDynamics;
const ETL = EffectsOfTLDynamics

# Load system 
file_name = "GFM_v_load.json";

impedance_csv = "../data/cable_data/dommel_data.csv"
capacitance_csv = "../data/cable_data/dommel_data_C.csv"
line_dict = default_2_bus_line_dict

# get load 
l_device = get_component(ElectricLoad, sys, "load1021")
transform_load_to_constant_impedance(l_device) # convert to constant impedance 
l_device.impedance_active_power = 1.0
l_device.impedance_reactive_power = 0.3

perturbation = LoadChange(0.025, l_device, :P_ref_impedance, 0.8) # define load step change 



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

z_km, y_km, Z_c_abs, z_km_ω, z_km_ω_5_to_1, Z_c_5_to_1_abs = get_line_parameters(impedance_csv, capacitance_csv, M, 1.0, 1.0)

### Define more data
l = 100 #km
line_dict["BUS 1-BUS 2-i_1"] = l
line_dict["BUS 1-BUS 2-i_1_static"] = l

perturbation_type = "LoadChange";

perturbation_params = get_default_perturbation(0.25, perturbation_type)
perturbation_params.load_change_params.load_name = "load1021";
perturbation_params.load_change_params.var_to_change = :P_ref_impedance;
perturbation_params.load_change_params.ref_value = 0.8;

p = ExpParams(
    nothing, 
    1, 
    100,
    20, 
    Z_c_abs, 
    z_km,
    y_km,
    z_km_ω,
    z_km_ω_5_to_1,
    Z_c_5_to_1_abs,
    line_dict,
    sim_p, 
    perturbation_type, 
    perturbation_params,
    0.0,
    0.0,
    1.0,
    1.0
);

sys = System("GFM_v_load.json")

sys = ETL.build_seg_model!(sys, p, false, "");

sys = build_new_impedance_model!(sys, p, false, "");


    sim = build_sim(sys, tspan, perturbation, dyn_lines, p);

algebraic  = "Algebraic"
dynamic = "Dynamic"
mssb = "MSSB"

results_alg, sim = run_experiment(file_name, algebraic, p);


# Build sim 
sim_dyn = PSID.Simulation(
           ResidualModel, # Type of formulation
           sys, # System
           mktempdir(), # Output directory
           time_span,
           l_change,
           all_lines_dynamic=true)

# Small sig analysis 
ss_dyn = small_signal_analysis(sim_dyn)
evals_dyn = ss_dyn.eigenvalues

execute!(
    sim_alg, #simulation structure
    IDA(), #Sundials DAE Solver
    dtmax = 0.02, #Arguments: Maximum timestep allowed
    );

results_alg = read_results(sim_alg)

execute!(
    sim_dyn, #simulation structure
    IDA(), #Sundials DAE Solver
    dtmax = 0.02, #Arguments: Maximum timestep allowed
    );

results_dyn = read_results(sim_dyn)


#get_state_series(results, ("generator-101-1", :));

# Plot eigs 
plot(real(evals), imag(evals), seriestype=:scatter, label="statpi")

plot!(real(evals), imag(evals), seriestype=:scatter, label="dynpi")

xlims!(-100,1)


volt = get_voltage_magnitude_series(results, 101);
plot(volt, xlabel = "time", ylabel = "Voltage [pu]", label = "V_1")

volt = get_voltage_magnitude_series(results, 102);
plot!(volt, xlabel = "time", ylabel = "Voltage [pu]", label = "V_2")

current = get_real_current_branch_flow(results, "BUS 1-BUS 2-i_1")
plot(current, xlabel="time", ylabel="Real I")

current = get_imaginary_current_branch_flow(results, "BUS 1-BUS 2-i_1")
plot!(current, xlabel="time", ylabel="Imag I")
