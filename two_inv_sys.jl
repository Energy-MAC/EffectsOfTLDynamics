# Want to build an inverter vs inverter system and find an equilibrium point using powerflow on the nonlinear system. Then take this equilibrium point and do a linearization around it separately. 
include("all_sims.jl")
using PowerSystems
using InfrastructureSystems
using PowerSimulationsDynamics
using PlotlyJS


file_name = "test_sys.json";
sys = System(joinpath(pwd(), file_name));

inv = get_component(DynamicInverter, sys, "generator-102-1"); # Get first inverter
inv2 = deepcopy(inv); # Create second inverter 
inv2.name = "generator-101-1"; # Rename second inverter

inv_static = get_component(StaticInjection, sys, "generator-102-1"); # Get static injection associated with first inverter
inv2_static = deepcopy(inv_static); # Copy for second inverter
inv2_static.name = "generator-101-1"; # Rename static injection to match dynamic component 
inv2_static.dynamic_injector = inv2; # assign second inverter dynamic injector to second inverter static injection object 
inv2_static.time_series_container = InfrastructureSystems.TimeSeriesContainer(); # To fix an error I was getting 


inf_bus = get_component(Source, sys, "InfBus"); # Get the infinite bus 

inv2_static.bus = inf_bus.bus; # Use the bus info from inf bus to update the bus of inv2 
# inv2_static.bus.bustype = BusTypes.SLACK; # May need to change bus type?

add_component!(sys, inv2_static) # add the second inverter 
remove_component!(sys,inf_bus) # remove the infinite bus source 

# # Build sim - algebraic lines
# sim = Simulation(
#     MassMatrixModel,
#         sys,
#         pwd(),
#         (0.0, 20.0);
#     )

# show_states_initial_value(sim);

# # Build sim - dynamic lines
# sim_dyn = Simulation(
#     MassMatrixModel,
#         sys,
#         pwd(),
#         (0.0, 20.0);
#         all_lines_dynamic = true
#     )

# show_states_initial_value(sim_dyn);

# Do time series simulation
tmax = 10.0;
tspan = (0.0, tmax);
dist = "CRC";
perturbation = choose_disturbance(sys, dist);

abstol = 1e-13;
reltol = 1e-10;
maxiters = Int(1e10);
p = ExpParams(0,0,0,0,0,0,0,abstol, reltol, maxiters);

# ALGEBRAIC 
sim_alg = build_sim(sys, tspan, perturbation, false);
execute_sim(sim_alg, p);
results_alg = results_sim(sim_alg);
stb_alg = small_signal_analysis(sim_alg);
summary_eigenvalues(stb_alg)


# DYN LINES 
sim_dyn = build_sim(sys, tspan, perturbation, true);
execute_sim(sim_dyn, p);
results_dyn = results_sim(sim_dyn);
stb_dyn = small_signal_analysis(sim_dyn);
summary_eigenvalues(stb_dyn)

vr_alg = get_state_series(results_alg, ("generator-102-1", :vr_filter));
plot(vr_alg, xlabel = "time", ylabel = "vr p.u.", label = "vr")

# PLOT RESULTS 

vr_dyn = get_state_series(results_dyn, ("generator-102-1", :vr_filter));
plot!(vr_dyn, xlabel = "time", ylabel = "vr p.u.", label = "vr_dyn")
xlims!(0.995,1.01)


plotlyjs()
plot()
plot!(real(stb_alg.eigenvalues), imag(stb_alg.eigenvalues), seriestype=:scatter, label="algebraic")
plot!(real(stb_dyn.eigenvalues), imag(stb_dyn.eigenvalues), seriestype=:scatter, label="dynamic")
xlabel!("Real")
ylabel!("Imag")
