using PowerSystems
using PowerSimulationsDynamics
using PowerSystemCaseBuilder

threebus_sys = build_system(PSIDSystems, "Three Bus Dynamic data Example System")
threebus_sys_dyn = deepcopy(threebus_sys);

dyn_branch = DynamicBranch(get_component(Line, threebus_sys_dyn,"BUS 2-BUS 3-i_1"))
add_component!(threebus_sys_dyn, dyn_branch)

g = get_component(DynamicInverter, threebus_sys_dyn, "generator-103-1")
dist = ControlReferenceChange(3.0, g, :V_ref, 0.95)

sim_dyn = Simulation(
           ResidualModel, #Type of model used
           threebus_sys_dyn, #system
           pwd(), #folder to output results
           (0.0, 5.0), #time span
           dist, #Type of perturbation
       )

show_states_initial_value(sim_dyn)

execute!(
           sim_dyn, #simulation structure
           IDA(), #Sundials DAE Solver
           dtmax = 0.02, #Maximum step size
       )

results_dyn = read_results(sim_dyn)
ilr_dyn = get_state_series(results_dyn, ("BUS 2-BUS 3-i_1", :Il_R))[2][1] # 2 is values, 1 is first timestep
ili_dyn = get_state_series(results_dyn, ("BUS 2-BUS 3-i_1", :Il_I))[2][1]
v1r_dyn = get_state_series(results_dyn, ("V_3", :R))[2][1]
v1i_dyn = get_state_series(results_dyn, ("V_3", :I))[2][1]
v2r_dyn = get_state_series(results_dyn, ("V_2", :R))[2][1]
v2i_dyn = get_state_series(results_dyn, ("V_2", :I))[2][1]

# Solve for r, x
# v1r - v2r - Rilr + omega*L*ili = 0 
# v1i - v2i - Rili - omega*L*ilr = 0

# R, L
A = [ilr_dyn -ili_dyn; ili_dyn ilr_dyn];
b = [v1r_dyn-v2r_dyn; v1i_dyn-v2i_dyn];

# Ax = b --> x = (A^-1)b
x = inv(A)*-b
r_line = x[1]
x_line = x[2]


dyn_line = first(get_components(DynamicBranch, threebus_sys_dyn))

dyn_line.branch.r
dyn_line.branch.x
