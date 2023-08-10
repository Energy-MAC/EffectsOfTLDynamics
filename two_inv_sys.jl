
include("all_sims.jl")
include("extra_functions.jl")
using PowerSystems
using InfrastructureSystems
using PowerSimulationsDynamics

PSID = PowerSimulationsDynamics;

sys = get_modified_twobus_sys();

t_max = 5.0;
tspan = (0.0, t_max);
dist = "CRC";
perturbation = choose_disturbance(sys, dist);

# Define line parameters 
Z_c = 380; # Ω
r_km = 0.05; # Ω/km
x_km = 0.488; # Ω/km
g_km = 0; # S/km
b_km = 3.371e-6; # S/km
abstol = 1e-13;
reltol = 1e-10;
maxiters = Int(1e10);

l = 100; # line length km
N = 2; # number of segments 

## - try to extract initial state values
p = ExpParams(N, l, Z_c, r_km, x_km, g_km, b_km,abstol, reltol, maxiters)
og_sys_copy = deepcopy(sys) # need to do this bc building seg model modifies sys 
sys_ms = build_seg_model(og_sys_copy, p)
sim = build_sim(sys_ms, tspan, perturbation, true);

show_states_initial_value(sim)

### ---- Do small signal analysis on MS lines 

plot_N_vs_xth_largest_λ_for_Ls(2, sys, Lrange, Nrange, Z_c, r_km, x_km, g_km, b_km,abstol, reltol, maxiters, perturbation, tspan)

plot_all_eigs_for_Ns(sys, Nrange, Z_c, r_km, x_km, g_km, b_km,abstol, reltol, maxiters, perturbation, tspan)

Lseg = 50;
Lrange = 100:50:500;
plot_nonzero_eig_fixed_Lseg_for_Ls(2, Lseg, sys, Lrange, Z_c, r_km, x_km, g_km, b_km,abstol, reltol, maxiters, perturbation, tspan)



function print_initial_states(sim)
    # get a list of states 
    inputs = sim.inputs
    global_state_map = PSID.make_global_state_map(inputs)

    for device in PSID.get_dynamic_injectors(inputs)
        states = PSY.get_states(device)
        name = PSY.get_name(device)
        println(name)
        println("====================")
        global_index = global_state_map[name]
        for s in states
            print(s, " ", round(sim.x0_init[global_index[s]]; digits = 4), "\n")
        end
        println("====================")
    end
    dyn_branches = PSID.get_dynamic_branches(inputs)
    if !isempty(dyn_branches)
        for br in dyn_branches
            states = PSY.get_states(br)
            name = PSY.get_name(br)
            global_index = global_state_map[name]
            x0_br = Dict{Symbol, Float64}()
            for (i, s) in enumerate(states)
                print(s, " ", round(sim.x0_init[global_index[s]]; digits = 5), "\n")
            end
        end
    end
    return

    # # get their value at t=0
    # init_val = x[2][1]; 

end





# # ALGEBRAIC 
# sim_alg = build_sim(sys, tspan, perturbation, false);
# show_states_initial_value(sim_alg)
# # execute_sim(sim_alg, p);
# # results_alg = results_sim(sim_alg);
# stb_alg = small_signal_analysis(sim_alg);
# summary_eigenvalues(stb_alg)

# # DYN LINES 
# sim_dyn = build_sim(sys, tspan, perturbation, true);
# show_states_initial_value(sim_dyn)
# # execute_sim(sim_dyn, p);
# # results_dyn = results_sim(sim_dyn);
# stb_dyn = small_signal_analysis(sim_dyn);
# summary_eigenvalues(stb_dyn)



# # PLOT RESULTS 
# vr_alg = get_state_series(results_alg, ("generator-102-1", :vr_filter));
# plot()
# plot!(vr_alg, xlabel = "time", ylabel = "vr p.u.", label = "vr")

# vr_dyn = get_state_series(results_dyn, ("generator-102-1", :vr_filter));
# plot!(vr_dyn, xlabel = "time", ylabel = "vr p.u.", label = "vr_dyn")

# vr_ms = get_state_series(results_ms, ("generator-102-1", :vr_filter))
# display(plot!(vr_ms, xlabel = "time", ylabel = "vr p.u.", label = "vr_$(N)"))

# end

# xlims!(0.995,1.3)


# # Plot eigenvalues 
# #plotlyjs()
# plot()
# plot!(real(stb_alg.eigenvalues), imag(stb_alg.eigenvalues), seriestype=:scatter, label="algebraic")
# plot!(real(stb_dyn.eigenvalues), imag(stb_dyn.eigenvalues), seriestype=:scatter, label="dynamic")
# xlabel!("Real")
# ylabel!("Imag")

