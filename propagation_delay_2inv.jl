# Plot voltage profiles at either end of the line to see if an actual propagation delay can be observed 

include("all_sims.jl")
include("extra_functions.jl")
using InfrastructureSystems
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
l = 200; # line length km
N = 20;

p = ExpParams(N, l, Z_c, r_km, x_km, g_km, b_km,abstol, reltol, maxiters)
sys_copy = deepcopy(sys) # need to do this bc building seg model modifies sys 
sys_ms = build_seg_model(sys_copy, p)
sim = build_sim(sys_ms, tspan, perturbation, true);
show_states_initial_value(sim)

execute_sim(sim, p);
results = results_sim(sim);

#v1 = get_state_series(results, ("generator-102-1", :vr_filter));
v1d = get_state_series(results, ("V_1", :R));
v1q = get_state_series(results, ("V_1", :I));

v2d = get_state_series(results, ("V_2", :R));
v2q = get_state_series(results, ("V_2", :I));

plot(grid=true)
#plot!(v1, xlabel = "time", ylabel = "vr p.u.", label = "vr")
plot!(v1d, xlabel="time", ylabel="p.u.", label="v1d, N="*string(N))
plot!(v2d, xlabel="time", label="v2d, N="*string(N))

xlims!(0,1)
ylims!(0.95, 1.05)
title!("N="*string(N)*", L="*string(l)*"km")
title!("L="*string(l)*"km")
ylabel!("voltage p.u.")
xlabel!("time")

save_current_fig()