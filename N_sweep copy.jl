include("./all_sims.jl") # runs all_sims file but loads all functions defined in that file.
include("extra_functions.jl")

file_name = "test_sys.json"
sys = make_sys_from_filename(file_name);
t_max = 30.0
dist = "CRC" # control reference change (on inverter)

Z_c = 380 # Ω
r_km = 0.05 # Ω/km
x_km = 0.488 # Ω/km
g_km = 0 # S/km
b_km = 3.371e-6 # S/km
abstol = 1e-13
reltol = 1e-10
maxiters = Int(1e10)

Lrange = [50, 300, 500, 1000] # km 
Nrange = 1:5; # Define segment range 

plot_N_vs_xth_largest_λ_for_Ls(3, sys, Lrange, Nrange, Z_c, r_km, x_km, g_km, b_km,abstol, reltol, maxiters, perturbation, tspan)

plot_all_eigs_for_Ns(sys, Nrange, Z_c, r_km, x_km, g_km, b_km, abstol, reltol, maxiters, perturbation, tspan)