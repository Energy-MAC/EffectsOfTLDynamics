include("./all_sims.jl") # runs all_sims file but loads all functions defined in that file.

# Sweep number of segments used in line model. 

file_name = "test_sys.json"
t_max = 30.0
dist = "CRC" # control reference change (on inverter)

Z_c = 380 # Ω
r_km = 0.05 # Ω/km
x_km = 0.488 # Ω/km
g_km = 0 # S/km
b_km = 3.371e-6 # S/km
l = 1000 #km

Nrange = [1,5,10,15,20, 25, 30, 35, 40, 45, 50]; # Define segment range 
stb = []
max_λ = [] # for storing max real eigenvalue

plot()
for N in Nrange;
    p = ExpParams(N, l, Z_c, r_km, x_km, g_km, b_km)
    sim = get_ms_dyn_sim(file_name, t_max, dist, p)
    #inv = get_component(DynamicInverter, sim.sys, "generator-102-1")
    #inv.inner_control.kpv = inv.inner_control.kpv*2 # increase gain 
    #inv.inner_control.kpc = inv.inner_control.kpc*2 # increase gain
    ss = small_signal_analysis(sim)
    push!(stb, ss)
    push!(max_λ, maximum(real(ss.eigenvalues)))
    display(plot!(real(ss.eigenvalues), imag(ss.eigenvalues), seriestype=:scatter, label="N="*string(N)))
end
xlims!(-1e-8,1e-8)

# Plot maximum real component of system eigenvalues as a function of N
plot(Nrange, max_λ, seriestype=:scatter, title="Line length= "*string(l)*"km", xlabel='N', ylabel="max real λ")


# Example of accessing participation factors 
# Each bus has a voltage state (Vr and Vi), same for each line segment (Ir, Ii), and the inverter has many states.
stb[1].participation_factors
# Inverter participation factors
pf1 = stb[1].participation_factors["generator-102-1"]

plot()
for (key, val) in pf1;
    print(key)
    display(plot!(val))
end

