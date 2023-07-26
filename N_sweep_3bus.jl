include("./all_sims_3bus.jl") # runs all_sims file but loads all functions defined in that file.
include("extra_functions.jl")

# Sweep number of segments used in line model. 
file_name = "threebus_sys.json"
t_max = 2.0
dist = "CRC" # control reference change (on inverter)

Z_c = 380 # Ω
r_km = 0.05 # Ω/km
x_km = 0.488 # Ω/km
g_km = 0 # S/km
b_km = 3.371e-6 # S/km
abstol = 1e-13
reltol = 1e-10
maxiters = Int(1e10)

lrange = [800]
Nrange = 1:10:100; # Define segment range 
stb = [] # small sig stability results

# plot()
# for N in Nrange;
#     p = ExpParams(N, l, Z_c, r_km, x_km, g_km, b_km,abstol, reltol, maxiters)
#     sim = get_ms_dyn_sim(file_name, t_max, dist, p)
#     #inv = get_component(DynamicInverter, sim.sys, "generator-102-1")
#     #inv.inner_control.kpv = inv.inner_control.kpv*2 # increase gain 
#     #inv.inner_control.kpc = inv.inner_control.kpc*2 # increase gain
#     ss = small_signal_analysis(sim)
#     push!(stb, ss)
#     # Plot eigenvalues 
#     display(plot!(real(ss.eigenvalues), imag(ss.eigenvalues), seriestype=:scatter, label="N="*string(N)))
# end
# #xlims!(-100,100)
# xlabel!("Real")
# ylabel!("Imag")

plot()
for l in lrange;
    max_λ = []
    for N in Nrange;
        p = ExpParams(N, l, Z_c, r_km, x_km, g_km, b_km,abstol, reltol, maxiters)
        sim = get_ms_dyn_sim(file_name, t_max, dist, p)
        #inv = get_component(DynamicInverter, sim.sys, "generator-102-1")
        #inv.inner_control.kpv = inv.inner_control.kpv*2 # increase gain 
        #inv.inner_control.kpc = inv.inner_control.kpc*2 # increase gain
        ss = small_signal_analysis(sim)
        push!(stb, ss)
        push!(max_λ, maximum(real(ss.eigenvalues)))
        # # Plot eigenvalues 
        # display(plot!(real(ss.eigenvalues), imag(ss.eigenvalues), seriestype=:scatter, label="N="*string(N)))
    end
    display(plot!(Nrange, max_λ, seriestype=:scatter, label=string(l)*"km"))
end

# #xlims!(-100,100)
xlabel!("N segments")
ylabel!("Max real λ")

# Save fig
savefig("/Users/ruthkravis/Downloads/plot.png")

# View eigenvalue summary to work out which state participates most in each mode
eig_sum = summary_eigenvalues(stb[6])