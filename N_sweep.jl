include("./all_sims.jl") # runs all_sims file but loads all functions defined in that file.
include("extra_functions.jl")
using PlotlyJS
file_name = "test_sys.json"
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

lrange = [50] # km 
Nrange = [1] # Define segment range 
stb = [] # small sig stability results
systems = []
sims = []

plotlyjs()
plot()

for l in lrange;
    max_λ = []
    second_λ = []
    for N in Nrange;
        p = ExpParams(N, l, Z_c, r_km, x_km, g_km, b_km,abstol, reltol, maxiters)
        sim = get_ms_dyn_sim(file_name, t_max, dist, p)
        sys = get_ms_dyn_sys(file_name, p)
        #inv = get_component(DynamicInverter, sim.sys, "generator-102-1")
        #inv.inner_control.kpv = inv.inner_control.kpv*2 # increase gain 
        #inv.inner_control.kpc = inv.inner_control.kpc*2 # increase gain
        ss = small_signal_analysis(sim)
        push!(stb, ss)
        push!(systems, sys)
        push!(sims, sim)
        push!(max_λ, maximum(real(ss.eigenvalues)))
        n_eigs = length(ss.eigenvalues)
        push!(second_λ, ss.eigenvalues[n_eigs-2])
        
        # # Plot eigenvalues 
        display(plot!(real(ss.eigenvalues), imag(ss.eigenvalues), seriestype=:scatter, label="N="*string(N)))
    end
    #display(plot!(Nrange, max_λ, seriestype=:scatter, label=string(l)*"km"))
    #display(plot!(Nrange, second_λ, seriestype=:scatter, label=string(l)*"km"))
end


xlims!(-4.254, -4.25)
xlabel!("Real")
ylabel!("Imag")
title!("1000 km line")
xlabel!("N segments")
ylabel!("Max real λ")
ylabel!("Second max real λ")

# Plot 

# View eigenvalue summary to work out which state participates most in each mode

eig_sum = summary_eigenvalues(stb[1])



savefig("/Users/ruthkravis/Downloads/plot.png")
