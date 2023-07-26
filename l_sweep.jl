include("./all_sims.jl") # runs all_sims file but loads all functions defined in that file.
include("extra_functions.jl")

# Sweep line length for fixed number of line segments. 

file_name = "test_sys.json"
t_max = 30.0
dist = "CRC" # control reference change (on inverter)

Z_c = 380 # Ω
r_km = 0.05 # Ω/km
x_km = 0.488 # Ω/km
g_km = 0 # S/km
b_km = 3.371e-6 # S/km
N = 20

lrange = [10, 300, 800]; # Define line length range (km)
stb = []
max_λ = [] # for storing max real eigenvalue

plot()
for l in lrange;
    p = ExpParams(N, l, Z_c, r_km, x_km, g_km, b_km, 0.0, 0.0, 1e10)
    sim = get_ms_dyn_sim(file_name, t_max, dist, p)
    ss = small_signal_analysis(sim)
    push!(stb, ss)
    push!(max_λ, maximum(real(ss.eigenvalues)))
    display(plot!(real(ss.eigenvalues), imag(ss.eigenvalues), seriestype=:scatter, label="l="*string(l)*"km"))
end
xlims!(-1e-7,1e-7)

# View eig summary
n_eigs = length(stb[1].eigenvalues)
for i in stb;
    eig_sum = summary_eigenvalues(i)
    display(eig_sum[n_eigs-2,:])
end


