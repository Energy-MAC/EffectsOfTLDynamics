cd(@__DIR__)
using PowerSystems
using PowerSimulationsDynamics
using Plots
using Revise
using EffectsOfTLDynamics
using CSV
using PlotlyJS

const ETL = EffectsOfTLDynamics;
const PSY = PowerSystems;
const PSID = PowerSimulationsDynamics;

function draw_oval_with_edge(x, y, a, b, edgecolor, edgestrokealpha, edgestrokewidth)
    n_points = 100  # Number of points to approximate the oval shape
    θ = LinRange(0, 2π, n_points)  # Create an array of angles

    # Calculate the x and y coordinates of points on the oval
    x_points = a * cos.(θ) .+ x
    y_points = b * sin.(θ) .+ y

    # Plot the oval's edge
    plot!(x_points, y_points, linecolor=edgecolor, linealpha=edgestrokealpha, linewidth=edgestrokewidth, line=:path, primary=false)
end

function get_eig_set(n)
    kq = param_df[n,:kq]
    Ta = param_df[n,:Ta]
    kd = param_df[n,:kd]
    kpv = param_df[n,:kpv]
    kiv = param_df[n,:kiv]
    kpc = param_df[n,:kpc]
    kic = param_df[n,:kic]
    gfm_share = 1.0 - gfl_share;

    # Create system  
    sys, p1, p3, perturbation1, perturbation3 = create_sys_with_params(file_name, kpv, kiv, Ta, kd, kq, kpc, kic, load_scale, line_scale, inv_share, gfm_share, gfm_bus, gfl_bus);

    sys1 = deepcopy(sys)
    sys2 = deepcopy(sys)
    sys3 = deepcopy(sys)
    sys4 = deepcopy(sys)

    tspan = (0.0, 0.25)

    sys_alg = ETL.build_new_impedance_model!(sys1, p1, false, "")
    sim_alg = build_sim(sys_alg, tspan, perturbation1, false, p1);
    eigs_alg = small_signal_analysis(sim_alg).eigenvalues

    sys_dyn = ETL.build_new_impedance_model!(sys2, p1, true, "")
    sim_dyn = build_sim(sys_dyn, tspan, perturbation1, true, p1);
    eigs_dyn = small_signal_analysis(sim_dyn).eigenvalues

    sys_mssb = ETL.build_seg_model!(sys3, p1, true, "")
    sim_mssb = build_sim(sys_mssb, tspan, perturbation1, true, p1);
    eigs_mssb = small_signal_analysis(sim_mssb).eigenvalues
    
    return eigs_alg, eigs_dyn, eigs_mssb, sim_alg, sim_dyn
end

include("9bus_analysis_functions.jl")
file_name = "../data/json_data/9bus_VSM3_SM1_GFL2.json"

# Load a df with parameter choices 
param_sample_file_name = "Nov6_1000_samples_lseg_10/samples.csv"
param_df = DataFrame(CSV.File(param_sample_file_name))
load_scale = 0.4;
line_scale = 1.0;
inv_share = 0.7;
gfl_share = 0.6;

colorp = Plots.palette(:default); # color palette 
st = :scatter # plot style 
ms = 10; # marker size 
msw = 1; # marker size width 

eigs_alg_stable, eigs_dyn_stable, eigs_mssb_stable, sim_alg_stable, sim_dyn_stable = get_eig_set(24); #all stable 
eigs_alg_unstable, eigs_dyn_unstable, eigs_mssb_unstable, sim_alg_unstable, sim_dyn_unstable = get_eig_set(50); # alg unstable 

alpha = 0.5; # opacity 

# summary_eigenvalues(small_signal_analysis(sim_alg))[end-3:end, :]
# summary_eigenvalues(small_signal_analysis(sim_dyn))[end-3:end, :]

#plot(real(eigs_mssb), imag(eigs_mssb), st=st,label=L"\mathrm{MSSB}", color=colorp[3], ms=ms, msw=msw)
Plots.plot(size=(1000,800))
plot!(real(eigs_alg_stable), imag(eigs_alg_stable), seriestype=st, label=L"\mathrm{statpi}-k_q=0.05", color=colorp[1], ms=ms, msw=msw)
plot!(real(eigs_alg_unstable), imag(eigs_alg_unstable), seriestype=st, label=L"\mathrm{statpi}-k_q=0.2", color=colorp[1], ms=ms, msw=msw, alpha=alpha)
plot!(real(eigs_dyn_stable), imag(eigs_dyn_stable), st=st,label=L"\mathrm{dynpi}-k_q=0.05", color=colorp[2], ms=ms,msw=msw)
plot!(real(eigs_dyn_unstable), imag(eigs_dyn_unstable), st=st,label=L"\mathrm{dynpi}-k_q=0.2", color=colorp[2], ms=ms,msw=msw, alpha=alpha)
plot!(real(eigs_alg_stable), imag(eigs_alg_stable), seriestype=st, color=colorp[1], ms=ms, msw=msw, primary=false)

ylabel!(L"\Im \ (\lambda)")
xlabel!(L"\Re \ (\lambda)")

# FONT SIZING
labelfontsize=25;
ticksize = 16;
legendsize = 16;
plot!(xguidefontsize=labelfontsize, yguidefontsize=labelfontsize, tickfontsize=ticksize,legendfontsize=legendsize,legend=:bottomright)
plot!(framestyle=:box)

xlims!(-7.5,0.2)
ylims!(-40,40)

lens!([-0.03, 0.03], [-8, 8];
inset=(1, bbox(0.62, 0.03, 0.35, 0.35)), 
subplot=2, ticks=true, framestyle=:box,
lw=2, ls=:dot, lc=:orange)

exp_folder = "Nov6_1000_samples_lseg_10/"
savefig(exp_folder*"eig_plot.png")
savefig(exp_folder*"eig_plot.svg")







#plot!(real(eigs_mssb), imag(eigs_mssb), st=st,label=L"\mathrm{MSSB-2}", color=colorp[3], ms=ms, msw=msw, alpha=alpha)
plot!(real(eigs_dyn), imag(eigs_dyn), st=st,label=L"\mathrm{dynpi-2}", color=colorp[2], ms=ms,msw=msw, alpha=alpha)
plot!(real(eigs_alg), imag(eigs_alg), seriestype=st, label=L"\mathrm{statpi-2}", color=colorp[1], ms=ms, msw=msw, alpha=alpha)

xlims!(-10,0.2)
ylims!(-40,40)

summary_eigenvalues(small_signal_analysis(sim_alg))[end-5:end, :]
summary_eigenvalues(small_signal_analysis(sim_dyn))[end-5:end, :]


lens!([-0.03, 0.03], [-8, 8];
inset=(1, bbox(0.62, 0.03, 0.35, 0.35)), 
subplot=2, ticks=true, framestyle=:box,
lw=2, ls=:dot, lc=:orange)


# # # MAKE ANNOTATIONS
# draw_oval_with_edge(-3, 0, 8.0, 1e4, :red, 0.7, 2) 
# annotate!(-5, -1.2e4, (L"GFL/SM\ states", 16, :red, :top))
# exp_folder = "Nov5_final_choice_v2/"
# savefig(exp_folder*"eig_plot.png")
# #xlims!(-12,1)
# # # draw_oval_with_edge(-95, 3.8e4, 100.0, 3.3e4, :red, 0.7, 2) 
# # # annotate!(-85, 6.8e4, (L"Line\ states", 16, :red, :top))


# # plot2 = Plots.plot(real(eigs_alg), imag(eigs_alg), seriestype=st, label=L"\mathrm{statpi}", color=colorp[1], ms=ms, msw=msw, size=(1000,800))

# # plot!(real(eigs_dyn), imag(eigs_dyn), st=st,label=L"\mathrm{dynpi}", color=colorp[2], ms=ms,msw=msw)
# # plot!(real(eigs_mssb), imag(eigs_mssb), st=st,label=L"\mathrm{MSSB}", color=colorp[3], ms=ms, msw=msw)
# # plot!(real(eigs_msmb), imag(eigs_msmb), st=st,label=L"\mathrm{MSMB}", color=colorp[4], ms=ms, msw=msw)

# # labelfontsize=25;
# # ticksize = 16;
# # legendsize = 16;
# # plot!(xguidefontsize=labelfontsize, yguidefontsize=labelfontsize, tickfontsize=ticksize,legendfontsize=legendsize,legend=:left)
# # plot!(framestyle=:box)

# # xlims!(-500,50)


summary_eigenvalues(small_signal_analysis(sim_alg))[end-5:end, :]
summary_eigenvalues(small_signal_analysis(sim_dyn))[end-5:end, :]

alg_pf = summary_participation_factors(small_signal_analysis(sim_alg))[!,[:Name, :λ_47]]
dyn_pf = summary_participation_factors(small_signal_analysis(sim_dyn))[!,[:Name, :λ_70]]

sum(sum.(alg_pf[!,:λ_47]))
sum(sum.(dyn_pf[!,:λ_70]))

# filter(:λ_47 => x -> x > 0.001, alg_pf)
# filter(:λ_70 => x -> x > 0.001, dyn_pf)